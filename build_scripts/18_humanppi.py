"""Step 18 (optional): residue-level interfaces of the predicted human protein interactions of Zhang et al.

Zhang J, Humphreys IR, ..., Baker D, Cong Q. Predicting protein-protein interactions in the human proteome.
Science 2025 (doi:10.1126/science.adt1630); data: Dryad doi:10.5061/dryad.15dv41p84 (CC0).

The study screened ~190 million human protein pairs with RoseTTAFold2-PPI and AlphaFold2 on deep omics MSAs and
predicted 17,849 interactions at an expected precision of 90% (29,257 at 80%). Most predicted pairs (~80%) were
screened because they were already in BioGRID, STRING (physical) or UniProt; ~20% are in no interaction database.
The 3D models (best_models.tar.gz, 67 GB) are of protein segments, with residues numbered as in UniProt.

For every predicted pair of the chosen set (humanppi.precision; with humanppi.require_known only pairs listed in
an interaction database), the models with confident contacts (those that come with a .contacts file) are read:
residue pairs of the two chains with heavy atoms within alphafold.contact_distance and both residues with pLDDT
(from the .npz file) >= alphafold.min_plddt; >= alphafold.min_interface_contacts residue pairs. Residues are
placed on the UniProt sequence by their numbers, or by alignment when the numbering does not fit the current
sequence. The interfaces of all segment pairs of a protein pair are merged.

Evidence tiers: humanppi (residue level), humanppi_pair (domain level); mouse: humanppi_ortholog, humanppi_pair_ortholog.

The two Dryad files are given in the configuration (humanppi.predictions, humanppi.models); Dryad does not allow
scripted downloads without an API token. No internet access is needed. The result is small and can be computed on
any machine that has the model archive and copied to <work_dir>/intermediate/humanppi/ of the build.

Writes <work_dir>/intermediate/humanppi/interfaces.parquet:
    species, uniprot_1, uniprot_2, pos_1, pos_2, n_pairs, models, rf_prob, af_prob, source, databases, pdb_template

usage: python 18_humanppi.py --config config.yaml [--max-pairs 100]
"""
import argparse
import io
import tarfile
from multiprocessing import Pool
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

from common import chain_positions, load_config, log, paths, threads, uniprot_sequences
from contacts import drop_low_confidence, residue_contacts

BATCH = 200  # protein pairs read from the archive before they are processed (bounds the memory)


def settings(cfg):
    s = dict(enabled=False, precision=90, require_known=False, min_identity=90,
             predictions='downloads/humanppi/final_predictions.tar.gz', models='downloads/humanppi/best_models.tar.gz')
    s.update(cfg.get('humanppi') or {})
    for key in ('predictions', 'models'):
        path = Path(s[key]).expanduser()
        s[key] = path if path.is_absolute() else cfg['work_dir'] / path
    s.update(contact_distance=cfg['alphafold']['contact_distance'], min_plddt=cfg['alphafold']['min_plddt'],
             min_contacts=cfg['alphafold']['min_interface_contacts'])
    return s


def predictions(s):
    '''The predicted pairs of the chosen set: Protein1, Protein2 and their scores and provenance.'''
    name = f'final_predictions/final_predictions_{s["precision"]}.tsv'
    with tarfile.open(s['predictions']) as archive:
        table = pd.read_csv(archive.extractfile(name), sep='\t', comment='#', dtype=str)
    log(f'{len(table)} predicted interactions at {s["precision"]}% expected precision '
        f'({(table["allDBs"] == "none").sum()} in no interaction database)')
    if s['require_known']:
        table = table[table['allDBs'] != 'none']
        log(f'{len(table)} of them in an interaction database (humanppi.require_known)')
    return table.reset_index(drop=True)


def model_interface(pdb_text, npz_bytes, sequences, s):
    '''Residue pairs (UniProt positions of protein 1, protein 2) in contact in one model, or None.'''
    st = gemmi.read_pdb_string(pdb_text)
    st.remove_hydrogens()
    model = st[0]
    if len(model) != 2:
        return None
    plddt = np.load(io.BytesIO(npz_bytes), allow_pickle=False)['plddt']
    if len(plddt) != sum(len(chain) for chain in model):
        return None
    maps, k = [], 0
    for chain, sequence in zip(model, sequences):
        for i, residue in enumerate(chain, 1):
            residue.label_seq = i  # index along the chain for residue_contacts
            for atom in residue:
                atom.b_iso = float(plddt[k])
            k += 1
        positions = chain_positions(chain, sequence, s['min_identity'])
        if positions is None:
            return None
        maps.append(positions)
    names = [chain.name for chain in model]
    drop_low_confidence(model, s['min_plddt'])
    pairs = set()
    for (c1, c2), residues in residue_contacts(model, s['contact_distance']).items():
        for a, b in residues:
            if c1 == names[1]:  # orient as protein 1, protein 2
                a, b = b, a
            if a in maps[0] and b in maps[1]:
                pairs.add((maps[0][a], maps[1][b]))
    return pairs


def pair_job(job):
    '''All confident models of one protein pair -> merged interface.'''
    (p1, p2), files, sequences, s = job
    pairs, used = set(), []
    for stem, parts in sorted(files.items()):
        if 'contacts' not in parts or 'pdb' not in parts or 'npz' not in parts:
            continue  # no confident inter-protein contacts in this model
        found = model_interface(parts['pdb'].decode(), parts['npz'], sequences, s)
        if found and len(found) >= s['min_contacts']:
            pairs |= found
            used.append(stem)
    if not pairs:
        return None
    pairs = sorted(pairs)
    return p1, p2, [a for a, _ in pairs], [b for _, b in pairs], len(pairs), ' '.join(used)


def read_pairs(s, wanted, sequences, max_pairs=None):
    '''Stream the model archive and yield one job per wanted protein pair (all its files).'''
    current, files, n = None, {}, 0
    with tarfile.open(s['models'], 'r|gz') as archive:
        for member in archive:
            parts = member.name.split('/')
            if len(parts) != 3 or not member.isfile():
                continue
            folder, filename = parts[1], parts[2]
            if folder != current:
                if current is not None and files:
                    yield files_job(current, files, wanted, sequences, s)
                    n += 1
                    if max_pairs and n >= max_pairs:
                        return
                current, files = folder, {}
            key = tuple(folder.split('_', 1))
            if frozenset(key) not in wanted:
                continue
            stem, suffix = filename.rsplit('.', 1)
            files.setdefault(stem, {})[suffix] = archive.extractfile(member).read()
        if current is not None and files:
            yield files_job(current, files, wanted, sequences, s)


def files_job(folder, files, wanted, sequences, s):
    p1, p2 = folder.split('_', 1)
    return (p1, p2), files, (sequences.get(p1, ''), sequences.get(p2, '')), s


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    parser.add_argument('--max-pairs', type=int, help='only the first protein pairs of the archive (for testing)')
    args = parser.parse_args()
    cfg = load_config(args.config)
    s = settings(cfg)
    if not s['enabled']:
        raise SystemExit('The humanPPI predictions are not enabled (humanppi.enabled in the configuration)')
    for key in ('predictions', 'models'):
        if not s[key].exists():
            raise SystemExit(f'humanppi.{key}: {s[key]} not found (download it from Dryad, doi:10.5061/dryad.15dv41p84)')
    table = predictions(s)
    wanted = {frozenset(p) for p in zip(table['Protein1'], table['Protein2'])}
    sequences, _ = uniprot_sequences(cfg, set(table['Protein1']) | set(table['Protein2']))
    log(f'{len(sequences)} of {len(set(table["Protein1"]) | set(table["Protein2"]))} proteins with a current UniProt sequence')

    rows, n_jobs = [], 0
    with Pool(threads()) as pool:
        batch = []
        for job in read_pairs(s, wanted, sequences, args.max_pairs):
            if not all(job[2]):
                continue  # a protein without a current UniProt sequence (merged or deleted entry)
            batch.append(job)
            if len(batch) >= BATCH:
                rows += [r for r in pool.map(pair_job, batch) if r]
                n_jobs += len(batch)
                batch = []
                if n_jobs % 2000 == 0:
                    log(f'  {n_jobs} protein pairs read, {len(rows)} with an interface')
        rows += [r for r in pool.map(pair_job, batch) if r]
        n_jobs += len(batch)
    result = pd.DataFrame(rows, columns=['uniprot_1', 'uniprot_2', 'pos_1', 'pos_2', 'n_pairs', 'models'])
    # orient uniprot_1 <= uniprot_2 as the other interface tables, and add the scores and provenance
    swap = result['uniprot_1'] > result['uniprot_2']
    result.loc[swap, ['uniprot_1', 'uniprot_2', 'pos_1', 'pos_2']] = result.loc[swap, ['uniprot_2', 'uniprot_1', 'pos_2', 'pos_1']].values
    info = table.assign(key=[frozenset(p) for p in zip(table['Protein1'], table['Protein2'])])
    info = info.set_index('key')[['RFprob', 'AFprob', 'Source', 'allDBs', 'PDBtemp']]
    keys = [frozenset(p) for p in zip(result['uniprot_1'], result['uniprot_2'])]
    for column, name in [('RFprob', 'rf_prob'), ('AFprob', 'af_prob'), ('Source', 'source'), ('allDBs', 'databases'), ('PDBtemp', 'pdb_template')]:
        result[name] = [info.at[k, column] for k in keys]
    result.insert(0, 'species', 'Human')
    out = paths(cfg)['intermediate'] / 'humanppi'
    out.mkdir(exist_ok=True)
    result.to_parquet(out / 'interfaces.parquet')
    log(f'{n_jobs} predicted pairs with models read; {len(result)} with a confident interface '
        f'({(result["databases"] == "none").sum()} in no interaction database); written: {out / "interfaces.parquet"}')
