"""Step 19 (optional): residue-level interfaces of Boltz-2 models of measured direct human interactions.

Ille AM, Markosian C, Burley SK, Pasqualini R, Arap W. Human protein interactome structure prediction at scale with
Boltz-2. bioRxiv 2025 (doi:10.1101/2025.07.03.663068); models: github.com/structural-interactome/human-interactome
(MIT licence).

The study modelled 1,394 binary human interactions from IntAct that were measured as direct binding (pull-down
with purified proteins), both proteins <= 1,000 residues, with Boltz-2. Some models are of isoforms (P05067-4) or
processed chains (P05067-PRO_0000000092). The models have per-residue pLDDT (B-factor column); interface scores
are not in the files.

Interfaces as for the AlphaFold models: residue pairs of the two chains with heavy atoms within
alphafold.contact_distance, both residues with pLDDT >= alphafold.min_plddt, >= alphafold.min_interface_contacts
residue pairs. A chain is placed on the UniProt entry it was modelled from if that entry can be placed on exons
(step 2; isoform entries identical to an Ensembl translation included), otherwise on the canonical entry, by
residue numbers or by alignment (>= boltz.min_identity % identity).

Evidence tiers: boltz (residue level), boltz_pair (domain level); mouse: boltz_ortholog, boltz_pair_ortholog.

Writes <work_dir>/intermediate/boltz/interfaces.parquet:
    species, uniprot_1, uniprot_2, pos_1, pos_2, n_pairs, model, mean_plddt

usage: python 19_boltz.py --config config.yaml
"""
import argparse
import re
import zipfile
from multiprocessing import Pool
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

from common import chain_positions, input_path, load_config, log, paths, placeable_proteins, threads, uniprot_sequences
from contacts import drop_low_confidence, residue_contacts

NAME = re.compile(r'^([A-Z0-9]+(?:-\d+)?(?:-PRO_\d+)?)_([A-Z0-9]+(?:-\d+)?(?:-PRO_\d+)?)_Boltz')


def settings(cfg):
    s = dict(enabled=False, models=None, min_identity=90)
    s.update(cfg.get('boltz') or {})
    path = Path(s['models']).expanduser() if s['models'] else input_path(cfg, 'boltz')
    s['models'] = path if path.is_absolute() else cfg['work_dir'] / path
    s.update(contact_distance=cfg['alphafold']['contact_distance'], min_plddt=cfg['alphafold']['min_plddt'],
             min_contacts=cfg['alphafold']['min_interface_contacts'])
    return s


def model_files(path):
    '''(file name, text) of every model: a folder of .cif files or the repository zip.'''
    if path.is_dir():
        for f in sorted(path.rglob('*.cif')):
            yield f.name, f.read_text()
    else:
        with zipfile.ZipFile(path) as archive:
            for name in sorted(n for n in archive.namelist() if n.endswith('.cif') and '/interaction_structures/' in n):
                yield name.rsplit('/', 1)[-1], archive.read(name).decode()


def read_model(text):
    '''The files repeat the _ma_qa_metric table, which the mmCIF parser rejects: the second copy is dropped.'''
    first = text.find('loop_\n_ma_qa_metric.id')
    second = text.find('loop_\n_ma_qa_metric.id', first + 1)
    if second > 0:
        end = text.find('\n#', second)
        text = text[:second] + text[end + 1:]
    st = gemmi.make_structure_from_block(gemmi.cif.read_string(text).sole_block())
    st.remove_hydrogens()
    return st


def target(entry, placeable, canonical):
    '''UniProt entry (and its sequence) a modelled protein is placed on: itself, or its canonical entry.'''
    if entry in placeable:
        return entry, placeable[entry]
    base = entry.split('-')[0]  # isoform or processed chain -> canonical accession
    sequence = placeable.get(base) or canonical.get(base)
    return (base, sequence) if sequence else (None, None)


def model_interface(job):
    name, text, entries, sequences, s = job
    try:
        st = read_model(text)
    except (RuntimeError, ValueError):
        return name, None, 'unreadable'
    model = st[0]
    if len(model) != 2:
        return name, None, 'not two chains'
    maps = []
    for chain, sequence in zip(model, sequences):
        for i, residue in enumerate(chain, 1):
            residue.label_seq = i
        positions = chain_positions(chain, sequence, s['min_identity'])
        if positions is None:
            return name, None, 'chain not placed on UniProt'
        maps.append(positions)
    plddt = [r[0].b_iso for chain in model for r in chain if len(r)]
    names = [chain.name for chain in model]
    drop_low_confidence(model, s['min_plddt'])
    pairs = set()
    for (c1, c2), residues in residue_contacts(model, s['contact_distance']).items():
        for a, b in residues:
            if c1 == names[1]:
                a, b = b, a
            if a in maps[0] and b in maps[1]:
                pairs.add((maps[0][a], maps[1][b]))
    if len(pairs) < s['min_contacts']:
        return name, None, 'no confident interface'
    pairs = sorted(pairs)
    return name, (entries[0], entries[1], [a for a, _ in pairs], [b for _, b in pairs], len(pairs), round(float(np.mean(plddt)), 1)), 'ok'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    args = parser.parse_args()
    cfg = load_config(args.config)
    s = settings(cfg)
    if not s['enabled']:
        raise SystemExit('The Boltz-2 models are not enabled (boltz.enabled in the configuration)')
    if not s['models'].exists():
        raise SystemExit(f'boltz.models: {s["models"]} not found (python 01_download.py --config ... --only boltz)')
    proteins = placeable_proteins(cfg, 'Human').drop_duplicates('uniprot')
    placeable = dict(zip(proteins['uniprot'], proteins['sequence']))
    files = list(model_files(s['models']))
    entries = {}
    for name, _ in files:
        m = NAME.match(name)
        if m:
            entries[name] = (m.group(1), m.group(2))
    canonical, _ = uniprot_sequences(cfg, {e.split('-')[0] for pair in entries.values() for e in pair})
    jobs = []
    for name, text in files:
        if name not in entries:
            continue
        (t1, q1), (t2, q2) = (target(e, placeable, canonical) for e in entries[name])
        if q1 and q2:
            jobs.append((name, text, (t1, t2), (q1, q2), s))
    log(f'{len(files)} model files, {len(jobs)} with both proteins placeable on UniProt sequences')
    with Pool(threads()) as pool:
        results = pool.map(model_interface, jobs, chunksize=4)
    log('models:', pd.Series([status for _, _, status in results]).value_counts().to_dict())
    rows = [(*row, name) for name, row, status in results if row]
    table = pd.DataFrame(rows, columns=['uniprot_1', 'uniprot_2', 'pos_1', 'pos_2', 'n_pairs', 'mean_plddt', 'model'])
    swap = table['uniprot_1'] > table['uniprot_2']
    table.loc[swap, ['uniprot_1', 'uniprot_2', 'pos_1', 'pos_2']] = table.loc[swap, ['uniprot_2', 'uniprot_1', 'pos_2', 'pos_1']].values
    table.insert(0, 'species', 'Human')
    out = paths(cfg)['intermediate'] / 'boltz'
    out.mkdir(exist_ok=True)
    table.to_parquet(out / 'interfaces.parquet')
    log(f'{len(table)} interfaces ({table["uniprot_1"].str.contains("-").sum() + table["uniprot_2"].str.contains("-").sum()} '
        f'sides on isoform entries); written: {out / "interfaces.parquet"}')
