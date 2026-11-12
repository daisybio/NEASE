"""Step 15 (optional; the download needs internet access): residue-level interfaces from Interactome3D.

Interactome3D (interactome3d.irbbarcelona.org) collects for the protein interactions of an organism experimental
structures (also of closely related proteins) and builds homology models: of the complex on a related structure
('Model') or of two interacting domains on a domain-domain template ('Dom_dom_model'). This step downloads the
interaction tables and structure tarballs of the configured set, maps the residues of every chain onto the UniProt
sequence (local alignment; the files are numbered from 1, not by UniProt position) and computes residue contacts
as for the PDB (heavy atoms within pdb.contact_distance, >= pdb.min_interface_contacts residue pairs).

Evidence tiers (network.evidence / network.interface_evidence):
    i3d_structure, i3d_model                 residue-level interfaces of experimental structures / models
    i3d_structure_pair, i3d_model_pair       domain pairs in contact in them (joint graph)
    *_ortholog (mouse)                       the human ones, transferred through orthologs

Writes <work_dir>/intermediate/interactome3d/interfaces.parquet:
    species, kind (i3d_structure / i3d_model), type, file, uniprot_1, uniprot_2, pos_1, pos_2, n_pairs, identity_1, identity_2

Run after step 2 and before step 13. Downloads go to <work_dir>/downloads/interactome3d/ and are skipped when present,
so they can be fetched first on a node with internet access:
usage: python 15_interactome3d.py --config config.yaml [--download-only] [--threads 8]
"""
import argparse
import re
import tarfile
from multiprocessing import Pool

import gemmi
import pandas as pd
import requests

from common import align_positions, download, load_config, log, paths, placeable_proteins
from contacts import residue_contacts

ORGANISM = {'Human': 'human', 'Mouse': 'mouse'}
PAGE_URL = 'https://interactome3d.irbbarcelona.org/downloadset.php?queryid={organism}&release=current&path={set}'
FILE_URL = 'https://interactome3d.irbbarcelona.org/user_data/{organism}/download/{set}/{file}'
COLUMNS = ['species', 'kind', 'type', 'file', 'uniprot_1', 'uniprot_2', 'pos_1', 'pos_2', 'n_pairs', 'identity_1', 'identity_2']


def settings(cfg):
    s = dict(enabled=False, species=['Human', 'Mouse'], set='representative', types=['Structure', 'Model', 'Dom_dom_model'],
             min_identity=50)
    s.update(cfg.get('interactome3d') or {})
    s['species'] = [sp for sp in s['species'] if sp in cfg['species']]
    return s


def data_folder(cfg, species, s):
    return paths(cfg)['downloads'] / 'interactome3d' / f'{ORGANISM[species]}_{s["set"]}'


def fetch(cfg, species, s):
    '''Download interactions.dat and the interactions_*.tgz tarballs listed on the download page of the set.'''
    fields = dict(organism=ORGANISM[species], set=s['set'])
    page = requests.get(cfg['sources'].get('i3d_page', PAGE_URL).format(**fields), timeout=300)
    page.raise_for_status()
    files = sorted(set(re.findall(rf'user_data/{fields["organism"]}/download/{s["set"]}/(interactions[\w.]*\.(?:dat|tgz))', page.text)))
    if 'interactions.dat' not in files:
        raise SystemExit(f'No Interactome3D files found on {page.url}')
    folder = data_folder(cfg, species, s)
    log(f'{species}: {len(files)} files ({folder})')
    for name in files:
        download(cfg['sources'].get('i3d_file', FILE_URL).format(file=name, **fields), folder / name, timeout=300)
    (folder / 'files.txt').write_text('\n'.join(files) + '\n')


def downloaded(cfg, species, s):
    '''True if all files of the set are there (the list is written when the download is complete).'''
    listing = data_folder(cfg, species, s) / 'files.txt'
    return listing.exists() and all((listing.parent / f).exists() for f in listing.read_text().split())


def one_letter(name):
    info = gemmi.find_tabulated_residue(name)
    return info.one_letter_code.upper() if info and info.one_letter_code.strip() else 'X'


def residue_map(residues, sequence, min_identity):
    '''Positions in the UniProt sequence of the residues of a chain (local alignment), and the identity in %.'''
    chain_seq = ''.join(one_letter(r.name) for r in residues)
    chain_pos, uniprot_pos, identical = align_positions(chain_seq, sequence)
    mapping = dict(zip(chain_pos, uniprot_pos))
    identity = 100 * identical / max(len(mapping), 1)
    if len(mapping) < 10 or identity < min_identity:
        return None, round(identity, 1)
    return mapping, round(identity, 1)


def file_interface(text, uniprot_1, uniprot_2, sequences, distance, min_identity):
    '''Residue pairs (UniProt positions) between chain A (uniprot_1) and chain B (uniprot_2) of an Interactome3D file.'''
    st = gemmi.read_pdb_string(text)
    st.remove_hydrogens()
    st.remove_ligands_and_waters()
    model = st[0]
    chains = {chain.name: chain for chain in model}
    if 'A' not in chains or 'B' not in chains:
        return None
    maps, identities = [], []
    for name, uniprot in (('A', uniprot_1), ('B', uniprot_2)):
        chain = chains[name]
        for i, residue in enumerate(chain, 1):  # number the residues along the chain for residue_contacts
            residue.label_seq = i
        mapping, identity = residue_map(chain, sequences[uniprot], min_identity)
        if mapping is None:
            return None
        maps.append(mapping)
        identities.append(identity)
    pairs = sorted({(maps[0][a], maps[1][b]) for a, b in residue_contacts(model, distance).get(('A', 'B'), ())
                    if a in maps[0] and b in maps[1]})
    return pairs, identities


def process_tarball(job):
    path, wanted, sequences, distance, min_identity, min_pairs = job
    rows, failed = [], 0
    with tarfile.open(path) as archive:
        for member in archive:
            name = member.name.split('/')[-1]
            if name not in wanted:
                continue
            kind, i3d_type, u1, u2 = wanted[name]
            try:
                result = file_interface(archive.extractfile(member).read().decode(), u1, u2, sequences, distance, min_identity)
            except Exception:  # damaged or unusual files: count and continue
                failed += 1
                continue
            if result is None or len(result[0]) < min_pairs:
                continue
            pairs, (i1, i2) = result
            if u1 > u2:
                u1, u2, i1, i2, pairs = u2, u1, i2, i1, sorted((b, a) for a, b in pairs)
            rows.append((kind, i3d_type, name, u1, u2, [a for a, _ in pairs], [b for _, b in pairs], len(pairs), i1, i2))
    return rows, failed


def interfaces(cfg, species, s, threads):
    folder = data_folder(cfg, species, s)
    table = pd.read_csv(folder / 'interactions.dat', sep='\t', dtype=str)
    table = table[table['TYPE'].isin(s['types'])]
    # UniProt sequences of the entries that can be placed on exons (identical to a translation or aligned to one)
    proteins = placeable_proteins(cfg, species)[['uniprot', 'sequence']]
    sequences = dict(proteins.drop_duplicates('uniprot').itertuples(index=False))
    known = table['PROT1'].isin(sequences) & table['PROT2'].isin(sequences)
    log(f'{species}: {len(table)} structures/models {table["TYPE"].value_counts().to_dict()}; {known.sum()} with both proteins '
        'placeable on Ensembl exons (step 2)')
    table = table[known]
    wanted = {f: ('i3d_structure' if t == 'Structure' else 'i3d_model', t, p1, p2)
              for f, t, p1, p2 in table[['FILENAME', 'TYPE', 'PROT1', 'PROT2']].itertuples(index=False)}
    needed = set(table['PROT1']) | set(table['PROT2'])
    sequences = {u: q for u, q in sequences.items() if u in needed}
    tarballs = sorted(folder.glob('interactions_*.tgz'))
    jobs = [(path, wanted, sequences, cfg['pdb']['contact_distance'], s['min_identity'], cfg['pdb']['min_interface_contacts'])
            for path in tarballs]
    rows, failed = [], 0
    with Pool(threads) as pool:
        for part, n_failed in pool.imap_unordered(process_tarball, jobs):
            rows += part
            failed += n_failed
    result = pd.DataFrame(rows, columns=COLUMNS[1:]).assign(species=species)[COLUMNS]
    log(f'{species}: {len(result)} interfaces {result["type"].value_counts().to_dict()}, {failed} files failed')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    parser.add_argument('--download-only', action='store_true', help='only download the Interactome3D files')
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    cfg = load_config(args.config)
    s = settings(cfg)
    if not s['enabled']:
        raise SystemExit('Interactome3D is not enabled (interactome3d.enabled in the configuration)')
    for species in s['species']:
        if args.download_only or not downloaded(cfg, species, s):
            fetch(cfg, species, s)
    if not args.download_only:
        out = paths(cfg)['intermediate'] / 'interactome3d'
        out.mkdir(exist_ok=True)
        table = pd.concat([interfaces(cfg, species, s, args.threads) for species in s['species']], ignore_index=True)
        table.to_parquet(out / 'interfaces.parquet')
        log(f'written: {out / "interfaces.parquet"}')
