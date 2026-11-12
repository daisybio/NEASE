"""Step 6: merge the shards of 05_pdb_contacts.py.

Writes to <work_dir>/intermediate/pdb/:
    ddi_pdb.parquet          domain-domain catalogue: Pfam pairs in contact (>= pdb.min_domain_contacts residue
                             pairs) in PDB structures of any species, with the number of entries
    interfaces_pdb.parquet   residue-level interfaces between human or mouse proteins (>= pdb.min_interface_contacts
                             residue pairs), residue pairs of all entries of a protein pair merged
    entries_log.parquet      status of every entry

usage: python 06_pdb_merge.py --config config.yaml
"""
import argparse

import pandas as pd

from common import load_config, log, paths, read_uniprot_entrez


def read_shards(folder, kind, shards):
    files = [folder / f'shard_{s:05d}.{kind}.parquet' for s in range(shards)]
    missing = [f.name for f in files if not f.exists()]
    if missing:
        raise SystemExit(f'{len(missing)} shards are missing (rerun these array tasks): {missing[:20]}')
    return pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    args = parser.parse_args()
    cfg = load_config(args.config)
    folder = paths(cfg)['intermediate'] / 'pdb'
    shards = cfg['pdb']['shards']

    logs = read_shards(folder / 'contacts', 'log', shards)
    logs.to_parquet(folder / 'entries_log.parquet')
    log('entries:', logs['status'].str.split(':').str[0].value_counts().to_dict())

    # domain-domain catalogue
    domains = read_shards(folder / 'contacts', 'domains', shards)
    domains = domains[domains['n_contacts'] >= cfg['pdb']['min_domain_contacts']]
    swap = domains['pfam_1'] > domains['pfam_2']
    domains.loc[swap, ['pfam_1', 'pfam_2', 'uniprot_1', 'uniprot_2']] = domains.loc[swap, ['pfam_2', 'pfam_1', 'uniprot_2', 'uniprot_1']].values
    ddi = domains.groupby(['pfam_1', 'pfam_2']).agg(n_entries=('pdb', 'nunique'), max_contacts=('n_contacts', 'max')).reset_index()
    ddi.to_parquet(folder / 'ddi_pdb.parquet')
    log(f'domain-domain catalogue: {len(ddi)} Pfam pairs ({ddi["pfam_1"].eq(ddi["pfam_2"]).sum()} homotypic) from '
        f'{domains["pdb"].nunique()} entries')

    # interfaces of human / mouse proteins
    interfaces = read_shards(folder / 'contacts', 'interfaces', shards)
    interfaces = interfaces[interfaces['n_pairs'] >= cfg['pdb']['min_interface_contacts']]
    species_of = {}
    for species in cfg['species']:
        species_of.update(dict.fromkeys(read_uniprot_entrez(cfg, species)['uniprot'], species))
    interfaces['species'] = interfaces['uniprot_1'].map(species_of)
    interfaces = interfaces[interfaces['species'].notna() & (interfaces['species'] == interfaces['uniprot_2'].map(species_of))]
    # one orientation per protein pair (uniprot_1 <= uniprot_2), then merge the residue pairs of all entries
    swap = interfaces['uniprot_1'] > interfaces['uniprot_2']
    interfaces.loc[swap, ['uniprot_1', 'uniprot_2', 'pos_1', 'pos_2']] = interfaces.loc[swap, ['uniprot_2', 'uniprot_1', 'pos_2', 'pos_1']].values
    rows = []
    for (species, u1, u2), group in interfaces.groupby(['species', 'uniprot_1', 'uniprot_2']):
        pairs = sorted({pair for p1, p2 in zip(group['pos_1'], group['pos_2']) for pair in zip(p1, p2)})
        entries = sorted(group['pdb'].unique())
        rows.append((species, u1, u2, [a for a, _ in pairs], [b for _, b in pairs], len(pairs), len(entries), ' '.join(entries[:20])))
    merged = pd.DataFrame(rows, columns=['species', 'uniprot_1', 'uniprot_2', 'pos_1', 'pos_2', 'n_pairs', 'n_entries', 'pdb'])
    merged.to_parquet(folder / 'interfaces_pdb.parquet')
    log('interfaces:', merged.groupby('species').size().to_dict(), 'protein pairs')
