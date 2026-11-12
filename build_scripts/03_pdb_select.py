"""Step 3: select the PDB entries to analyse and split them into shards (array tasks of 05_pdb_contacts.py).

Entries with a Pfam domain on a protein chain (for the domain-domain catalogue, any species) and entries with
a chain of a human/mouse protein (for residue-level interfaces) are selected from SIFTS.

Writes to <work_dir>/intermediate/pdb/:
    entries.tsv            pdb, shard
    pfam_map.parquet       Pfam domains of the chains (SIFTS; label_seq_id ranges)
    uniprot_map.parquet    chain segments -> UniProt positions (SIFTS; label_seq_id ranges)

usage: python 03_pdb_select.py --config config.yaml [--max-entries 500]
"""
import argparse

import pandas as pd

from common import input_path, load_config, log, paths, read_uniprot_entrez


def read_sifts(path, columns):
    table = pd.read_csv(path, sep='\t', comment='#', usecols=list(columns), dtype=str).rename(columns=columns)
    for column in table.columns:
        if column.endswith(('_start', '_end', '_beg')):
            table[column] = pd.to_numeric(table[column], errors='coerce')
    return table.dropna().astype({c: 'int64' for c in table.columns if c.endswith(('_start', '_end', '_beg'))})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    parser.add_argument('--max-entries', type=int, help='random subset of entries (for testing)')
    args = parser.parse_args()
    cfg = load_config(args.config)
    out = paths(cfg)['intermediate'] / 'pdb'
    out.mkdir(exist_ok=True)

    pfam = read_sifts(input_path(cfg, 'sifts_pfam'), {'PDB': 'pdb', 'CHAIN': 'chain', 'PDB_START': 'label_start',
                                                      'PDB_END': 'label_end', 'PFAM_ACCESSION': 'pfam', 'UNIPROT_ACCESSION': 'uniprot'})
    uniprot = read_sifts(input_path(cfg, 'sifts_uniprot'), {'PDB': 'pdb', 'CHAIN': 'chain', 'SP_PRIMARY': 'uniprot',
                                                            'RES_BEG': 'label_beg', 'RES_END': 'label_end', 'SP_BEG': 'sp_beg'})
    log(f'SIFTS: {pfam["pdb"].nunique()} entries with Pfam domains, {uniprot["pdb"].nunique()} entries mapped to UniProt')

    species_proteins = set()
    for species in cfg['species']:
        species_proteins |= set(read_uniprot_entrez(cfg, species)['uniprot'])
    with_species = set(uniprot.loc[uniprot['uniprot'].isin(species_proteins), 'pdb'])
    entries = sorted(set(pfam['pdb']) | with_species)
    log(f'{len(entries)} entries ({len(with_species)} with {"/".join(cfg["species"])} proteins)')
    if args.max_entries:
        entries = sorted(pd.Series(entries).sample(min(args.max_entries, len(entries)), random_state=0))
        log(f'test subset: {len(entries)} entries')

    shards = cfg['pdb']['shards']
    pd.DataFrame({'pdb': entries, 'shard': [i % shards for i in range(len(entries))]}).to_csv(out / 'entries.tsv', sep='\t', index=False)
    keep = set(entries)
    pfam[pfam['pdb'].isin(keep)].reset_index(drop=True).to_parquet(out / 'pfam_map.parquet')
    uniprot[uniprot['pdb'].isin(keep)].reset_index(drop=True).to_parquet(out / 'uniprot_map.parquet')
    log(f'{shards} shards -> {out}')
