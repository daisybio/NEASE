"""Step 11: physical protein-protein interactions (NCBI gene IDs) from BioGRID, IntAct, HuRI, BioPlex and STRING.

Only physical interactions are used (see config ppi). Mouse: mouse interactions plus the human interactions of
genes with one-to-one orthologs (ppi.mouse_from_human).

Writes <work_dir>/intermediate/<species>/ppi.parquet:
    entrez_1, entrez_2, sources, n_sources, n_publications, human_sources (mouse: transferred from human)

usage: python 11_ppi.py --config config.yaml
"""
import argparse
import zipfile

import pandas as pd

from common import SPECIES, input_path, load_config, log, read_uniprot_entrez, species_dir

COLUMNS = ['entrez_1', 'entrez_2', 'source', 'publication']


def biogrid(cfg, taxa):
    keep = ['Entrez Gene Interactor A', 'Entrez Gene Interactor B', 'Experimental System', 'Experimental System Type',
            'Publication Source', 'Organism ID Interactor A', 'Organism ID Interactor B']
    with zipfile.ZipFile(input_path(cfg, 'biogrid')) as archive:
        name = [n for n in archive.namelist() if n.endswith('.txt')][0]
        parts = []
        with archive.open(name) as f:
            for chunk in pd.read_csv(f, sep='\t', usecols=keep, dtype=str, chunksize=500000):
                chunk = chunk[(chunk['Experimental System Type'] == 'physical')
                              & ~chunk['Experimental System'].isin(cfg['ppi']['biogrid_exclude_systems'])
                              & (chunk['Organism ID Interactor A'] == chunk['Organism ID Interactor B'])
                              & chunk['Organism ID Interactor A'].isin(taxa)]
                parts.append(pd.DataFrame({'taxid': chunk['Organism ID Interactor A'], 'entrez_1': chunk['Entrez Gene Interactor A'],
                                           'entrez_2': chunk['Entrez Gene Interactor B'], 'source': 'biogrid',
                                           'publication': chunk['Publication Source']}))
    table = pd.concat(parts, ignore_index=True)
    return {taxid: t.drop(columns='taxid') for taxid, t in table.groupby('taxid')}


def intact(cfg, taxa, uniprot_entrez):
    names = {0: 'id_1', 1: 'id_2', 8: 'publications', 9: 'taxid_1', 10: 'taxid_2', 11: 'type', 15: 'expansion', 35: 'negative'}
    with zipfile.ZipFile(input_path(cfg, 'intact')) as archive:
        name = [n for n in archive.namelist() if n.endswith('intact.txt')][0]
        parts = []
        with archive.open(name) as f:
            for chunk in pd.read_csv(f, sep='\t', usecols=list(names), header=0, dtype=str, chunksize=500000, quoting=3):
                chunk.columns = [names[i] for i in sorted(names)]
                for side in ('1', '2'):
                    chunk[f'taxid_{side}'] = chunk[f'taxid_{side}'].str.extract(r'taxid:(-?\d+)', expand=False)
                    chunk[f'uniprot_{side}'] = chunk[f'id_{side}'].str.extract(r'^uniprotkb:([A-Z0-9]+)', expand=False)
                chunk['type'] = chunk['type'].str.extract(r'(MI:\d{4})', expand=False)
                keep = ((chunk['negative'] != 'true') & chunk['type'].isin(cfg['ppi']['intact_types'])
                        & (chunk['taxid_1'] == chunk['taxid_2']) & chunk['taxid_1'].isin(taxa)
                        & chunk['uniprot_1'].notna() & chunk['uniprot_2'].notna())
                if cfg['ppi']['intact_exclude_expanded']:
                    keep &= chunk['expansion'].isin(['-', None])
                chunk = chunk[keep]
                parts.append(pd.DataFrame({'taxid': chunk['taxid_1'], 'uniprot_1': chunk['uniprot_1'], 'uniprot_2': chunk['uniprot_2'],
                                           'publication': chunk['publications'].str.extract(r'pubmed:(\d+)', expand=False)}))
    table = pd.concat(parts, ignore_index=True)
    out = {}
    for taxid, t in table.groupby('taxid'):
        mapping = uniprot_entrez[taxid]
        t = t.merge(mapping.rename(columns={'uniprot': 'uniprot_1', 'entrez': 'entrez_1'}), on='uniprot_1')
        t = t.merge(mapping.rename(columns={'uniprot': 'uniprot_2', 'entrez': 'entrez_2'}), on='uniprot_2')
        out[taxid] = t.assign(source='intact', publication='PUBMED:' + t['publication'])[COLUMNS]
    return out


def huri(cfg):
    path = input_path(cfg, 'huri')
    if not path.exists():
        log('  HuRI not downloaded, skipped (its interactions are also part of IntAct)')
        return pd.DataFrame(columns=COLUMNS)
    pairs = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1], names=['gene_1', 'gene_2'], dtype=str)
    genes = pd.read_parquet(species_dir(cfg, 'Human') / 'genes.parquet', columns=['gene_id', 'entrez']).dropna()
    for side in ('1', '2'):
        pairs = pairs.merge(genes.rename(columns={'gene_id': f'gene_{side}', 'entrez': f'entrez_{side}'}), on=f'gene_{side}')
    return pairs.assign(source='huri', publication=None)[COLUMNS]


def bioplex(cfg):
    parts = []
    for name in ('bioplex_293t', 'bioplex_hct116'):
        table = pd.read_csv(input_path(cfg, name), sep='\t', dtype={'GeneA': str, 'GeneB': str})
        table = table[table['pInt'] >= cfg['ppi']['bioplex_min_pint']]
        parts.append(pd.DataFrame({'entrez_1': table['GeneA'], 'entrez_2': table['GeneB'], 'source': 'bioplex', 'publication': None}))
    return pd.concat(parts, ignore_index=True)


def string(cfg, species):
    links = pd.read_csv(input_path(cfg, 'string_physical', species), sep=' ', usecols=['protein1', 'protein2', 'experimental'])
    links = links[links['experimental'] >= cfg['ppi']['string_min_experimental']]
    aliases = pd.read_csv(input_path(cfg, 'string_aliases', species), sep='\t', dtype=str)
    aliases.columns = ['protein', 'alias', 'source']
    aliases = aliases[aliases['source'].str.contains('entrez', case=False) & aliases['alias'].str.isdigit()]
    mapping = aliases[['protein', 'alias']].drop_duplicates()
    for side in ('1', '2'):
        links = links.merge(mapping.rename(columns={'protein': f'protein{side}', 'alias': f'entrez_{side}'}), on=f'protein{side}')
    return links.assign(source='string', publication=None)[COLUMNS]


def combine(evidence):
    table = pd.concat(evidence, ignore_index=True).dropna(subset=['entrez_1', 'entrez_2'])
    table = table[(table['entrez_1'] != '-') & (table['entrez_2'] != '-')]
    swap = table['entrez_1'] > table['entrez_2']
    table.loc[swap, ['entrez_1', 'entrez_2']] = table.loc[swap, ['entrez_2', 'entrez_1']].values
    grouped = table.groupby(['entrez_1', 'entrez_2'])
    return pd.DataFrame({'sources': grouped['source'].agg(lambda s: ','.join(sorted(set(s)))),
                         'n_sources': grouped['source'].nunique(),
                         'n_publications': grouped['publication'].nunique()}).reset_index()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    args = parser.parse_args()
    cfg = load_config(args.config)
    sources = cfg['ppi']['sources']
    taxa = {str(SPECIES[s]['taxid']): s for s in cfg['species']}

    evidence = {s: [] for s in cfg['species']}
    if 'biogrid' in sources:
        log('BioGRID')
        for taxid, table in biogrid(cfg, list(taxa)).items():
            evidence[taxa[taxid]].append(table)
    if 'intact' in sources:
        log('IntAct')
        uniprot_entrez = {str(SPECIES[s]['taxid']): read_uniprot_entrez(cfg, s) for s in cfg['species']}
        for taxid, table in intact(cfg, list(taxa), uniprot_entrez).items():
            evidence[taxa[taxid]].append(table)
    if 'huri' in sources and 'Human' in evidence:
        log('HuRI')
        evidence['Human'].append(huri(cfg))
    if 'bioplex' in sources and 'Human' in evidence:
        log('BioPlex')
        evidence['Human'].append(bioplex(cfg))
    if 'string' in sources:
        for species in cfg['species']:
            log(f'STRING {species}')
            evidence[species].append(string(cfg, species))

    networks = {}
    for species in cfg['species']:
        ppi = combine(evidence[species])
        ppi = ppi[(ppi['n_sources'] >= cfg['ppi']['min_sources']) | (ppi['n_publications'] >= cfg['ppi']['min_publications'])]
        networks[species] = ppi.assign(human_sources='')
        log(f'{species}: {len(ppi)} interactions,', ppi['sources'].str.split(',').explode().value_counts().to_dict())

    if 'Mouse' in networks and cfg['ppi']['mouse_from_human']:
        orthologs = pd.read_parquet(species_dir(cfg, 'Mouse').parent / 'orthologs.parquet', columns=['human_entrez', 'mouse_entrez']).dropna()
        human = networks['Human'][['entrez_1', 'entrez_2', 'sources']]
        for side in ('1', '2'):
            human = human.merge(orthologs.rename(columns={'human_entrez': f'entrez_{side}', 'mouse_entrez': f'mouse_{side}'}), on=f'entrez_{side}')
        human = pd.DataFrame({'entrez_1': human[['mouse_1', 'mouse_2']].min(axis=1), 'entrez_2': human[['mouse_1', 'mouse_2']].max(axis=1),
                              'human_sources': human['sources']})
        human = human.groupby(['entrez_1', 'entrez_2'])['human_sources'].agg(lambda s: ','.join(sorted(set(','.join(s).split(','))))).reset_index()
        mouse = networks['Mouse'].drop(columns='human_sources').merge(human, on=['entrez_1', 'entrez_2'], how='outer')
        mouse['sources'] = mouse['sources'].fillna('')
        mouse['human_sources'] = mouse['human_sources'].fillna('')
        mouse[['n_sources', 'n_publications']] = mouse[['n_sources', 'n_publications']].fillna(0).astype(int)
        networks['Mouse'] = mouse
        log(f'Mouse with human interactions: {len(mouse)} interactions ({(mouse["sources"] == "").sum()} only from human)')

    for species, ppi in networks.items():
        ppi.reset_index(drop=True).to_parquet(species_dir(cfg, species) / 'ppi.parquet')
