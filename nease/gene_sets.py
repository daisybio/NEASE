"""Gene sets for NEASE enrichment from MSigDB collections and KEGG (events.enrich_msigdb()).

The MSigDB collections are downloaded from MSigDB and the KEGG pathways from the KEGG REST API
the first time they are used, and kept in a local cache (see cache_dir()).

NEASE's edge enrichment test needs, for every gene set, the total degree of its
genes in the network used to find affected edges, and the number of edges of
that network. The built-in pathway tables ship these numbers precomputed; for
other gene sets (and the built-in mouse table) they are computed here from the
bundled networks. Both networks are gene-level (nodes are NCBI gene IDs):

    - PPI/DDI: protein interactions explained by a domain-domain interaction.
      Built from the joint graph (data/network/graph_*.pkl), whose nodes are
      domains ('<NCBI gene ID>/<Pfam ID>'), by merging the domains of each gene.
      Two domains of the same protein interacting give a self-interaction.
    - structural PPI: PPI/DDI plus
        * domain-motif interactions from ELM (database/ELM_interactions*): motif gene - domain gene
        * co-resolved interactions from PDB (database/pdb*): gene with the interface
          residue (mapped from its Ensembl ID, or symbol, to its NCBI gene ID) - co-resolved partner

The degree of a gene set is the sum of the degrees of its genes; as in the
built-in tables, a self-interaction adds 2 to the degree of its gene.
"""
import datetime
import json
import os
import shutil
import urllib.request
from functools import lru_cache

import networkx as nx
import pandas as pd

from .load import network, elm_interactions, pdb, database_mapping


here = os.path.dirname(__file__)

# Column names used for the computed degrees (same as the built-in human table)
DEGREE_STRUCTURAL = 'Degree in the structural PPI'
DEGREE_DDI = 'Degree in the PPI/DDI'

# MSigDB release used by NEASE, and its symbol tables bundled with NEASE (NCBI gene ID -> symbol)
MSIGDB_RELEASE = '2026.2'
MSIGDB_URL = 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/{release}.{species}/{file}.v{release}.{species}.json'
MSIGDB_SPECIES = {'Human': 'Hs', 'Mouse': 'Mm'}
SYMBOL_TABLES = {'Human': os.path.join(here, 'data/gene_sets/msigdb_ncbi_symbols_human.tsv'),
                 'Mouse': os.path.join(here, 'data/gene_sets/msigdb_ncbi_symbols_mouse.tsv')}

# Gene set databases of enrich_msigdb(): name -> (MSigDB collection, MSigDB file of the collection)
MSIGDB_DATABASES = {
    'Human': {'Hallmark': ('H', 'h.all'),
              'Positional': ('C1', 'c1.all'),
              'CGP': ('C2:CGP', 'c2.cgp'),
              'BioCarta': ('C2:CP:BIOCARTA', 'c2.cp.biocarta'),
              'KEGG_MEDICUS': ('C2:CP:KEGG_MEDICUS', 'c2.cp.kegg_medicus'),
              'PID': ('C2:CP:PID', 'c2.cp.pid'),
              'Reactome': ('C2:CP:REACTOME', 'c2.cp.reactome'),
              'WikiPathways': ('C2:CP:WIKIPATHWAYS', 'c2.cp.wikipathways'),
              'miRDB': ('C3:MIR:MIRDB', 'c3.mir.mirdb'),
              'GTRD': ('C3:TFT:GTRD', 'c3.tft.gtrd'),
              '3CA': ('C4:3CA', 'c4.3ca'),
              'CGN': ('C4:CGN', 'c4.cgn'),
              'CM': ('C4:CM', 'c4.cm'),
              'GO:BP': ('C5:GO:BP', 'c5.go.bp'),
              'GO:CC': ('C5:GO:CC', 'c5.go.cc'),
              'GO:MF': ('C5:GO:MF', 'c5.go.mf'),
              'HPO': ('C5:HPO', 'c5.hpo'),
              'C6:Onco': ('C6', 'c6.all'),
              'C7:ImmuneSigDB': ('C7:IMMUNESIGDB', 'c7.immunesigdb'),
              'C7:VAX': ('C7:VAX', 'c7.vax'),
              'C8:Cell_type': ('C8', 'c8.all'),
              'C9:Comp_Perturb': ('C9', 'c9.all')},
    'Mouse': {'Hallmark': ('MH', 'mh.all'),
              'Positional': ('M1', 'm1.all'),
              'CGP': ('M2:CGP', 'm2.cgp'),
              'BioCarta': ('M2:CP:BIOCARTA', 'm2.cp.biocarta'),
              'Reactome': ('M2:CP:REACTOME', 'm2.cp.reactome'),
              'WikiPathways': ('M2:CP:WIKIPATHWAYS', 'm2.cp.wikipathways'),
              'miRDB': ('M3:MIRDB', 'm3.mirdb'),
              'GTRD': ('M3:GTRD', 'm3.gtrd'),
              'GO:BP': ('M5:GO:BP', 'm5.go.bp'),
              'GO:CC': ('M5:GO:CC', 'm5.go.cc'),
              'GO:MF': ('M5:GO:MF', 'm5.go.mf'),
              'MP:General': ('M5:MPO:MP', 'm5.mpo.mp'),
              'MPT:Tumor': ('M5:MPO:MPT', 'm5.mpo.mpt'),
              'M7:Immuno': ('M7', 'm7.all'),
              'M8:Cell_type': ('M8', 'm8.all')}}

# KEGG pathways from the KEGG REST API (https://www.kegg.jp/kegg/rest/keggapi.html); KEGG gene IDs of
# human (hsa) and mouse (mmu) genes are NCBI gene IDs
KEGG_ORGANISMS = {'Human': 'hsa', 'Mouse': 'mmu'}
KEGG_URL = 'https://rest.kegg.jp'

# All gene set databases of enrich_msigdb()
GENE_SET_DATABASES = {organism: list(databases) + ['KEGG'] for organism, databases in MSIGDB_DATABASES.items()}


@lru_cache(maxsize=None)
def network_degrees(organism):
    '''
    Gene-level degrees and edge counts of the PPI/DDI and structural PPI networks.

    Returns
    -------
    dict: {DEGREE_DDI: (degree dict, number of edges), DEGREE_STRUCTURAL: (degree dict, number of edges)}
    '''
    # PPI/DDI: collapse domain nodes ('<entrez>/<pfam>') of the joint graph to genes
    ddi = nx.Graph()
    ddi.add_edges_from((u.split('/')[0], v.split('/')[0]) for u, v in network[organism].edges())

    structural = ddi.copy()

    # domain-motif interactions (ELM)
    elm = elm_interactions[organism]
    structural.add_edges_from(zip(elm['Interator gene 1'].astype(str), elm['Interator gene 2'].astype(str)))

    # co-resolved PDB interactions: gene with the residue (Ensembl ID or symbol) -> partner (Entrez ID)
    mapping = database_mapping[organism]
    ensembl_to_entrez = dict(zip(mapping['Gene stable ID'], mapping['NCBI gene ID'].astype(str)))
    symbol_to_entrez = dict(zip(mapping['Gene name'], mapping['NCBI gene ID'].astype(str)))
    residues = pdb[organism].dropna(subset=['entrezgene'])
    for gene, symbol, partner in zip(residues['Gene stable ID'], residues['symbol'], residues['entrezgene']):
        gene = ensembl_to_entrez.get(gene, symbol_to_entrez.get(symbol))
        if gene is not None:
            structural.add_edge(gene, str(partner))

    # networkx counts self-interactions twice in the degree, as in the built-in tables
    return {DEGREE_DDI: (dict(ddi.degree()), ddi.number_of_edges()),
            DEGREE_STRUCTURAL: (dict(structural.degree()), structural.number_of_edges())}


def gene_set_table(gene_sets, organism):
    '''
    Build a pathway table (same format as the built-in ones) with degrees computed on the bundled networks.

    Parameters
    ----------
    gene_sets: iterable of (name, ID, source, list of Entrez gene IDs)
    organism: 'Human' or 'Mouse'

    Returns
    -------
    pd.DataFrame, dict {degree column: number of network edges}
    '''
    degrees = network_degrees(organism)
    rows = []
    for name, set_id, source, genes in gene_sets:
        genes = list(dict.fromkeys(str(g) for g in genes))
        row = {'pathway': name, 'external_id': set_id, 'source': source, 'entrez_gene_ids': genes}
        for column, (degree, _) in degrees.items():
            row[column] = sum(degree.get(g, 0) for g in genes)
        rows.append(row)
    table = pd.DataFrame(rows, columns=['pathway', 'external_id', 'source', 'entrez_gene_ids', DEGREE_STRUCTURAL, DEGREE_DDI])
    return table, {column: n for column, (_, n) in degrees.items()}


def read_symbol_table(organism, chip=None):
    '''
    Gene symbol -> list of NCBI gene IDs.

    chip: optional path to an MSigDB "NCBI_Gene_ID" chip file (e.g. Human_NCBI_Gene_ID_MSigDB.v2026.2.Hs.chip)
          matching the release of the gene sets. By default the table bundled with NEASE (MSigDB 2026.2) is used.
    '''
    path = chip if chip is not None else SYMBOL_TABLES[organism]
    table = pd.read_csv(path, sep='\t', dtype=str, comment='#', usecols=[0, 1])
    table.columns = ['NCBI gene ID', 'Gene Symbol']
    table = table.dropna()
    return table.groupby('Gene Symbol')['NCBI gene ID'].apply(list).to_dict()


def cache_dir():
    '''
    Folder of the downloaded gene sets: $NEASE_CACHE, by default ~/.cache/nease ($XDG_CACHE_HOME/nease).
    Delete it (or its kegg folder) to download the gene sets again, e.g. to get the current KEGG pathways.
    '''
    if os.environ.get('NEASE_CACHE'):
        return os.environ['NEASE_CACHE']
    return os.path.join(os.environ.get('XDG_CACHE_HOME') or os.path.join(os.path.expanduser('~'), '.cache'), 'nease')


def download(url, path):
    # download once; a partial download never replaces the file
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        print('Downloading ' + url)
        part = f'{path}.{os.getpid()}.part'
        with urllib.request.urlopen(url, timeout=600) as response, open(part, 'wb') as f:
            shutil.copyfileobj(response, f)
        if os.path.getsize(part) == 0:
            os.remove(part)
            raise ValueError(f'Empty answer from {url}')
        os.replace(part, path)
    return path


def msigdb_file(database, organism):
    # local path of the MSigDB file of a database, downloaded if needed
    _, name = MSIGDB_DATABASES[organism][database]
    species = MSIGDB_SPECIES[organism]
    url = MSIGDB_URL.format(release=MSIGDB_RELEASE, species=species, file=name)
    return download(url, os.path.join(cache_dir(), 'msigdb', f'{MSIGDB_RELEASE}.{species}', os.path.basename(url)))


def read_msigdb_json(gene_sets, organism, databases, symbol_to_entrez):
    '''
    Gene sets of the selected databases from MSigDB JSON bundles.

    Parameters
    ----------
    gene_sets: list of str
        Paths of MSigDB JSON bundles (one collection, e.g. 'c2.cp.reactome.v2026.2.Hs.json', or several, e.g. 'msigdb.v2026.2.Hs.json').
    organism: 'Human' or 'Mouse'
        Human (Hs) bundles can only be used with human data and mouse (Mm) bundles with mouse data.
    databases: list of database names of MSIGDB_DATABASES
    symbol_to_entrez: dict gene symbol -> list of NCBI gene IDs (read_symbol_table)

    Returns
    -------
    list of (name, ID, database, list of NCBI gene IDs), set of gene symbols
    '''
    species = {'Human': '/human/', 'Mouse': '/mouse/'}[organism]
    # MSigDB collection -> database name
    selected = {MSIGDB_DATABASES[organism][database][0]: database for database in databases}

    sets, seen, symbols = [], set(), set()
    for path in gene_sets:
        with open(path) as f:
            bundle = json.load(f)
        for name, info in bundle.items():
            if species not in info.get('msigdbURL', species):
                raise ValueError(f'{os.path.basename(path)} does not contain {organism.lower()} gene sets '
                                 f'(found {info.get("msigdbURL")}). Use the MSigDB {MSIGDB_SPECIES[organism]} bundles for {organism} data.')
            set_id = info.get('systematicName') or name
            if set_id in seen or info['collection'] not in selected:
                continue
            seen.add(set_id)
            genes = []
            for symbol in info['geneSymbols']:
                symbols.add(symbol)
                genes.extend(symbol_to_entrez.get(symbol, []))
            sets.append((name, set_id, selected[info['collection']], genes))
    return sets, symbols


def kegg_gene_sets(organism):
    '''
    KEGG pathways of the organism from the KEGG REST API (downloaded once, then read from the cache).

    Returns
    -------
    list of (name, ID, 'KEGG', list of NCBI gene IDs), date of the download
    '''
    code = KEGG_ORGANISMS[organism]
    folder = os.path.join(cache_dir(), 'kegg')
    names = download(f'{KEGG_URL}/list/pathway/{code}', os.path.join(folder, f'{code}_pathways.tsv'))
    links = download(f'{KEGG_URL}/link/{code}/pathway', os.path.join(folder, f'{code}_pathway_genes.tsv'))

    # 'path:mmu00010<TAB>mmu:103988' -> pathway ID -> NCBI gene IDs
    genes = {}
    with open(links) as f:
        for line in f:
            if line.strip():
                pathway, gene = line.rstrip('\n').split('\t')
                genes.setdefault(pathway.split(':')[-1], []).append(gene.split(':')[-1])

    # 'mmu00010<TAB>Glycolysis / Gluconeogenesis - Mus musculus (house mouse)'
    sets = []
    with open(names) as f:
        for line in f:
            if line.strip():
                pathway, name = line.rstrip('\n').split('\t')
                sets.append((name.rsplit(' - ', 1)[0], pathway, 'KEGG', genes.get(pathway, [])))

    date = datetime.date.fromtimestamp(os.path.getmtime(links)).isoformat()
    return sets, date


def load_gene_sets(databases, organism, gene_sets=None, chip=None):
    '''
    Pathway table (gene_set_table) of the selected gene set databases.

    Parameters
    ----------
    databases: list of names of GENE_SET_DATABASES[organism]
    organism: 'Human' or 'Mouse'
    gene_sets: str or list of str, optional
        MSigDB JSON bundle(s) to read the MSigDB databases from instead of downloading them,
        e.g. the complete 'msigdb.v2026.2.Hs.json'.
    chip: str, optional
        MSigDB NCBI Gene ID chip file used to map gene symbols to NCBI gene IDs (see read_symbol_table).

    Returns
    -------
    pd.DataFrame (pathway table), dict {degree column: number of network edges}, dict (summary)
    '''
    msigdb = [database for database in databases if database in MSIGDB_DATABASES[organism]]
    network_genes = network_degrees(organism)[DEGREE_STRUCTURAL][0]

    sets, summary = [], {}
    if msigdb:
        if gene_sets is None:
            gene_sets = [msigdb_file(database, organism) for database in msigdb]
        elif isinstance(gene_sets, (str, os.PathLike)):
            gene_sets = [gene_sets]
        symbol_to_entrez = read_symbol_table(organism, chip)
        msigdb_sets, symbols = read_msigdb_json(gene_sets, organism, msigdb, symbol_to_entrez)
        sets += msigdb_sets
        summary['gene symbols'] = len(symbols)
        summary['gene symbols in the network'] = sum(any(g in network_genes for g in symbol_to_entrez.get(s, [])) for s in symbols)

    if 'KEGG' in databases:
        kegg_sets, summary['KEGG download date'] = kegg_gene_sets(organism)
        sets += kegg_sets

    table, n = gene_set_table(sets, organism)
    summary['gene sets'] = len(table)
    # databases in the selected order
    summary['databases'] = [database for database in databases if database in set(table['source'])]
    return table, n, summary
