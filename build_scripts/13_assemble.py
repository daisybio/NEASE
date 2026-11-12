"""Step 13: assemble the NEASE network files (per species) from the evidence selected in config network.

Domain-level network (joint graph): nodes '<NCBI gene ID>/<Pfam ID>' for the domains of NEASE's exon -> domain
table; an edge links two domains of interacting proteins when
    pdb_pair / af_pair (+ _ortholog)   the domains are in contact in a structure / confident model of the pair
    i3d_structure_pair / i3d_model_pair (+ _ortholog)   the same in an Interactome3D structure / homology model
    pdb_homolog_pair (+ _ortholog)     the same in a PDB structure of closely related proteins (other species or
                                       paralogs, >= homologs.min_identity % identical; step 17)
    humanppi_pair (+ _ortholog)        the same in a confident model of the predicted human interactions of
                                       Zhang et al., Science 2025 (step 18)
    boltz_pair (+ _ortholog)           the same in a Boltz-2 model of a measured direct interaction (Ille et al.
                                       2025; step 19)
    ddi_*                              the proteins interact physically and the Pfam pair is in a domain-domain
                                       catalogue (ddi_pdb: computed from the PDB; ddi_3did2020; ddi_ppidm_*)
Every edge carries the attribute 'evidence' (list). Residue-level interfaces (network.interface_evidence) give
NEASE's "pdb" table, the motif tables come from 12_elm.py. Pathway degrees and network sizes are recomputed.

Writes <work_dir>/output/<name>/ in the layout of nease/data (network/, database/, pathways/), a
network_info_<species>.json with the network sizes and the evidence used, and details/ with all tables
including their evidence. nease_data/network/network_size_<species>.json holds the network sizes for NEASE's
enrichment test: NEASE uses them instead of the sizes of its built-in network when this folder is loaded
(NEASE_DATA_DIR=<work_dir>/output/<name>/nease_data, or the files copied into nease/data together).

usage: python 13_assemble.py --config config.yaml [--name nease_data] [--evidence pdb_pair ddi_pdb ...]
"""
import argparse
import json
import zipfile
from datetime import datetime

import networkx as nx
import pandas as pd

from common import (forward_pfam, gene_domains, input_path, load_config, log, nease_file, paths, pfam_forward,
                    read_uniprot_entrez, save_pickle, SPECIES, species_dir)
from exon_mapping import ResidueExonMapper

DEGREE_DDI, DEGREE_STRUCTURAL = 'Degree in the PPI/DDI', 'Degree in the structural PPI'
INTERFACE_KINDS = ('pdb', 'af', 'i3d_structure', 'i3d_model', 'pdb_homolog', 'humanppi', 'boltz')
PAIR_TIER = {kind + suffix: f'{kind}_pair{suffix}' for kind in INTERFACE_KINDS for suffix in ('', '_ortholog')}
PPIDM = {'ddi_ppidm_gold': 'Gold', 'ddi_ppidm_silver': 'Silver', 'ddi_ppidm_bronze': 'Bronze'}


def ddi_catalogues(cfg, evidence, forward):
    '''Pfam pairs of the selected domain-domain catalogues, both orientations.'''
    parts = []
    if 'ddi_pdb' in evidence:
        ddi = pd.read_parquet(paths(cfg)['intermediate'] / 'pdb' / 'ddi_pdb.parquet')
        ddi = ddi[ddi['n_entries'] >= cfg['network']['min_ddi_entries']]
        parts.append(ddi[['pfam_1', 'pfam_2']].assign(evidence='ddi_pdb'))
    with zipfile.ZipFile(input_path(cfg, 'ppidm')) as archive:
        if 'ddi_3did2020' in evidence:
            with archive.open('PPIDMdata_2021/EvaluationData/3did2020.csv') as f:
                parts.append(pd.read_csv(f, usecols=['D1', 'D2']).set_axis(['pfam_1', 'pfam_2'], axis=1).assign(evidence='ddi_3did2020'))
        for tier, name in PPIDM.items():
            if tier in evidence:
                with archive.open(f'PPIDMdata_2021/OutputData/PPIDM_{name}DDIs.csv') as f:
                    parts.append(pd.read_csv(f, sep=';', usecols=['D1', 'D2']).set_axis(['pfam_1', 'pfam_2'], axis=1).assign(evidence=tier))
    if not parts:
        return pd.DataFrame(columns=['pfam_1', 'pfam_2', 'evidence'])
    ddi = forward_pfam(pd.concat(parts, ignore_index=True), ['pfam_1', 'pfam_2'], forward)
    both = pd.concat([ddi, ddi.rename(columns={'pfam_1': 'pfam_2', 'pfam_2': 'pfam_1'})], ignore_index=True)
    both = both.drop_duplicates().reset_index(drop=True)
    # how much every catalogue adds to the others (Pfam pairs, one orientation)
    pairs = ddi.assign(pair=[frozenset(p) for p in zip(ddi['pfam_1'], ddi['pfam_2'])]).groupby('evidence')['pair'].agg(set)
    for tier, tier_pairs in pairs.items():
        others = set().union(*(p for t, p in pairs.items() if t != tier))
        log(f'  {tier}: {len(tier_pairs)} Pfam pairs, {len(tier_pairs - others)} not in the other selected catalogues')
    return both


def interfaces(cfg, species, evidence):
    '''Residue-pair interfaces between proteins of the species: uniprot_1, uniprot_2, pos_1, pos_2, evidence.'''
    parts = []
    human_parts = []
    folder = paths(cfg)['intermediate']
    for kind in INTERFACE_KINDS:
        wanted = kind in evidence or f'{kind}_ortholog' in evidence
        if not wanted:
            continue
        if kind == 'pdb':
            table = pd.read_parquet(folder / 'pdb' / 'interfaces_pdb.parquet')
        elif kind in ('humanppi', 'boltz'):
            path = folder / kind / 'interfaces.parquet'
            if not path.exists():
                log(f'  no {kind} interfaces found (step {18 if kind == "humanppi" else 19} not run?)')
                continue
            table = pd.read_parquet(path)
        elif kind == 'pdb_homolog':
            path = folder / 'pdb' / 'interfaces_homolog.parquet'
            if not path.exists():
                log('  no interfaces of related proteins found (step 17 not run?)')
                continue
            table = pd.read_parquet(path)
        elif kind == 'af':
            files = sorted((folder / 'alphafold' / 'contacts').glob('task_*.parquet'))
            if not files:
                log('  no AlphaFold interfaces found (steps 7-9 not run?)')
                continue
            table = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
        else:
            path = folder / 'interactome3d' / 'interfaces.parquet'
            if not path.exists():
                log(f'  no Interactome3D interfaces found for {kind} (step 15 not run?)')
                continue
            table = pd.read_parquet(path)
            table = table[table['kind'] == kind]
        table = table[['species', 'uniprot_1', 'uniprot_2', 'pos_1', 'pos_2']].assign(evidence=kind)
        if kind in evidence:
            parts.append(table[table['species'] == species])
        if species == 'Mouse' and f'{kind}_ortholog' in evidence:
            human_parts.append(table[table['species'] == 'Human'].assign(evidence=f'{kind}_ortholog'))
    if human_parts:
        parts.append(transfer_interfaces(cfg, pd.concat(human_parts, ignore_index=True)))
    if not parts:
        return pd.DataFrame(columns=['uniprot_1', 'uniprot_2', 'pos_1', 'pos_2', 'evidence'])
    return pd.concat(parts, ignore_index=True).drop(columns='species', errors='ignore').reset_index(drop=True)


def transfer_interfaces(cfg, human):
    '''Human interfaces on the mouse orthologs (residue pairs with both positions aligned).'''
    maps = pd.read_parquet(paths(cfg)['intermediate'] / 'position_maps.parquet')
    maps = {h: (m, dict(zip(hp, mp))) for h, m, hp, mp in maps[['human_uniprot', 'mouse_uniprot', 'human_pos', 'mouse_pos']].itertuples(index=False)}
    rows = []
    for u1, u2, pos_1, pos_2, evidence in human[['uniprot_1', 'uniprot_2', 'pos_1', 'pos_2', 'evidence']].itertuples(index=False):
        if u1 in maps and u2 in maps:
            (m1, a1), (m2, a2) = maps[u1], maps[u2]
            pairs = [(a1[p], a2[q]) for p, q in zip(pos_1, pos_2) if p in a1 and q in a2]
            if pairs:
                rows.append((m1, m2, [p for p, _ in pairs], [q for _, q in pairs], evidence))
    log(f'  {len(rows)} interfaces transferred from human')
    return pd.DataFrame(rows, columns=['uniprot_1', 'uniprot_2', 'pos_1', 'pos_2', 'evidence'])


def protein_domains(cfg, species, forward):
    '''
    UniProt accession -> NCBI gene ID, Pfam domain and its range in the UniProt sequence (NEASE's table, via the
    identical transcripts, or via the aligned transcripts of entries without an identical one).
    '''
    table = pd.read_pickle(nease_file(cfg, species, 'exon_domains'), compression='gzip')
    table = table[['Transcript stable ID', 'NCBI gene ID', 'Pfam ID', 'Pfam start', 'Pfam end']].drop_duplicates()
    table.columns = ['transcript_id', 'entrez', 'pfam', 'start', 'end']
    table['entrez'] = table['entrez'].astype(str)
    proteins = pd.read_parquet(species_dir(cfg, species) / 'proteins.parquet', columns=['uniprot', 'transcript_id'])
    parts = [proteins.merge(table, on='transcript_id').drop(columns='transcript_id')]
    path = species_dir(cfg, species) / 'uniprot_maps.parquet'
    if path.exists():
        maps = pd.read_parquet(path, columns=['uniprot', 'transcript_id', 'uni_pos', 'ens_pos'])
        rows = []
        for uniprot, uni_pos, ens_pos, entrez, pfam, start, end in maps.merge(table, on='transcript_id')[
                ['uniprot', 'uni_pos', 'ens_pos', 'entrez', 'pfam', 'start', 'end']].itertuples(index=False):
            inside = uni_pos[(ens_pos >= start) & (ens_pos <= end)]  # domain range on the transcript -> UniProt positions
            if len(inside):
                rows.append((uniprot, entrez, pfam, int(inside.min()), int(inside.max())))
        parts.append(pd.DataFrame(rows, columns=['uniprot', 'entrez', 'pfam', 'start', 'end']))
    table = pd.concat(parts, ignore_index=True).drop_duplicates()
    return forward_pfam(table, ['pfam'], forward).drop_duplicates()


def pair_domain_edges(cfg, contacts, domains):
    '''Domain pairs in contact (>= pdb.min_domain_contacts residue pairs) in the interfaces.'''
    rows = []
    domains_1 = domains.rename(columns={'uniprot': 'uniprot_1', 'entrez': 'entrez_1', 'pfam': 'pfam_1', 'start': 'start_1', 'end': 'end_1'})
    domains_2 = domains.rename(columns={'uniprot': 'uniprot_2', 'entrez': 'entrez_2', 'pfam': 'pfam_2', 'start': 'start_2', 'end': 'end_2'})
    contacts = contacts.reset_index(drop=True)
    for start in range(0, len(contacts), 5000):
        chunk = contacts.iloc[start:start + 5000].reset_index().rename(columns={'index': 'interface'})
        pairs = chunk[['interface', 'uniprot_1', 'uniprot_2', 'pos_1', 'pos_2', 'evidence']].explode(['pos_1', 'pos_2'])
        pairs = pairs.merge(domains_1, on='uniprot_1')
        pairs = pairs[(pairs['pos_1'] >= pairs['start_1']) & (pairs['pos_1'] <= pairs['end_1'])]
        pairs = pairs.merge(domains_2, on='uniprot_2')
        pairs = pairs[(pairs['pos_2'] >= pairs['start_2']) & (pairs['pos_2'] <= pairs['end_2'])]
        key = ['interface', 'entrez_1', 'pfam_1', 'entrez_2', 'pfam_2', 'evidence']
        counts = pairs.drop_duplicates(key + ['pos_1', 'pos_2']).groupby(key).size().reset_index(name='n')
        rows.append(counts[counts['n'] >= cfg['pdb']['min_domain_contacts']])
    edges = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=['entrez_1', 'pfam_1', 'entrez_2', 'pfam_2', 'evidence'])
    edges['evidence'] = edges['evidence'].map(PAIR_TIER)
    return edges[['entrez_1', 'pfam_1', 'entrez_2', 'pfam_2', 'evidence']].drop_duplicates()


def catalogue_edges(ppi, domains, ddi):
    '''Physical interactions x domains of both proteins x domain-domain catalogues.'''
    rows = []
    d1 = domains.rename(columns={'entrez': 'entrez_1', 'pfam': 'pfam_1'})
    d2 = domains.rename(columns={'entrez': 'entrez_2', 'pfam': 'pfam_2'})
    for start in range(0, len(ppi), 200000):
        chunk = ppi.iloc[start:start + 200000][['entrez_1', 'entrez_2']]
        chunk = chunk.merge(d1, on='entrez_1').merge(d2, on='entrez_2').merge(ddi, on=['pfam_1', 'pfam_2'])
        rows.append(chunk)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=['entrez_1', 'pfam_1', 'entrez_2', 'pfam_2', 'evidence'])


def residue_table(cfg, species, contacts, uniprot_entrez):
    '''NEASE "pdb" table: exons with interface residues and the NCBI gene ID of the partner.'''
    sides = []
    for a, b in (('1', '2'), ('2', '1')):
        side = contacts[[f'uniprot_{a}', f'uniprot_{b}', f'pos_{a}', 'evidence']].copy()
        side.columns = ['uniprot', 'partner', 'position', 'evidence']
        side['position'] = side['position'].map(lambda p: sorted(set(p)))
        sides.append(side.explode('position'))
    residues = pd.concat(sides, ignore_index=True).dropna()
    residues['position'] = residues['position'].astype(int)
    residues = residues.merge(uniprot_entrez.rename(columns={'uniprot': 'partner', 'entrez': 'entrezgene'}), on='partner')
    exons = ResidueExonMapper(cfg, species).map(residues[['uniprot', 'position', 'entrezgene', 'evidence']].drop_duplicates())
    table = exons.groupby(['symbol', 'gene_id', 'exon_id', 'coding_start', 'coding_end', 'entrezgene'], dropna=False, as_index=False)['evidence'].agg(
        lambda s: ','.join(sorted(set(s))))
    return table.rename(columns={'gene_id': 'Gene stable ID', 'exon_id': 'Exon with residue ID', 'coding_start': 'Genomic coding start',
                                 'coding_end': 'Genomic coding end'})


def degrees(joint, elm_interactions, pdb_table, mapping):
    '''
    Gene-level networks exactly as NEASE builds them (nease.gene_sets.network_degrees), so that the network sizes and
    pathway degrees written here are the ones NEASE computes for this network: PPI/DDI (joint graph) and structural
    (+ ELM, residue interfaces). mapping: NEASE's exon -> domain table (gene of the residue: Ensembl ID or symbol).
    '''
    ddi = nx.Graph()
    ddi.add_edges_from((u.split('/')[0], v.split('/')[0]) for u, v in joint.edges())
    structural = ddi.copy()
    structural.add_edges_from(zip(elm_interactions['Interator gene 1'].astype(str), elm_interactions['Interator gene 2'].astype(str)))
    ensembl_to_entrez = dict(zip(mapping['Gene stable ID'], mapping['NCBI gene ID'].astype(str)))
    symbol_to_entrez = dict(zip(mapping['Gene name'], mapping['NCBI gene ID'].astype(str)))
    residues = pdb_table.dropna(subset=['entrezgene'])
    for gene, symbol, partner in zip(residues['Gene stable ID'], residues['symbol'], residues['entrezgene']):
        gene = ensembl_to_entrez.get(gene, symbol_to_entrez.get(symbol))
        if gene is not None:
            structural.add_edge(gene, str(partner))
    return {DEGREE_DDI: (dict(ddi.degree()), ddi.number_of_edges()), DEGREE_STRUCTURAL: (dict(structural.degree()), structural.number_of_edges())}


def assemble(cfg, species, evidence, interface_evidence, out, forward):
    log(f'== {species}')
    folder = species_dir(cfg, species)
    details = out / 'details' / species
    details.mkdir(parents=True, exist_ok=True)
    uniprot_entrez = read_uniprot_entrez(cfg, species)
    domains = gene_domains(cfg, species, forward)

    # structures and models of protein pairs
    contacts = interfaces(cfg, species, set(PAIR_TIER) & (set(interface_evidence) | {k for k, v in PAIR_TIER.items() if v in evidence}))
    log(f'  {len(contacts)} protein pair interfaces', contacts['evidence'].value_counts().to_dict())
    pair_edges = pair_domain_edges(cfg, contacts[contacts['evidence'].map(PAIR_TIER).isin(evidence)],
                                   protein_domains(cfg, species, forward))

    # physical interactions (and pairs seen in structures)
    ppi = pd.read_parquet(folder / 'ppi.parquet')
    implies = cfg['network']['structure_implies_ppi']
    # true: every structure/model tier; a list: these tiers (and their _ortholog transfers to mouse)
    implying = contacts if implies is True else contacts[contacts['evidence'].str.replace('_ortholog', '').isin(implies or [])]
    if len(implying):
        pairs = implying[['uniprot_1', 'uniprot_2', 'evidence']]
        for side in ('1', '2'):
            pairs = pairs.merge(uniprot_entrez.rename(columns={'uniprot': f'uniprot_{side}', 'entrez': f'entrez_{side}'}), on=f'uniprot_{side}')
        pairs = pd.DataFrame({'entrez_1': pairs[['entrez_1', 'entrez_2']].min(axis=1), 'entrez_2': pairs[['entrez_1', 'entrez_2']].max(axis=1),
                              'structure': pairs['evidence']})
        structure = pairs.groupby(['entrez_1', 'entrez_2'])['structure'].agg(lambda s: ','.join(sorted(set(s)))).reset_index()
        ppi = ppi.merge(structure, on=['entrez_1', 'entrez_2'], how='outer')
        ppi['structure'] = ppi['structure'].fillna('')
    ppi.to_parquet(details / 'ppi.parquet')

    ddi = ddi_catalogues(cfg, evidence, forward)
    edges = pd.concat([pair_edges, catalogue_edges(ppi, domains, ddi)], ignore_index=True)
    # node names and one row per edge
    edges['node_1'] = edges['entrez_1'] + '/' + edges['pfam_1']
    edges['node_2'] = edges['entrez_2'] + '/' + edges['pfam_2']
    swap = edges['node_1'] > edges['node_2']
    edges.loc[swap, ['node_1', 'node_2']] = edges.loc[swap, ['node_2', 'node_1']].values
    edges = edges.groupby(['node_1', 'node_2'])['evidence'].agg(lambda s: sorted(set(s))).reset_index()
    edges.to_parquet(details / 'joint_edges.parquet')
    joint = nx.Graph()
    joint.add_edges_from((u, v, {'evidence': e}) for u, v, e in edges.itertuples(index=False))
    log(f'  joint graph: {joint.number_of_nodes()} domains, {joint.number_of_edges()} edges',
        edges['evidence'].explode().value_counts().to_dict())

    ppi_graph = nx.Graph()
    ppi_graph.add_edges_from(zip(ppi['entrez_1'], ppi['entrez_2']))

    # residue level (interfaces) and motifs
    selected = contacts[contacts['evidence'].isin(interface_evidence)]
    pdb_table = residue_table(cfg, species, selected, uniprot_entrez) if len(selected) else pd.DataFrame(
        columns=['symbol', 'Gene stable ID', 'Exon with residue ID', 'Genomic coding start', 'Genomic coding end', 'entrezgene', 'evidence'])
    pdb_table.to_parquet(details / 'pdb.parquet')
    elm = pd.read_parquet(folder / 'elm.parquet')
    elm.to_parquet(details / 'elm.parquet')
    elm_interactions = pd.read_parquet(folder / 'elm_interactions.parquet')
    elm_interactions.to_parquet(details / 'elm_interactions.parquet')
    elm_interactions = pd.DataFrame({'Interator gene 1': elm_interactions['gene_1'].astype('int64'), 'Elm id of gene 1': elm_interactions['elm'],
                                     'Interator gene 2': elm_interactions['gene_2'].astype('int64'), 'Domain of gene 2': elm_interactions['pfam']})
    log(f'  residue level: {len(pdb_table)} exon-partner pairs ({pdb_table["Exon with residue ID"].nunique()} exons); '
        f'motifs: {elm["Exon stable ID"].nunique()} exons, {len(elm_interactions)} motif-domain interactions')

    # pathway degrees and network sizes
    network_degrees = degrees(joint, elm_interactions, pdb_table, pd.read_pickle(nease_file(cfg, species, 'exon_domains'), compression='gzip'))
    pathways = pd.read_pickle(nease_file(cfg, species, 'pathways'))
    for column, (degree, _) in network_degrees.items():
        pathways[column] = pathways['entrez_gene_ids'].map(lambda g: sum(degree.get(str(x), 0) for x in g))
    sizes = {column: n for column, (_, n) in network_degrees.items()}
    log(f'  network sizes: {sizes}')

    files = SPECIES[species]['nease']
    data = out / 'nease_data'
    save_pickle(joint, data / files['graph'])
    save_pickle(ppi_graph, data / files['ppi'])
    save_pickle(elm[['ELMIdentifier', 'ensembl.transcript', 'Gene stable ID', 'Exon stable ID', 'Genomic coding start',
                     'Genomic coding end', 'entrezgene']].drop_duplicates().reset_index(drop=True), data / files['elm'])
    save_pickle(elm_interactions.drop_duplicates().reset_index(drop=True), data / files['elm_interactions'])
    save_pickle(pdb_table[['symbol', 'Gene stable ID', 'Exon with residue ID', 'Genomic coding start', 'Genomic coding end',
                           'entrezgene']].reset_index(drop=True), data / files['pdb'])
    save_pickle(pathways, data / files['pathways'])
    # network sizes for NEASE's enrichment test (nease/load.py: used when they belong to the loaded joint graph)
    (data / 'network' / f'network_size_{species.lower()}.json').write_text(json.dumps(
        {'joint_graph_edges': joint.number_of_edges(), 'network_edges': sizes}, indent=2))
    info = {'species': species, 'created': datetime.now().isoformat(timespec='seconds'), 'network_edges': sizes,
            'evidence': evidence, 'interface_evidence': interface_evidence,
            'joint_graph': {'nodes': joint.number_of_nodes(), 'edges': joint.number_of_edges(),
                            'edges_by_evidence': edges['evidence'].explode().value_counts().to_dict()},
            'ppi_edges': ppi_graph.number_of_edges()}
    (out / f'network_info_{species}.json').write_text(json.dumps(info, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    parser.add_argument('--name', default='nease_data', help='output folder name (e.g. one per evidence selection)')
    parser.add_argument('--evidence', nargs='+', help='overrides network.evidence')
    parser.add_argument('--interface-evidence', nargs='+', help='overrides network.interface_evidence')
    args = parser.parse_args()
    cfg = load_config(args.config)
    out = paths(cfg)['output'] / args.name
    forward = pfam_forward(cfg)
    manifest = paths(cfg)['downloads'] / 'manifest.tsv'
    for species in cfg['species']:
        evidence = args.evidence or cfg['network']['evidence']
        interface_evidence = args.interface_evidence or cfg['network']['interface_evidence']
        if species == 'Human':
            evidence = [e for e in evidence if not e.endswith('_ortholog')]
            interface_evidence = [e for e in interface_evidence if not e.endswith('_ortholog')]
        assemble(cfg, species, evidence, interface_evidence, out, forward)
    if manifest.exists():
        (out / 'downloads_manifest.tsv').write_text(manifest.read_text())
    log(f'done: {out}')
