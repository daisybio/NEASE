"""Step 14: compare a built network with the one shipped with NEASE (and optionally test it against measured data).

Reports per species: size of every layer, overlap with the shipped network (gene and domain level), which
evidence recovers the shipped edges, the network sizes and the correlation of pathway degrees.

--benchmark: a TSV of measured effects of exon skipping on interactions (e.g. from isoform-resolved interaction
screens such as Yang et al., Cell 2016) with the columns
    species, gene (Ensembl gene ID), start, end (genomic coordinates of the exon), partner (NCBI gene ID),
    lost (1: the isoform without the exon lost the interaction, 0: it kept it)
For the shipped and the new network, an interaction counts as predicted to be affected when the exon overlaps a
domain linked to the partner in the joint graph, a motif or a residue interface with the partner. A useful
network predicts lost interactions more often than kept ones.

Writes <work_dir>/output/<name>/validation.md

usage: python 14_validate.py --config config.yaml [--name nease_data] [--benchmark benchmark.tsv]
"""
import argparse
import json
import pickle

import pandas as pd

from common import SPECIES, load_config, log, nease_file, overlapping, paths

SHIPPED_SIZES = {'Human': {'Degree in the PPI/DDI': 52467, 'Degree in the structural PPI': 60235},
                 'Mouse': {'Degree in the PPI/DDI': 13348, 'Degree in the structural PPI': 31144}}


def load(folder, species):
    files = SPECIES[species]['nease']
    with open(folder / files['graph'], 'rb') as f:
        joint = pickle.load(f)
    with open(folder / files['ppi'], 'rb') as f:
        ppi = pickle.load(f)
    return dict(joint=joint, ppi=ppi, elm=pd.read_pickle(folder / files['elm']),
                elm_interactions=pd.read_pickle(folder / files['elm_interactions']),
                pdb=pd.read_pickle(folder / files['pdb']), pathways=pd.read_pickle(folder / files['pathways']))


def edge_set(graph, level):
    if level == 'gene':
        return {frozenset((u.split('/')[0], v.split('/')[0])) for u, v in graph.edges()}
    return {frozenset((u, v)) for u, v in graph.edges()}


def overlap(old, new):
    both = len(old & new)
    return f'{len(old)} shipped, {len(new)} new, {both} shared ({100 * both / max(len(old), 1):.1f}% of shipped kept, Jaccard {both / max(len(old | new), 1):.2f})'


def predicted(network, exon_domains, rows):
    '''For every benchmark row: is the interaction with the partner predicted to be affected by the exon?'''
    joint = network['joint']
    elm, elm_interactions, pdb = network['elm'], network['elm_interactions'], network['pdb']
    motif_partners = elm_interactions.groupby(['Interator gene 1', 'Elm id of gene 1'])['Interator gene 2'].agg(lambda s: set(map(str, s)))
    out = []
    for row in rows.itertuples(index=False):
        hit = False
        domains = exon_domains[(exon_domains['Gene stable ID'] == row.gene) & (exon_domains['Genomic coding start'] <= row.end)
                               & (exon_domains['Genomic coding end'] >= row.start)]
        for node in set(domains['NCBI gene ID'].astype(str) + '/' + domains['Pfam ID']):
            if node in joint and any(n.split('/')[0] == row.partner for n in joint.neighbors(node)):
                hit = True
        motifs = elm[(elm['Gene stable ID'] == row.gene) & (elm['Genomic coding start'] <= row.end) & (elm['Genomic coding end'] >= row.start)]
        for gene, motif in zip(motifs['entrezgene'], motifs['ELMIdentifier']):
            hit |= row.partner in motif_partners.get((int(gene), motif), set())
        hit |= ((pdb['Gene stable ID'] == row.gene) & (pdb['Genomic coding start'] <= row.end) & (pdb['Genomic coding end'] >= row.start)
                & (pdb['entrezgene'].astype(str) == row.partner)).any()
        out.append(hit)
    return pd.Series(out, index=rows.index)


def benchmark_report(cfg, species, networks, path):
    table = pd.read_csv(path, sep='\t', dtype={'partner': str})
    table = table[table['species'] == species]
    if table.empty:
        return []
    exon_domains = pd.read_pickle(nease_file(cfg, species, 'exon_domains'), compression='gzip')
    lines = [f'\n### Benchmark ({len(table)} exon-partner pairs, {int(table["lost"].sum())} lost)\n',
             '| network | lost interactions predicted | kept interactions predicted |', '|---|---|---|']
    for name, network in networks.items():
        hit = predicted(network, exon_domains, table)
        lost, kept = hit[table['lost'] == 1], hit[table['lost'] == 0]
        lines.append(f'| {name} | {lost.sum()}/{len(lost)} ({100 * lost.mean():.1f}%) | {kept.sum()}/{len(kept)} ({100 * kept.mean():.1f}%) |')
    return lines


def report(cfg, species, folder, benchmark):
    shipped = load(cfg['nease_repo'] / 'nease' / 'data', species)
    new = load(folder / 'nease_data', species)
    info = json.loads((folder / f'network_info_{species}.json').read_text())
    lines = [f'## {species}\n']
    lines.append(f'- Joint graph, domain level: {overlap(edge_set(shipped["joint"], "domain"), edge_set(new["joint"], "domain"))}')
    lines.append(f'- Joint graph, gene level: {overlap(edge_set(shipped["joint"], "gene"), edge_set(new["joint"], "gene"))}')
    lines.append(f'- Protein interaction network: {overlap(edge_set(shipped["ppi"], "gene"), edge_set(new["ppi"], "gene"))}')
    old_edges = edge_set(shipped['joint'], 'domain')
    recovered = pd.Series([e for u, v, d in new['joint'].edges(data=True) if frozenset((u, v)) in old_edges for e in d.get('evidence', [])])
    lines.append(f'- Evidence of the shipped domain edges that are recovered: {recovered.value_counts().to_dict()}')
    lines.append(f'- Evidence of all new domain edges: {info["joint_graph"]["edges_by_evidence"]}')
    old_motifs = set(zip(shipped['elm']['Exon stable ID'], shipped['elm']['ELMIdentifier']))
    new_motifs = set(zip(new['elm']['Exon stable ID'], new['elm']['ELMIdentifier']))
    lines.append(f'- Motif exons (exon, ELM class): {overlap(old_motifs, new_motifs)}')
    key = ['Interator gene 1', 'Elm id of gene 1', 'Interator gene 2', 'Domain of gene 2']
    lines.append(f'- Motif-domain interactions: {overlap(set(map(tuple, shipped["elm_interactions"][key].astype(str).values)), set(map(tuple, new["elm_interactions"][key].astype(str).values)))}')
    old_pdb, new_pdb = (t.assign(partner=t['entrezgene'].astype(str).str.replace(r'\.0$', '', regex=True)) for t in (shipped['pdb'], new['pdb']))
    kept = overlapping(old_pdb, new_pdb, ['Gene stable ID', 'partner'])
    lines.append(f'- Residue-level interfaces (exon, partner): {len(old_pdb)} shipped rows, {len(new_pdb)} new rows; {int(kept.sum())} '
                 f'shipped rows ({100 * kept.mean():.1f}%) have a new row with the same gene and partner and overlapping '
                 'coordinates (how NEASE matches them; exon IDs differ between Ensembl releases)')
    lines.append(f'- Network sizes (edges) shipped: {SHIPPED_SIZES[species]}, new: {info["network_edges"]}')
    size_file = folder / 'nease_data' / 'network' / f'network_size_{species.lower()}.json'
    if size_file.exists():
        sizes = json.loads(size_file.read_text())
        ok = sizes['joint_graph_edges'] == new['joint'].number_of_edges() and sizes['network_edges'] == info['network_edges']
        lines.append(f'- Network sizes NEASE uses for this network ({size_file.name}): {sizes["network_edges"]}'
                     + ('' if ok else ' -- WARNING: does not match the joint graph / network_info; rerun step 13'))
    else:
        lines.append(f'- WARNING: {size_file.name} missing (network built before the size file existed: rerun step 13). '
                     'NEASE would use the sizes of its built-in network for the enrichment test.')
    for column in SHIPPED_SIZES[species]:
        # Spearman correlation as the Pearson correlation of the ranks (pandas' method='spearman' needs scipy)
        degrees = pd.concat([shipped['pathways'][column], new['pathways'][column]], axis=1).dropna()
        corr = degrees.rank().corr().iloc[0, 1]
        lines.append(f'- Pathway degrees "{column}": Spearman correlation shipped vs new {corr:.3f}')
    if benchmark:
        lines += benchmark_report(cfg, species, {'shipped': shipped, 'new': new}, benchmark)
    return lines


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    parser.add_argument('--name', default='nease_data')
    parser.add_argument('--benchmark')
    args = parser.parse_args()
    cfg = load_config(args.config)
    folder = paths(cfg)['output'] / args.name
    lines = [f'# Validation of {folder}\n']
    for species in cfg['species']:
        lines += report(cfg, species, folder, args.benchmark) + ['']
    (folder / 'validation.md').write_text('\n'.join(lines))
    print('\n'.join(lines))
    log(f'written: {folder / "validation.md"}')
