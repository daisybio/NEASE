"""Step 10: human-mouse one-to-one orthologs and residue position maps.

Few mouse protein interactions are measured directly; human interactions, structures and motifs are transferred
to mouse through one-to-one orthologs (Ensembl Compara homologies of the Ensembl release). Residue positions are
transferred with global alignments (pairs with >= orthologs.min_identity % identity) of every human UniProt entry of
the gene (interfaces and motifs are given on the entry a structure or motif refers to, which is often not the one of
the canonical transcript) to one mouse protein per gene (the canonical transcript's entry, otherwise Swiss-Prot,
otherwise the longest). UniProt entries identical to an Ensembl translation or aligned to one (step 2) are used.

Writes to <work_dir>/intermediate/:
    orthologs.parquet       human_gene, mouse_gene, human_entrez, mouse_entrez, perc_id
    position_maps.parquet   human_uniprot, mouse_uniprot, identity, human_pos, mouse_pos (aligned positions)

usage: python 10_orthologs.py --config config.yaml
"""
import argparse
from multiprocessing import Pool

import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

from common import input_path, load_config, log, matrix_letters, paths, placeable_proteins, species_dir


def aligner():
    a = Align.PairwiseAligner(mode='global', substitution_matrix=substitution_matrices.load('BLOSUM62'))
    a.open_gap_score, a.extend_gap_score = -11, -1
    return a


ALIGNER = None


def align(job):
    global ALIGNER
    if ALIGNER is None:
        ALIGNER = aligner()
    human_uniprot, mouse_uniprot, human_seq, mouse_seq = job
    human_seq, mouse_seq = matrix_letters(human_seq), matrix_letters(mouse_seq)  # e.g. U (selenocysteine) -> X
    alignment = ALIGNER.align(human_seq, mouse_seq)[0]
    human_pos, mouse_pos, identical = [], [], 0
    for (h0, h1), (m0, m1) in zip(*alignment.aligned):
        for h, m in zip(range(h0, h1), range(m0, m1)):
            human_pos.append(h + 1)
            mouse_pos.append(m + 1)
            identical += human_seq[h] == mouse_seq[m]
    return human_uniprot, mouse_uniprot, round(100 * identical / max(len(human_seq), len(mouse_seq)), 1), human_pos, mouse_pos


def gene_proteins(cfg, species, one_per_gene):
    '''
    UniProt entries (not isoform entries) of every gene that can be placed on exons (identical to an Ensembl
    translation, or aligned to one). one_per_gene: only the first of identical, canonical transcript's entry,
    Swiss-Prot, longest (genes whose canonical transcript has no identical UniProt entry keep a protein).
    '''
    proteins = placeable_proteins(cfg, species)
    proteins = proteins[~proteins['uniprot'].str.contains('-')]
    proteins = proteins.sort_values(['gene_id', 'identical', 'canonical', 'swissprot', 'length', 'uniprot'],
                                    ascending=[True, False, False, False, False, True], kind='stable')
    proteins = proteins.drop_duplicates(['gene_id', 'uniprot'])
    if one_per_gene:
        proteins = proteins.drop_duplicates('gene_id')
    return proteins[['gene_id', 'uniprot', 'sequence']]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    args = parser.parse_args()
    cfg = load_config(args.config)
    out = paths(cfg)['intermediate']
    if 'Mouse' not in cfg['species']:
        log('only needed for the mouse network: nothing to do')
        raise SystemExit(0)

    parts = []
    # homologies of the mouse genes (the human-mouse pairs are only listed in the mouse file)
    columns = ['gene_stable_id', 'homology_type', 'homology_gene_stable_id', 'homology_species', 'homology_identity']
    for chunk in pd.read_csv(input_path(cfg, 'compara_homologies'), sep='\t', usecols=columns, dtype=str, chunksize=1_000_000):
        chunk = chunk[(chunk['homology_species'] == 'homo_sapiens') & (chunk['homology_type'] == 'ortholog_one2one')]
        parts.append(chunk.rename(columns={'gene_stable_id': 'mouse_gene', 'homology_gene_stable_id': 'human_gene', 'homology_identity': 'perc_id'}))
    homologs = pd.concat(parts, ignore_index=True)[['human_gene', 'mouse_gene', 'perc_id']].drop_duplicates()

    for species, column in [('Human', 'human'), ('Mouse', 'mouse')]:
        genes = pd.read_parquet(species_dir(cfg, species) / 'genes.parquet', columns=['gene_id', 'entrez']).dropna().drop_duplicates()
        homologs = homologs.merge(genes.rename(columns={'gene_id': f'{column}_gene', 'entrez': f'{column}_entrez'}), on=f'{column}_gene', how='left')
    homologs.reset_index(drop=True).to_parquet(out / 'orthologs.parquet')
    log(f'{homologs["human_gene"].nunique()} one-to-one orthologs '
        f'({homologs.dropna()[["human_entrez", "mouse_entrez"]].drop_duplicates().shape[0]} NCBI gene ID pairs)')

    human = gene_proteins(cfg, 'Human', one_per_gene=False).rename(columns={'gene_id': 'human_gene', 'uniprot': 'human_uniprot', 'sequence': 'human_seq'})
    mouse = gene_proteins(cfg, 'Mouse', one_per_gene=True).rename(columns={'gene_id': 'mouse_gene', 'uniprot': 'mouse_uniprot', 'sequence': 'mouse_seq'})
    pairs = homologs[['human_gene', 'mouse_gene']].drop_duplicates().merge(human, on='human_gene').merge(mouse, on='mouse_gene')
    jobs = list(pairs[['human_uniprot', 'mouse_uniprot', 'human_seq', 'mouse_seq']].itertuples(index=False, name=None))
    log(f'aligning {len(jobs)} protein pairs')
    with Pool(cfg['orthologs']['threads']) as pool:
        maps = pd.DataFrame(pool.map(align, jobs, chunksize=50), columns=['human_uniprot', 'mouse_uniprot', 'identity', 'human_pos', 'mouse_pos'])
    maps = maps[maps['identity'] >= cfg['orthologs']['min_identity']]
    maps.reset_index(drop=True).to_parquet(out / 'position_maps.parquet')
    log(f'{len(maps)} position maps for {maps["human_uniprot"].nunique()} human UniProt entries of '
        f'{pairs.loc[pairs["human_uniprot"].isin(maps["human_uniprot"]), "human_gene"].nunique()} genes '
        f'(median identity {maps["identity"].median()}%)')
