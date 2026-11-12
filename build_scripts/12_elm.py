"""Step 12: linear motif (ELM) tables in NEASE's format.

    elm               exons that contain a motif instance (true-positive ELM instances, mapped to exons
                      through the transcripts identical to the UniProt sequence)
    ELM_interactions  motif (gene 1, ELM class) - domain (gene 2, Pfam) interactions:
                        curated ELM interactions, and (network.elm_infer_from_ppi) motif instances x physical
                        interaction partners that carry a Pfam domain binding the motif class
Mouse: mouse instances and interactions plus the human ones transferred through orthologs.

Writes <work_dir>/intermediate/<species>/elm.parquet and elm_interactions.parquet (with an evidence column)

usage: python 12_elm.py --config config.yaml
"""
import argparse
import io

import pandas as pd

from common import (SPECIES, gene_domains, input_path, load_config, log, paths, pfam_forward, read_uniprot_entrez,
                    species_dir)
from exon_mapping import ResidueExonMapper


def read_elm_tsv(path):
    with open(path) as f:
        text = ''.join(line for line in f if not line.startswith('#'))
    return pd.read_csv(io.StringIO(text), sep='\t', dtype=str)


def instances(cfg, species):
    table = read_elm_tsv(input_path(cfg, 'elm_instances', species))
    table = table[table['InstanceLogic'] == 'true positive']
    return pd.DataFrame({'ELMIdentifier': table['ELMIdentifier'], 'uniprot': table['Primary_Acc'],
                         'start': table['Start'].astype(int), 'end': table['End'].astype(int), 'evidence': 'elm'})


def transfer_instances(cfg, human):
    '''Human motif instances on mouse proteins (aligned positions; all residues of the motif must be aligned).'''
    maps = pd.read_parquet(paths(cfg)['intermediate'] / 'position_maps.parquet')
    rows = []
    for m in maps.merge(human, left_on='human_uniprot', right_on='uniprot').itertuples(index=False):
        aligned = dict(zip(m.human_pos, m.mouse_pos))
        mouse = [aligned.get(p) for p in range(m.start, m.end + 1)]
        if None not in mouse:
            rows.append((m.ELMIdentifier, m.mouse_uniprot, min(mouse), max(mouse), 'elm_ortholog'))
    return pd.DataFrame(rows, columns=['ELMIdentifier', 'uniprot', 'start', 'end', 'evidence'])


def motif_exons(cfg, species, motifs):
    mapper = ResidueExonMapper(cfg, species)
    residues = motifs.assign(position=[list(range(s, e + 1)) for s, e in zip(motifs['start'], motifs['end'])]).explode('position')
    residues['position'] = residues['position'].astype(int)
    exons = mapper.map(residues[['ELMIdentifier', 'uniprot', 'position', 'evidence']]).drop(columns='position').drop_duplicates()
    genes = pd.read_parquet(species_dir(cfg, species) / 'genes.parquet', columns=['gene_id', 'entrez']).dropna()
    exons = exons.merge(genes, on='gene_id')
    exons = exons.groupby(['ELMIdentifier', 'transcript_id', 'gene_id', 'exon_id', 'coding_start', 'coding_end', 'entrez', 'uniprot'],
                          as_index=False)['evidence'].agg(lambda s: ','.join(sorted(set(s))))
    return exons.rename(columns={'transcript_id': 'ensembl.transcript', 'gene_id': 'Gene stable ID', 'exon_id': 'Exon stable ID',
                                 'coding_start': 'Genomic coding start', 'coding_end': 'Genomic coding end', 'entrez': 'entrezgene'})


def curated_interactions(cfg, species, uniprot_entrez):
    table = pd.read_csv(input_path(cfg, 'elm_interactions'), sep='\t', dtype=str)
    taxid = str(SPECIES[species]['taxid'])
    # taxa are given as "9606"(Homo sapiens)
    taxon = lambda column: table[column].str.extract(r'(\d+)', expand=False)
    table = table[(taxon('taxonomyElm') == taxid) & (taxon('taxonomyDomain') == taxid) & table['Domain'].str.startswith('PF', na=False)]
    table = table.merge(uniprot_entrez.rename(columns={'uniprot': 'interactorElm', 'entrez': 'gene_1'}), on='interactorElm')
    table = table.merge(uniprot_entrez.rename(columns={'uniprot': 'interactorDomain', 'entrez': 'gene_2'}), on='interactorDomain')
    return pd.DataFrame({'gene_1': table['gene_1'], 'elm': table['Elm'], 'gene_2': table['gene_2'], 'pfam': table['Domain'], 'evidence': 'elm_curated'})


def inferred_interactions(cfg, species, motifs_by_gene, forward):
    '''Motif instances x physical interaction partners with a Pfam domain that binds the motif class.'''
    binding = read_elm_tsv(input_path(cfg, 'elm_domains')).rename(columns={'ELM identifier': 'elm', 'Interaction Domain Id': 'pfam'})
    binding = binding[binding['pfam'].str.startswith('PF', na=False)][['elm', 'pfam']].drop_duplicates()
    ppi = pd.read_parquet(species_dir(cfg, species) / 'ppi.parquet', columns=['entrez_1', 'entrez_2'])
    partners = pd.concat([ppi, ppi.rename(columns={'entrez_1': 'entrez_2', 'entrez_2': 'entrez_1'})], ignore_index=True)
    partners.columns = ['gene_1', 'gene_2']
    table = motifs_by_gene.merge(partners, on='gene_1').merge(binding, on='elm')
    domains = gene_domains(cfg, species, forward).rename(columns={'entrez': 'gene_2'})
    table = table.merge(domains, on=['gene_2', 'pfam'])
    return table[['gene_1', 'elm', 'gene_2', 'pfam']].assign(evidence='elm_inferred')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    args = parser.parse_args()
    cfg = load_config(args.config)
    forward = pfam_forward(cfg)

    human_motifs = instances(cfg, 'Human') if 'Human' in cfg['species'] else None
    human_interactions = None
    for species in cfg['species']:
        uniprot_entrez = read_uniprot_entrez(cfg, species)
        motifs = instances(cfg, species)
        if species == 'Mouse':
            motifs = pd.concat([motifs, transfer_instances(cfg, human_motifs)], ignore_index=True)
        exons = motif_exons(cfg, species, motifs)
        exons.to_parquet(species_dir(cfg, species) / 'elm.parquet')
        log(f'{species}: {len(motifs)} motif instances -> {exons["Exon stable ID"].nunique()} exons of '
            f'{exons["Gene stable ID"].nunique()} genes', exons['evidence'].value_counts().to_dict())

        interactions = [curated_interactions(cfg, species, uniprot_entrez)]
        if species == 'Mouse' and human_interactions is not None:
            orthologs = pd.read_parquet(paths(cfg)['intermediate'] / 'orthologs.parquet', columns=['human_entrez', 'mouse_entrez']).dropna()
            curated = human_interactions[human_interactions['evidence'] == 'elm_curated']
            for side in ('1', '2'):
                curated = curated.merge(orthologs.rename(columns={'human_entrez': f'gene_{side}', 'mouse_entrez': f'mouse_{side}'}), on=f'gene_{side}')
            interactions.append(pd.DataFrame({'gene_1': curated['mouse_1'], 'elm': curated['elm'], 'gene_2': curated['mouse_2'],
                                              'pfam': curated['pfam'], 'evidence': 'elm_curated_ortholog'}))
        if cfg['network']['elm_infer_from_ppi']:
            motifs_by_gene = exons[['entrezgene', 'ELMIdentifier']].drop_duplicates()
            motifs_by_gene.columns = ['gene_1', 'elm']
            interactions.append(inferred_interactions(cfg, species, motifs_by_gene, forward))
        interactions = pd.concat(interactions, ignore_index=True)
        interactions = interactions.groupby(['gene_1', 'elm', 'gene_2', 'pfam'], as_index=False)['evidence'].agg(lambda s: ','.join(sorted(set(s))))
        if species == 'Human':
            human_interactions = interactions.assign(evidence=interactions['evidence'].str.split(',')).explode('evidence')
        interactions.to_parquet(species_dir(cfg, species) / 'elm_interactions.parquet')
        log(f'{species}: {len(interactions)} motif-domain interactions', interactions['evidence'].str.split(',').explode().value_counts().to_dict())
