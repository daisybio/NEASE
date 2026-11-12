"""Build NEASE's exon -> Pfam domain table (nease/data/database/Human or Mouse) from Ensembl BioMart.

Every row is an exon of a protein-coding transcript whose coding part overlaps a Pfam
domain of that transcript's protein. The overlap is computed on protein positions:
CDS start / 3 <= Pfam end and CDS end / 3 >= Pfam start (same rule as the original tables).
Genes on the primary chromosomes with an NCBI gene ID are kept; a gene with several
NCBI gene IDs gets one row per ID.

usage: python scripts/build_exon_domain_table.py Human|Mouse <output> [--mart https://jun2026.archive.ensembl.org]

The default mart is the Ensembl 116 archive. With https://feb2023.archive.ensembl.org
(Ensembl 109) the script reproduces the mouse table shipped with NEASE <= 1.3.1 row for row.
The table is saved as a gzip-compressed pickle (read it with pd.read_pickle(path, compression='gzip')).
"""
import argparse
import io
import time

import pandas as pd
import requests

DATASETS = {'Human': 'hsapiens_gene_ensembl', 'Mouse': 'mmusculus_gene_ensembl'}
CHROMOSOMES = {'Human': [str(c) for c in range(1, 23)] + ['X', 'Y'],
               'Mouse': [str(c) for c in range(1, 20)] + ['X', 'Y']}

STRUCTURE = {'chromosome_name': 'Chromosome/scaffold name', 'strand': 'Strand',
             'ensembl_exon_id': 'Exon stable ID', 'rank': 'Exon rank in transcript',
             'genomic_coding_start': 'Genomic coding start', 'genomic_coding_end': 'Genomic coding end',
             'cds_start': 'CDS start', 'cds_end': 'CDS end',
             'ensembl_transcript_id': 'Transcript stable ID', 'ensembl_gene_id': 'Gene stable ID'}
DOMAINS = {'ensembl_transcript_id': 'Transcript stable ID', 'pfam': 'Pfam ID',
           'pfam_start': 'Pfam start', 'pfam_end': 'Pfam end'}
GENES = {'ensembl_gene_id': 'Gene stable ID', 'external_gene_name': 'Gene name',
         'description': 'Gene description', 'entrezgene_id': 'NCBI gene ID'}

COLUMNS = ['Chromosome/scaffold name', 'Strand', 'Exon stable ID', 'Exon rank in transcript',
           'Genomic coding start', 'Genomic coding end', 'CDS start', 'CDS end', 'Transcript stable ID',
           'Pfam ID', 'Pfam start', 'Pfam end', 'Gene stable ID', 'Gene name', 'Gene description', 'NCBI gene ID']


def biomart(mart, dataset, attributes, chromosome):
    query = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
             '<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" completionStamp="1" datasetConfigVersion="0.6">'
             f'<Dataset name="{dataset}" interface="default">'
             f'<Filter name="chromosome_name" value="{chromosome}"/><Filter name="transcript_biotype" value="protein_coding"/>'
             + ''.join(f'<Attribute name="{a}"/>' for a in attributes) + '</Dataset></Query>')
    for attempt in range(10):
        try:
            r = requests.get(f'{mart}/biomart/martservice', params={'query': query}, timeout=1800)
            # BioMart reports errors (and an unavailable service) with HTTP 200; complete answers end with [success]
            if r.ok and r.text.rstrip().endswith('[success]'):
                text = r.text.rstrip()[:-len('[success]')]
                return pd.read_csv(io.StringIO(text), sep='\t', header=None, names=list(attributes), dtype=str)
            error = r.text[:200]
        except requests.RequestException as e:
            error = str(e)
        print(f'  BioMart request failed ({error.strip()}), retrying')
        time.sleep(30 * (attempt + 1))
    raise RuntimeError(f'BioMart query failed for chromosome {chromosome}')


def build(organism, mart):
    dataset = DATASETS[organism]
    structure, domains, genes = [], [], []
    for chromosome in CHROMOSOMES[organism]:
        print(f'chromosome {chromosome}')
        structure.append(biomart(mart, dataset, STRUCTURE, chromosome).rename(columns=STRUCTURE))
        domains.append(biomart(mart, dataset, DOMAINS, chromosome).rename(columns=DOMAINS))
        genes.append(biomart(mart, dataset, GENES, chromosome).rename(columns=GENES))
    structure = pd.concat(structure).dropna(subset=['Genomic coding start', 'CDS start'])
    domains = pd.concat(domains).dropna(subset=['Pfam ID'])
    genes = pd.concat(genes).dropna(subset=['NCBI gene ID']).drop_duplicates()

    table = structure.merge(domains, on='Transcript stable ID')
    for column in ['Strand', 'Exon rank in transcript', 'Genomic coding start', 'Genomic coding end']:
        table[column] = table[column].astype('int64')
    for column in ['CDS start', 'CDS end', 'Pfam start', 'Pfam end']:
        table[column] = table[column].astype('float64')
    overlap = (table['CDS start'] / 3 <= table['Pfam end']) & (table['CDS end'] / 3 >= table['Pfam start'])
    table = table[overlap].merge(genes, on='Gene stable ID')
    table['NCBI gene ID'] = table['NCBI gene ID'].astype('int64')
    table['Gene description'] = table['Gene description'].fillna('')
    table = table[COLUMNS].drop_duplicates()
    # object columns instead of pandas 3 string columns, so that the pickle also loads with pandas 2
    table = table.astype({c: object for c in table.columns if pd.api.types.is_string_dtype(table[c])})
    table.columns = table.columns.astype(object)
    return table.sort_values(COLUMNS[:2] + ['Gene stable ID', 'Transcript stable ID', 'Exon rank in transcript', 'Pfam start', 'Pfam ID', 'NCBI gene ID'],
                             kind='stable').reset_index(drop=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('organism', choices=DATASETS)
    parser.add_argument('output')
    parser.add_argument('--mart', default='https://jun2026.archive.ensembl.org')
    args = parser.parse_args()
    table = build(args.organism, args.mart)
    table.to_pickle(args.output, compression='gzip')
    print(f'{len(table)} rows, {table["Gene stable ID"].nunique()} genes, '
          f'{table["Transcript stable ID"].nunique()} transcripts, {table["Pfam ID"].nunique()} Pfam domains -> {args.output}')
