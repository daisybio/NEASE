"""Step 2: exon and protein tables of the Ensembl release (per species).

Interface residues and motifs are given as positions on UniProt sequences. They are mapped to exons through
the Ensembl transcripts whose protein is identical to the UniProt sequence (100% identity in the Ensembl
UniProt cross-references, protein length consistent with the coding sequence). UniProt entries without such a
transcript (often one residue differs, or UniProt's canonical isoform is not an Ensembl translation) are aligned
to the translations they are linked to (Ensembl cross-references, otherwise the canonical transcript of the gene)
and kept when >= uniprot.min_identity % of the UniProt sequence is aligned to identical residues.

Writes to <work_dir>/intermediate/<species>/:
    exon_cds.parquet   coding part of every exon of protein-coding transcripts: genomic coding start/end and
                       position in the coding sequence (cds_start, cds_end; nucleotides, 1-based)
    proteins.parquet   UniProt accession -> identical Ensembl transcripts (gene, symbol, canonical, sequence)
    genes.parquet      Ensembl gene -> symbol, NCBI gene IDs
    uniprot_maps.parquet   UniProt entries without an identical translation -> aligned transcripts: uniprot,
                       transcript_id, gene_id, symbol, canonical, swissprot, identity, uni_pos, ens_pos, sequence

usage: python 02_ensembl.py --config config.yaml [--species Human]
"""
import argparse
import gzip
from multiprocessing import Pool

import pandas as pd

from common import align_positions, input_path, load_config, log, read_uniprot_entrez, species_dir, threads, uniprot_sequences

ATTRIBUTES = {'gene_id': r'gene_id "([^"]+)"', 'transcript_id': r'transcript_id "([^"]+)"',
              'exon_number': r'exon_number "(\d+)"', 'exon_id': r'exon_id "([^"]+)"', 'gene_name': r'gene_name "([^"]+)"',
              'gene_biotype': r'gene_biotype "([^"]+)"', 'transcript_biotype': r'transcript_biotype "([^"]+)"',
              'protein_id': r'protein_id "([^"]+)"'}
COLUMNS = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']


def read_gtf(path):
    parts = []
    for chunk in pd.read_csv(path, sep='\t', comment='#', header=None, names=COLUMNS, dtype={'chromosome': str},
                             chunksize=500000):
        chunk = chunk[chunk['feature'].isin(['gene', 'transcript', 'exon', 'CDS'])]
        for name, pattern in ATTRIBUTES.items():
            chunk[name] = chunk['attributes'].str.extract(pattern, expand=False)
        chunk['canonical'] = chunk['attributes'].str.contains('tag "Ensembl_canonical"', regex=False)
        chunk['start_nf'] = chunk['attributes'].str.contains('tag "cds_start_NF"', regex=False)
        parts.append(chunk.drop(columns=['attributes', 'source', 'score']))
    return pd.concat(parts, ignore_index=True)


def read_peptides(path):
    '''Ensembl protein FASTA -> transcript ID, protein ID, sequence (IDs without versions).'''
    rows, header, seq = [], None, []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if header:
                    rows.append((*header, ''.join(seq)))
                fields = dict(x.split(':', 1) for x in line[1:].split()[1:] if ':' in x)
                header, seq = (fields['transcript'].split('.')[0], line[1:].split()[0].split('.')[0]), []
            else:
                seq.append(line.strip())
    if header:
        rows.append((*header, ''.join(seq)))
    return pd.DataFrame(rows, columns=['transcript_id', 'protein_id', 'sequence'])


def build(cfg, species):
    out = species_dir(cfg, species)
    gtf = read_gtf(input_path(cfg, 'ensembl_gtf', species))
    log(f'{species}: {len(gtf)} GTF records')

    transcripts = gtf[(gtf['feature'] == 'transcript') & (gtf['transcript_biotype'] == 'protein_coding')]
    transcripts = transcripts[['transcript_id', 'gene_id', 'gene_name', 'canonical', 'start_nf']]

    exons = gtf.loc[gtf['feature'] == 'exon', ['transcript_id', 'exon_number', 'exon_id']]
    cds = gtf.loc[gtf['feature'] == 'CDS', ['transcript_id', 'exon_number', 'chromosome', 'strand', 'start', 'end', 'frame']]
    cds = cds[cds['transcript_id'].isin(transcripts['transcript_id'])].merge(exons, on=['transcript_id', 'exon_number'])
    cds['exon_rank'] = cds['exon_number'].astype(int)
    cds = cds.sort_values(['transcript_id', 'exon_rank'], kind='stable')
    length = cds['end'] - cds['start'] + 1
    cds['cds_end'] = length.groupby(cds['transcript_id']).cumsum()
    cds['cds_start'] = cds['cds_end'] - length + 1

    # coding sequences that start in frame (positions of residues are counted from the first codon)
    first = cds.groupby('transcript_id').first()
    complete = set(first.index[first['frame'].astype(str) == '0']) - set(transcripts.loc[transcripts['start_nf'], 'transcript_id'])
    cds = cds[cds['transcript_id'].isin(complete)].merge(transcripts[['transcript_id', 'gene_id']], on='transcript_id')
    exon_cds = cds.rename(columns={'start': 'coding_start', 'end': 'coding_end'})[
        ['transcript_id', 'gene_id', 'exon_id', 'exon_rank', 'chromosome', 'strand', 'coding_start', 'coding_end', 'cds_start', 'cds_end']]
    exon_cds.reset_index(drop=True).to_parquet(out / 'exon_cds.parquet')
    log(f'{species}: {exon_cds["transcript_id"].nunique()} coding transcripts, {exon_cds["exon_id"].nunique()} exons')

    # proteins: UniProt accessions with identical Ensembl translations
    peptides = read_peptides(input_path(cfg, 'ensembl_pep', species))
    cds_length = cds.groupby('transcript_id')['cds_end'].max()
    peptides['cds_length'] = peptides['transcript_id'].map(cds_length)
    n = peptides['sequence'].str.len()
    peptides = peptides[(peptides['cds_length'] == 3 * n) | (peptides['cds_length'] == 3 * n + 3)]

    all_xrefs = pd.read_csv(input_path(cfg, 'ensembl_uniprot', species), sep='\t', dtype=str)
    xrefs = all_xrefs[(all_xrefs['source_identity'] == '100') & (all_xrefs['xref_identity'] == '100')]
    proteins = xrefs.rename(columns={'xref': 'uniprot', 'transcript_stable_id': 'transcript_id'})[['uniprot', 'transcript_id', 'db_name']]
    proteins['swissprot'] = proteins['db_name'].str.contains('SWISSPROT|isoform', case=False)
    proteins = proteins.merge(peptides[['transcript_id', 'protein_id', 'sequence']], on='transcript_id')
    proteins = proteins.merge(transcripts[['transcript_id', 'gene_id', 'gene_name', 'canonical']], on='transcript_id')
    proteins = proteins.rename(columns={'gene_name': 'symbol'})
    proteins['length'] = proteins['sequence'].str.len()
    proteins = proteins[['uniprot', 'transcript_id', 'protein_id', 'gene_id', 'symbol', 'canonical', 'swissprot', 'length', 'sequence']]
    proteins = pd.concat([proteins, isoform_proteins(cfg, species, all_xrefs, peptides, transcripts)], ignore_index=True)
    proteins.drop_duplicates(['uniprot', 'transcript_id']).reset_index(drop=True).to_parquet(out / 'proteins.parquet')
    log(f'{species}: {proteins["uniprot"].nunique()} UniProt accessions with identical transcripts '
        f'({proteins["transcript_id"].nunique()} transcripts, {proteins["gene_id"].nunique()} genes)')

    # genes: symbol and NCBI gene IDs
    genes = gtf.loc[gtf['feature'] == 'gene', ['gene_id', 'gene_name', 'gene_biotype']].rename(columns={'gene_name': 'symbol', 'gene_biotype': 'biotype'})
    entrez = pd.read_csv(input_path(cfg, 'ensembl_entrez', species), sep='\t', dtype=str)
    entrez = entrez[entrez['db_name'] == 'EntrezGene'][['gene_stable_id', 'xref']].drop_duplicates()
    entrez.columns = ['gene_id', 'entrez']
    genes = genes.merge(entrez, on='gene_id', how='left')
    genes.reset_index(drop=True).to_parquet(out / 'genes.parquet')
    log(f'{species}: {genes["gene_id"].nunique()} genes, {genes["entrez"].notna().sum()} with NCBI gene IDs')

    maps = uniprot_maps(cfg, species, all_xrefs, peptides, transcripts, genes, set(proteins['uniprot']))
    maps.to_parquet(out / 'uniprot_maps.parquet')


def isoform_proteins(cfg, species, xrefs, peptides, transcripts):
    '''
    UniProt isoform entries (Q674X7-2) identical to the Ensembl translation they are linked to. Ensembl gives these
    links without identity values, so the sequences (Swiss-Prot isoforms; '-1' is the canonical sequence) are compared.
    '''
    links = xrefs.loc[xrefs['db_name'] == 'Uniprot_isoform', ['xref', 'transcript_stable_id']].drop_duplicates()
    links.columns = ['uniprot', 'transcript_id']
    links = links.merge(peptides[['transcript_id', 'protein_id', 'sequence']], on='transcript_id')
    canonical = links['uniprot'].str.split('-').str[0]
    sequences, _ = uniprot_sequences(cfg, set(links['uniprot']) | set(canonical))
    isoform_seq = [sequences.get(u) or (sequences.get(c) if u.endswith('-1') else None) for u, c in zip(links['uniprot'], canonical)]
    links = links[[s is not None and s == q for s, q in zip(isoform_seq, links['sequence'])]]
    links = links.merge(transcripts[['transcript_id', 'gene_id', 'gene_name', 'canonical']], on='transcript_id').rename(columns={'gene_name': 'symbol'})
    links = links.assign(swissprot=True, length=links['sequence'].str.len())
    log(f'{species}: {links["uniprot"].nunique()} UniProt isoform entries identical to their linked translations')
    return links[['uniprot', 'transcript_id', 'protein_id', 'gene_id', 'symbol', 'canonical', 'swissprot', 'length', 'sequence']]


def align_job(job):
    uniprot, transcript, uniprot_seq, ensembl_seq = job
    uni_pos, ens_pos, identical = align_positions(uniprot_seq, ensembl_seq)
    return uniprot, transcript, identical, uni_pos, ens_pos


def uniprot_maps(cfg, species, xrefs, peptides, transcripts, genes, identical):
    '''UniProt entries without an identical Ensembl translation, aligned to the translations they are linked to.'''
    min_identity = (cfg.get('uniprot') or {}).get('min_identity', 95)
    linked = xrefs.rename(columns={'xref': 'uniprot', 'transcript_stable_id': 'transcript_id'})[['uniprot', 'transcript_id']]
    linked = linked[~linked['uniprot'].isin(identical) & ~linked['uniprot'].str.contains('-')]
    # also the canonical transcript of the gene (UniProt ID mapping): the cross-referenced transcripts can be
    # incomplete (no usable translation)
    ids = read_uniprot_entrez(cfg, species)
    ids = ids[~ids['uniprot'].isin(identical) & ~ids['uniprot'].str.contains('-')]
    canonical = transcripts.loc[transcripts['canonical'], ['transcript_id', 'gene_id']].merge(genes[['gene_id', 'entrez']].dropna(), on='gene_id')
    candidates = pd.concat([linked, ids.merge(canonical, on='entrez')[['uniprot', 'transcript_id']]]).drop_duplicates()
    candidates = candidates.merge(peptides[['transcript_id', 'sequence']], on='transcript_id')
    sequences, swissprot = uniprot_sequences(cfg, set(candidates['uniprot']))
    candidates['uniprot_seq'] = candidates['uniprot'].map(sequences)
    candidates = candidates.dropna(subset=['uniprot_seq'])
    # a translation shorter than min_identity % of the UniProt sequence cannot reach min_identity
    candidates = candidates[candidates['sequence'].str.len() >= candidates['uniprot_seq'].str.len() * min_identity / 100]
    jobs = list(candidates[['uniprot', 'transcript_id', 'uniprot_seq', 'sequence']].itertuples(index=False, name=None))
    log(f'{species}: aligning {len(jobs)} UniProt entries without an identical translation to {candidates["transcript_id"].nunique()} transcripts')
    with Pool(threads()) as pool:
        results = pool.map(align_job, jobs, chunksize=20)
    rows = []
    for uniprot, transcript, n_identical, uni_pos, ens_pos in results:
        identity = 100 * n_identical / len(sequences[uniprot])
        if identity >= min_identity:
            rows.append((uniprot, transcript, round(identity, 1), uni_pos, ens_pos, sequences[uniprot]))
    maps = pd.DataFrame(rows, columns=['uniprot', 'transcript_id', 'identity', 'uni_pos', 'ens_pos', 'sequence'])
    maps = maps.merge(transcripts[['transcript_id', 'gene_id', 'gene_name', 'canonical']], on='transcript_id').rename(columns={'gene_name': 'symbol'})
    maps['swissprot'] = maps['uniprot'].isin(swissprot)
    log(f'{species}: {maps["uniprot"].nunique()} UniProt entries without an identical translation aligned to '
        f'{maps["transcript_id"].nunique()} transcripts (>= {min_identity}% identical; {maps["uniprot"][maps["swissprot"]].nunique()} Swiss-Prot)')
    return maps[['uniprot', 'transcript_id', 'gene_id', 'symbol', 'canonical', 'swissprot', 'identity', 'uni_pos', 'ens_pos', 'sequence']].reset_index(drop=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    parser.add_argument('--species', nargs='+')
    args = parser.parse_args()
    cfg = load_config(args.config)
    for species in args.species or cfg['species']:
        build(cfg, species)
