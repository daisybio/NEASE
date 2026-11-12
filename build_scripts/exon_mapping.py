"""Map residues of UniProt sequences to the exons that encode them (tables of 02_ensembl.py)."""
import pandas as pd

from common import species_dir


class ResidueExonMapper:
    def __init__(self, cfg, species):
        folder = species_dir(cfg, species)
        self.proteins = pd.read_parquet(folder / 'proteins.parquet', columns=['uniprot', 'transcript_id', 'gene_id', 'symbol'])
        self.identical = set(self.proteins['uniprot'])
        # UniProt entries without an identical translation, aligned to transcripts (step 2; absent in older builds)
        path = folder / 'uniprot_maps.parquet'
        self.maps = pd.read_parquet(path, columns=['uniprot', 'transcript_id', 'gene_id', 'symbol', 'uni_pos', 'ens_pos']) if path.exists() else None
        exons = pd.read_parquet(folder / 'exon_cds.parquet')
        self.exons = exons.sort_values('cds_start', kind='stable').reset_index(drop=True)

    def map(self, residues):
        '''
        residues: DataFrame with the columns uniprot and position (1-based) plus any other (hashable) columns.
        Returns one row per input row and exon that encodes the residue (in any transcript identical to the
        UniProt sequence, or aligned to it; a codon split by an intron maps to both exons), with the added columns
        transcript_id, gene_id, symbol, exon_id, coding_start, coding_end.
        '''
        table = residues.merge(self.proteins, on='uniprot').assign(_tpos=lambda t: t['position'])
        if self.maps is not None:
            other = residues[~residues['uniprot'].isin(self.identical)]
            maps = self.maps[self.maps['uniprot'].isin(set(other['uniprot']))]
            if len(maps):
                # position in the UniProt sequence -> position in the aligned translation
                maps = maps.explode(['uni_pos', 'ens_pos']).rename(columns={'uni_pos': 'position', 'ens_pos': '_tpos'})
                maps = maps.astype({'position': 'int64', '_tpos': 'int64'})
                other = other.astype({'position': 'int64'}).merge(maps, on=['uniprot', 'position'])
                table = pd.concat([table, other], ignore_index=True)
        hits = []
        # first and last nucleotide of the codon
        for offset in (2, 0):
            left = table.assign(nt=3 * table['_tpos'] - offset).sort_values('nt', kind='stable')
            right = self.exons[['transcript_id', 'exon_id', 'coding_start', 'coding_end', 'cds_start', 'cds_end']]
            merged = pd.merge_asof(left, right, left_on='nt', right_on='cds_start', by='transcript_id', direction='backward')
            hits.append(merged[merged['nt'] <= merged['cds_end']])
        out = pd.concat(hits, ignore_index=True).drop(columns=['nt', 'cds_start', 'cds_end', '_tpos'])
        return out.drop_duplicates().reset_index(drop=True)
