"""Step 16 (diagnostic): why residue-level interactions of the network shipped with NEASE are missing from the build.

NEASE's residue table ("pdb") lists exons that contain interface residues, with the NCBI gene ID of the partner.
A shipped row counts as kept when a new row has the same gene and partner and overlapping coordinates (this is how
NEASE matches exons to the table; exon IDs are not compared, they differ between Ensembl releases). For every
shipped row that is not kept, the first reason that applies:

    not_coding_in_ensembl          the region does not overlap a coding exon of the gene in the Ensembl release of the build
    pair_kept_other_exons          the gene pair is in the new table, with other exons (interface or exon mapping differs)
    gene_not_in_nease_mapping      the gene of the exon has no NCBI gene ID in NEASE's exon table
    no_uniprot_gene / _partner     no UniProt accession for the gene / the partner (UniProt ID mapping)
    interface_no_identical_transcript   a PDB interface exists, but its UniProt entry of the gene has no identical
                                   Ensembl translation (e.g. one residue differs), so its residues cannot be placed on exons
    interface_not_mapped           a PDB interface exists and could be mapped (unexpected: check the evidence selection)
    below_min_contacts             the proteins touch in a PDB structure, with fewer than pdb.min_interface_contacts
                                   residue pairs
    same_entry_no_contact          both proteins are in a PDB entry that was analysed, but they do not touch there
                                   (the biological assembly is used: contacts between copies in the crystal are not
                                   counted; or the chains are far apart)
    entry_excluded:<reason>        both proteins are only in PDB entries that were excluded (method, resolution,
                                   too large, single chain in the assembly, missing, error)
    human_interface_not_transferred   (mouse) the human orthologs have a PDB interface that was not transferred
                                   (no one-to-one orthologs, protein identity < orthologs.min_identity, or the
                                   interface residues do not align)
    related_structure_only         (step 17 run) only a structure of closely related proteins (other species,
                                   paralogs) has the interface: the tier pdb_homolog is not selected (or its
                                   residues are on other exons)
    no_pdb_entry_with_both         no PDB entry contains both proteins (SIFTS): the shipped row came from a
                                   structure of more distant relatives (homology), another isoform or an older mapping

Rows whose gene pair is in NEASE's shipped PPI network are reported separately: they are the losses that matter
most (for mouse, most shipped rows are not supported by any known interaction).

Writes <work_dir>/output/<name>/losses.md and losses_<species>.tsv.gz (one row per lost shipped row)

usage: python 16_diagnose_losses.py --config config.yaml [--name nease_data]
"""
import argparse
import pickle
from collections import Counter

import pandas as pd

from common import load_config, log, nease_file, overlapping, paths, placeable_proteins, read_uniprot_entrez, species_dir


def entrez_text(values):
    return values.astype(str).str.replace(r'\.0$', '', regex=True)


def gene_to_entrez(cfg, species):
    '''Gene of a residue row -> NCBI gene ID, as NEASE does it (Ensembl gene ID, otherwise the symbol).'''
    mapping = pd.read_pickle(nease_file(cfg, species, 'exon_domains'), compression='gzip')
    ensembl = dict(zip(mapping['Gene stable ID'], entrez_text(mapping['NCBI gene ID'])))
    symbol = dict(zip(mapping['Gene name'], entrez_text(mapping['NCBI gene ID'])))
    names = dict(zip(entrez_text(mapping['NCBI gene ID']), mapping['Gene name']))
    return lambda gene, sym: ensembl.get(gene, symbol.get(sym)), names


def residue_rows(table, to_entrez):
    table = table.copy()
    table['partner'] = entrez_text(table['entrezgene'])
    table['gene'] = [to_entrez(g, s) for g, s in zip(table['Gene stable ID'], table['symbol'])]
    return table


class Structures:
    '''What the PDB steps saw for protein pairs: merged interfaces, all chain contacts, entries and their status.'''

    def __init__(self, cfg):
        folder = paths(cfg)['intermediate'] / 'pdb'
        interfaces = pd.read_parquet(folder / 'interfaces_pdb.parquet', columns=['species', 'uniprot_1', 'uniprot_2'])
        self.interfaces = {s: set(zip(t['uniprot_1'], t['uniprot_2'])) for s, t in interfaces.groupby('species')}
        path = folder / 'interfaces_homolog.parquet'
        related = pd.read_parquet(path, columns=['species', 'uniprot_1', 'uniprot_2']) if path.exists() else pd.DataFrame(columns=['species', 'uniprot_1', 'uniprot_2'])
        self.related = {s: set(zip(t['uniprot_1'], t['uniprot_2'])) for s, t in related.groupby('species')}
        # only proteins of the species of the build (the PDB tables cover all species)
        proteins = set()
        for species in cfg['species']:
            proteins |= set(read_uniprot_entrez(cfg, species)['uniprot'])
        parts = []
        for shard in range(cfg['pdb']['shards']):
            path = folder / 'contacts' / f'shard_{shard:05d}.interfaces.parquet'
            if path.exists():
                part = pd.read_parquet(path, columns=['uniprot_1', 'uniprot_2', 'n_pairs'])
                part = part[part['uniprot_1'].isin(proteins) & part['uniprot_2'].isin(proteins)]
                low, high = part[['uniprot_1', 'uniprot_2']].min(axis=1), part[['uniprot_1', 'uniprot_2']].max(axis=1)
                parts.append(pd.DataFrame({'u1': low, 'u2': high, 'n': part['n_pairs']}).groupby(['u1', 'u2'])['n'].max())
        contacts = pd.concat(parts).groupby(level=[0, 1]).max() if parts else pd.Series(dtype='int64')
        self.max_contacts = contacts.to_dict()
        chains = pd.read_parquet(folder / 'uniprot_map.parquet', columns=['pdb', 'uniprot']).drop_duplicates()
        chains = chains[chains['uniprot'].isin(proteins)]
        self.entries = chains.groupby('uniprot')['pdb'].agg(set).to_dict()
        status = pd.read_parquet(folder / 'entries_log.parquet', columns=['pdb', 'status'])
        self.status = dict(zip(status['pdb'], status['status'].str.split(':').str[0]))
        log(f'PDB: {sum(map(len, self.interfaces.values()))} interfaces, {len(self.max_contacts)} protein pairs in contact, '
            f'{len(self.entries)} proteins in {len(self.status)} entries')

    def reason(self, species, pairs, identical):
        '''Why there is no usable PDB interface for any of the protein pairs (UniProt accessions of gene, partner).'''
        keys = [(min(a, b), max(a, b)) for a, b in pairs]
        found = [(a, b) for (a, b), k in zip(pairs, keys) if k in self.interfaces.get(species, ())]
        if found:
            # the residues of the gene are placed on exons only if the UniProt entry of the interface itself has an
            # identical Ensembl translation (other entries of the gene do not help)
            return 'interface_not_mapped' if any(a in identical for a, _ in found) else 'interface_no_identical_transcript'
        if any(k in self.related.get(species, ()) for k in keys):
            return 'related_structure_only'
        if any(self.max_contacts.get(k, 0) > 0 for k in keys):
            return 'below_min_contacts'
        shared = set()
        for a, b in pairs:
            shared |= self.entries.get(a, set()) & self.entries.get(b, set())
        if shared:
            statuses = Counter(self.status.get(e, 'not analysed') for e in shared)
            if statuses.get('ok'):
                return 'same_entry_no_contact'
            return f'entry_excluded:{statuses.most_common(1)[0][0]}'
        return None


def report(cfg, species, folder, structures):
    to_entrez, names = gene_to_entrez(cfg, species)
    with open(nease_file(cfg, species, 'ppi'), 'rb') as f:
        shipped_ppi = {frozenset(map(str, e)) for e in pickle.load(f).edges()}
    shipped = residue_rows(pd.read_pickle(nease_file(cfg, species, 'pdb')), to_entrez)
    new = residue_rows(pd.read_parquet(folder / 'details' / species / 'pdb.parquet'), to_entrez)
    new_pairs = {frozenset((g, p)) for g, p in zip(new['gene'], new['partner']) if g}
    shipped['kept'] = overlapping(shipped, new, ['Gene stable ID', 'partner'])
    shipped['in_shipped_ppi'] = [bool(g) and frozenset((g, p)) in shipped_ppi for g, p in zip(shipped['gene'], shipped['partner'])]

    exons = pd.read_parquet(species_dir(cfg, species) / 'exon_cds.parquet', columns=['gene_id', 'coding_start', 'coding_end'])
    shipped['coding'] = overlapping(shipped, exons.rename(columns={'gene_id': 'Gene stable ID'}), ['Gene stable ID'],
                                    other_start='coding_start', other_end='coding_end')
    identical = set(placeable_proteins(cfg, species)['uniprot'])  # identical to a translation, or aligned to one
    uniprots = read_uniprot_entrez(cfg, species).groupby('entrez')['uniprot'].agg(list).to_dict()
    human_uniprots, orthologs = None, {}
    if species == 'Mouse':
        human_uniprots = read_uniprot_entrez(cfg, 'Human').groupby('entrez')['uniprot'].agg(list).to_dict()
        table = pd.read_parquet(paths(cfg)['intermediate'] / 'orthologs.parquet', columns=['human_entrez', 'mouse_entrez']).dropna()
        orthologs = table.groupby('mouse_entrez')['human_entrez'].agg(list).to_dict()

    pair_reason = {}

    def classify_pair(gene, partner):
        if (gene, partner) in pair_reason:
            return pair_reason[(gene, partner)]
        if not gene:
            reason = 'gene_not_in_nease_mapping'
        elif gene not in uniprots:
            reason = 'no_uniprot_gene'
        elif partner not in uniprots:
            reason = 'no_uniprot_partner'
        else:
            pairs = [(a, b) for a in uniprots[gene] for b in uniprots[partner]]
            reason = structures.reason(species, pairs, identical)
            if reason is None and species == 'Mouse':
                human_pairs = [(a, b) for hg in orthologs.get(gene, []) for hp in orthologs.get(partner, [])
                               for a in human_uniprots.get(hg, []) for b in human_uniprots.get(hp, [])]
                if any((min(a, b), max(a, b)) in structures.interfaces.get('Human', ()) for a, b in human_pairs):
                    reason = 'human_interface_not_transferred'
            reason = reason or 'no_pdb_entry_with_both'
        pair_reason[(gene, partner)] = reason
        return reason

    lost = shipped[~shipped['kept']].copy()
    reasons = []
    for coding, gene, partner in zip(lost['coding'], lost['gene'], lost['partner']):
        if not coding:
            reasons.append('not_coding_in_ensembl')
        elif gene and frozenset((gene, partner)) in new_pairs:
            reasons.append('pair_kept_other_exons')
        else:
            reasons.append(classify_pair(gene, partner))
    lost['reason'] = reasons
    lost['partner_symbol'] = lost['partner'].map(names)
    lost[['symbol', 'Gene stable ID', 'Exon with residue ID', 'gene', 'partner', 'partner_symbol', 'in_shipped_ppi', 'reason']].to_csv(
        folder / f'losses_{species.lower()}.tsv.gz', sep='\t', index=False)

    n, n_ppi = len(shipped), int(shipped['in_shipped_ppi'].sum())
    kept, kept_ppi = int(shipped['kept'].sum()), int((shipped['kept'] & shipped['in_shipped_ppi']).sum())
    lines = [f'## {species}\n',
             f'- Shipped residue rows (exon, partner): {n}; kept {kept} ({100 * kept / max(n, 1):.1f}%), lost {n - kept}',
             f'- Of these, rows whose gene pair is in the shipped PPI network: {n_ppi}; kept {kept_ppi} '
             f'({100 * kept_ppi / max(n_ppi, 1):.1f}%)',
             f'- New residue rows: {len(new)}\n',
             '| reason | lost rows | lost rows in shipped PPI | lost gene pairs | examples (gene - partner) |', '|---|---|---|---|---|']
    pairs = lost.drop_duplicates(['gene', 'partner', 'reason'])
    for reason, count in lost['reason'].value_counts().items():
        rows = lost[lost['reason'] == reason]
        examples = rows.groupby(['symbol', 'partner_symbol']).size().sort_values(ascending=False).head(4)
        lines.append(f'| {reason} | {count} | {int(rows["in_shipped_ppi"].sum())} | {int((pairs["reason"] == reason).sum())} | '
                     + ', '.join(f'{g}-{p}' for g, p in examples.index) + ' |')
    log(f'{species}: kept {kept}/{n}, lost by reason', lost['reason'].value_counts().to_dict())
    return lines


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    parser.add_argument('--name', default='nease_data', help='output folder of 13_assemble.py')
    args = parser.parse_args()
    cfg = load_config(args.config)
    folder = paths(cfg)['output'] / args.name
    structures = Structures(cfg)
    lines = [f'# Shipped residue-level interactions missing from {folder}\n',
             'See 16_diagnose_losses.py for the reasons. Per-row tables: losses_<species>.tsv.gz\n']
    for species in cfg['species']:
        lines += report(cfg, species, folder, structures) + ['']
    (folder / 'losses.md').write_text('\n'.join(lines))
    print('\n'.join(lines))
    log(f'written: {folder / "losses.md"}')
