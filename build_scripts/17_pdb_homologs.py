"""Step 17 (optional): residue-level interfaces from PDB structures of closely related proteins.

Steps 5-6 use a structure for a human/mouse protein pair only if its chains are those proteins themselves (SIFTS).
Many complexes are solved for the same proteins of another species (rat, cow, rabbit actin, Xenopus histones, ...)
or for a near-identical paralog. This step maps every chain protein of the PDB interfaces (step 5) to the human
and mouse proteins it is closely related to and transfers the interfaces:

  1. sequences of the chain proteins (UniProt: Swiss-Prot of all species and the reference proteomes)
  2. candidate pairs (chain protein, human/mouse protein) that share >= homologs.min_shared_words identical
     10-residue words (low-complexity words are ignored), then local alignments of the parts of the chain protein
     that are in the structures (so that repeats, e.g. the nine ubiquitins of polyubiquitin, align at the right
     copy): a part is kept with >= homologs.min_identity % identity over >= homologs.min_aligned aligned residues
  3. every interface (>= pdb.min_interface_contacts residue pairs) whose two chains map to two proteins of the same
     species, other than the chains' own proteins, with >= homologs.min_mapped of its residue pairs mapped; the
     residue pairs of all structures of a protein pair are merged

Evidence tiers (network.evidence / interface_evidence): pdb_homolog_pair (domain level), pdb_homolog (residue level).
Interfaces of structures of the proteins themselves stay in the tier pdb. Like the other structure tiers, a
transferred interface counts as a physical interaction when network.structure_implies_ppi is true.

Writes to <work_dir>/intermediate/pdb/:
    homolog_maps.parquet        query (chain protein), target (human/mouse protein), species, start, end (aligned
                                range of the query), identity, aligned, q_pos, t_pos (aligned positions)
    interfaces_homolog.parquet  species, uniprot_1, uniprot_2, pos_1, pos_2, n_pairs, n_entries, pdb, min_identity

Run after steps 2 and 5 and before step 13. usage: python 17_pdb_homologs.py --config config.yaml
"""
import argparse
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
import pandas as pd

from common import align_positions, load_config, log, paths, placeable_proteins, read_uniprot_entrez, threads, uniprot_sequences

WORD = 10
ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'
LOOKUP = np.full(256, len(ALPHABET), dtype=np.int64)
for _i, _a in enumerate(ALPHABET):
    LOOKUP[ord(_a)] = _i
POWERS = (len(ALPHABET) + 1) ** np.arange(WORD - 1, -1, -1, dtype=np.int64)
MIN_DISTINCT = 4             # words with fewer distinct residues are low-complexity (QQQQQ..., EEKEEKEE...)
MAX_TARGETS_PER_WORD = 1000  # guard against huge repeat families (ubiquitin, histones, actins have up to ~100)
SEGMENT_MARGIN = 20          # residues around the interface residues of a chain protein that are aligned
SPAN_STEP = 25               # aligned ranges are rounded to this (similar structures share an alignment)


def settings(cfg):
    s = dict(enabled=False, min_identity=95, min_aligned=50, min_mapped=0.8, min_shared_words=5)
    s.update(cfg.get('homologs') or {})
    return s


def words(sequence, step=1):
    '''Codes of the 10-residue words of a sequence (every step-th word; words with unusual residues are skipped).'''
    if len(sequence) < WORD:
        return np.empty(0, dtype=np.int64)
    codes = LOOKUP[np.frombuffer(sequence.encode('ascii', 'replace'), dtype=np.uint8)]
    windows = np.lib.stride_tricks.sliding_window_view(codes, WORD)[::step]
    windows = windows[(windows < len(ALPHABET)).all(axis=1)]
    ordered = np.sort(windows, axis=1)
    windows = windows[(np.diff(ordered, axis=1) != 0).sum(axis=1) + 1 >= MIN_DISTINCT]
    return np.unique(windows @ POWERS)


def word_table(sequences, step):
    ids, codes = [], []
    for i, sequence in enumerate(sequences):
        w = words(sequence, step)
        codes.append(w)
        ids.append(np.full(len(w), i, dtype=np.int32))
    return pd.DataFrame({'word': np.concatenate(codes) if codes else np.empty(0, dtype=np.int64),
                         'id': np.concatenate(ids) if ids else np.empty(0, dtype=np.int32)})


def candidates(query_seq, target_seq, min_shared):
    '''(query, target) accession pairs sharing >= min_shared words (every 3rd word of the query).'''
    queries, targets = list(query_seq), list(target_seq)
    t = word_table([target_seq[u] for u in targets], 1)
    frequency = t.groupby('word')['id'].transform('size')
    t = t[frequency <= MAX_TARGETS_PER_WORD]
    q = word_table([query_seq[u] for u in queries], 3)
    pairs = q.merge(t, on='word', suffixes=('_q', '_t')).groupby(['id_q', 'id_t']).size()
    pairs = pairs[pairs >= min_shared].reset_index()
    out = pd.DataFrame({'query': np.array(queries, dtype=object)[pairs['id_q']], 'target': np.array(targets, dtype=object)[pairs['id_t']]})
    return out[out['query'] != out['target']].reset_index(drop=True)


def span(positions, length, min_length):
    '''
    Range of the chain protein aligned for one interface: its interface residues +- SEGMENT_MARGIN, rounded to
    SPAN_STEP (interfaces of similar structures share an alignment), at least min_length long. 1-based inclusive.
    Every interface is mapped with the alignment of its own range: repeats (e.g. the nine ubiquitins of polyubiquitin,
    numbered on different copies in different structures) then align at the copy the structure uses.
    '''
    lo, hi = int(min(positions)) - SEGMENT_MARGIN, int(max(positions)) + SEGMENT_MARGIN
    lo, hi = max(1, (lo - 1) // SPAN_STEP * SPAN_STEP + 1), min(length, -(-hi // SPAN_STEP) * SPAN_STEP)
    missing = min_length - (hi - lo + 1)
    if missing > 0:
        lo, hi = max(1, lo - (missing + 1) // 2), min(length, hi + (missing + 1) // 2)
    return lo, hi


def add_spans(interfaces, query_seq, min_length):
    '''Columns start_1, end_1, start_2, end_2: the aligned range of each side (-1: protein without sequence).'''
    for side in ('1', '2'):
        spans = [span(p, len(query_seq[u]), min_length) if u in query_seq and len(p) else (-1, -1)
                 for u, p in zip(interfaces[f'uniprot_{side}'], interfaces[f'pos_{side}'])]
        interfaces[f'start_{side}'] = [a for a, _ in spans]
        interfaces[f'end_{side}'] = [b for _, b in spans]
    return interfaces


def align_job(job):
    query, target, start, end, q_part, t_seq = job
    q_pos, t_pos, identical = align_positions(q_part, t_seq)
    return query, target, start, end, identical, [p + start - 1 for p in q_pos], t_pos


def homolog_maps(cfg, s, query_seq, target_seq, target_species, spans):
    '''Alignments of the ranges of the chain proteins in structures to the related human/mouse proteins.'''
    pairs = candidates(query_seq, target_seq, s['min_shared_words'])
    log(f'{len(pairs)} candidate pairs (chain protein, human/mouse protein) sharing >= {s["min_shared_words"]} words')
    jobs = [(q, t, a, b, query_seq[q][a - 1:b], target_seq[t])
            for q, t in zip(pairs['query'], pairs['target']) for a, b in sorted(spans.get(q, ()))]
    log(f'{len(jobs)} alignments of the ranges of the chain proteins in structures')
    with Pool(threads()) as pool:
        results = pool.map(align_job, jobs, chunksize=20)
    rows = []
    for query, target, start, end, identical, q_pos, t_pos in results:
        aligned = len(q_pos)
        if aligned >= s['min_aligned'] and 100 * identical / max(aligned, 1) >= s['min_identity']:
            rows.append((query, target, target_species[target], start, end, round(100 * identical / aligned, 1), aligned, q_pos, t_pos))
    maps = pd.DataFrame(rows, columns=['query', 'target', 'species', 'start', 'end', 'identity', 'aligned', 'q_pos', 't_pos'])
    # other UniProt entries of the chain's own gene are not related proteins (their structures are tier pdb)
    genes = {}
    for species in cfg['species']:
        for uniprot, entrez in read_uniprot_entrez(cfg, species).itertuples(index=False):
            genes.setdefault(uniprot, set()).add(entrez)
    gene = lambda u: genes.get(u) or genes.get(u.split('-')[0], set())  # isoform entries (Q674X7-2): their canonical entry
    same_gene = [bool(gene(q) & gene(t)) for q, t in zip(maps['query'], maps['target'])]
    log(f'{sum(same_gene)} range maps to other UniProt entries of the same gene left out')
    maps = maps[~pd.Series(same_gene, index=maps.index, dtype=bool)].reset_index(drop=True)
    pairs = maps[['query', 'target', 'species']].drop_duplicates()
    log(f'{len(maps)} range maps: {pairs["query"].nunique()} chain proteins -> {pairs["target"].nunique()} human/mouse proteins '
        f'(>= {s["min_identity"]}% identity over >= {s["min_aligned"]} residues)', pairs.groupby('species').size().to_dict())
    return maps


def lookup_array(q_pos, t_pos):
    '''Array a with a[query position] = target position (0: not aligned).'''
    a = np.zeros(int(max(q_pos)) + 1, dtype=np.int64)
    a[np.asarray(q_pos, dtype=np.int64)] = t_pos
    return a


def map_positions(positions, lookup):
    if lookup is None:
        return positions
    return np.where(positions < len(lookup), lookup[np.minimum(positions, len(lookup) - 1)], 0)


def transfer(cfg, s, interfaces, maps, target_species):
    '''Interfaces of the chain proteins on the related human/mouse proteins (not on the chains' own proteins).'''
    # chain protein -> target -> (species, {range: (identity, position lookup)})
    options = defaultdict(dict)
    for query, target, species, start, end, identity, q_pos, t_pos in maps[
            ['query', 'target', 'species', 'start', 'end', 'identity', 'q_pos', 't_pos']].itertuples(index=False):
        options[query].setdefault(target, (species, {}))[1][(start, end)] = (identity, lookup_array(q_pos, t_pos))

    def side(protein, rng):
        '''(target, species, identity, lookup) for one side of an interface: related proteins and the protein itself.'''
        out = [(target, species, *ranges[rng]) for target, (species, ranges) in options.get(protein, {}).items() if rng in ranges]
        if protein in target_species:
            out.append((protein, target_species[protein], 100.0, None))
        return out

    min_contacts = cfg['pdb']['min_interface_contacts']
    merged = defaultdict(set)
    entries = defaultdict(set)
    identity = {}
    for (u1, u2, a1, b1, a2, b2), group in interfaces.groupby(['uniprot_1', 'uniprot_2', 'start_1', 'end_1', 'start_2', 'end_2']):
        side_1, side_2 = side(u1, (a1, b1)), side(u2, (a2, b2))
        if not side_1 or not side_2:
            continue
        p1 = np.concatenate(group['pos_1'].to_numpy()).astype(np.int64)
        p2 = np.concatenate(group['pos_2'].to_numpy()).astype(np.int64)
        pairs = np.unique(np.stack([p1, p2], axis=1), axis=0)
        for t1, species_1, identity_1, map_1 in side_1:
            for t2, species_2, identity_2, map_2 in side_2:
                if species_1 != species_2 or (map_1 is None and map_2 is None):
                    continue
                a, b = map_positions(pairs[:, 0], map_1), map_positions(pairs[:, 1], map_2)
                ok = (a > 0) & (b > 0)
                if ok.sum() < min_contacts or ok.sum() < s['min_mapped'] * len(pairs):
                    continue
                a, b = a[ok], b[ok]
                if t1 > t2:
                    t1_, t2_, a, b = t2, t1, b, a
                else:
                    t1_, t2_ = t1, t2
                key = (species_1, t1_, t2_)
                merged[key].update(zip(a.tolist(), b.tolist()))
                entries[key].update(group['pdb'])
                identity[key] = min(identity.get(key, 100.0), identity_1, identity_2)
    rows = []
    for key, residue_pairs in merged.items():
        residue_pairs = sorted(residue_pairs)
        found = sorted(entries[key])
        rows.append((*key, [a for a, _ in residue_pairs], [b for _, b in residue_pairs], len(residue_pairs), len(found),
                     ' '.join(found[:20]), identity[key]))
    return pd.DataFrame(rows, columns=['species', 'uniprot_1', 'uniprot_2', 'pos_1', 'pos_2', 'n_pairs', 'n_entries', 'pdb', 'min_identity'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    args = parser.parse_args()
    cfg = load_config(args.config)
    s = settings(cfg)
    if not s['enabled']:
        raise SystemExit('Interfaces of related proteins are not enabled (homologs.enabled in the configuration)')
    folder = paths(cfg)['intermediate'] / 'pdb'

    # interfaces of all species (step 5), oriented uniprot_1 <= uniprot_2
    files = [folder / 'contacts' / f'shard_{shard:05d}.interfaces.parquet' for shard in range(cfg['pdb']['shards'])]
    missing = [f.name for f in files if not f.exists()]
    if missing:
        raise SystemExit(f'{len(missing)} shards of step 5 are missing (rerun these array tasks): {missing[:20]}')
    parts = []
    for path in files:
        part = pd.read_parquet(path, columns=['pdb', 'uniprot_1', 'uniprot_2', 'pos_1', 'pos_2', 'n_pairs'])
        part = part[part['n_pairs'] >= cfg['pdb']['min_interface_contacts']]
        swap = part['uniprot_1'] > part['uniprot_2']
        part.loc[swap, ['uniprot_1', 'uniprot_2', 'pos_1', 'pos_2']] = part.loc[swap, ['uniprot_2', 'uniprot_1', 'pos_2', 'pos_1']].values
        parts.append(part.drop(columns='n_pairs'))
    interfaces = pd.concat(parts, ignore_index=True)
    log(f'{len(interfaces)} chain interfaces of {interfaces["pdb"].nunique()} entries')

    # human / mouse proteins that can be placed on exons, and the chain proteins
    target_seq, target_species = {}, {}
    for species in cfg['species']:
        proteins = placeable_proteins(cfg, species).drop_duplicates('uniprot')
        target_seq.update(zip(proteins['uniprot'], proteins['sequence']))
        target_species.update(dict.fromkeys(proteins['uniprot'], species))
    chain_proteins = set(interfaces['uniprot_1']) | set(interfaces['uniprot_2'])
    query_seq, _ = uniprot_sequences(cfg, chain_proteins)
    for protein in chain_proteins - set(query_seq):
        if protein in target_seq:
            query_seq[protein] = target_seq[protein]
    log(f'{len(chain_proteins)} chain proteins, {len(query_seq)} with a UniProt sequence; {len(target_seq)} human/mouse proteins')

    interfaces = add_spans(interfaces, query_seq, s['min_aligned'] + 10)
    spans = defaultdict(set)
    for side in ('1', '2'):
        for protein, a, b in zip(interfaces[f'uniprot_{side}'], interfaces[f'start_{side}'], interfaces[f'end_{side}']):
            if a > 0:
                spans[protein].add((a, b))
    maps = homolog_maps(cfg, s, query_seq, target_seq, target_species, spans)
    maps.to_parquet(folder / 'homolog_maps.parquet')
    interfaces = interfaces[interfaces['uniprot_1'].isin(set(maps['query'])) | interfaces['uniprot_2'].isin(set(maps['query']))]
    result = transfer(cfg, s, interfaces, maps, target_species)
    result.to_parquet(folder / 'interfaces_homolog.parquet')
    log('interfaces of related proteins:', result.groupby('species').size().to_dict(), 'protein pairs')
    log(f'written: {folder / "interfaces_homolog.parquet"}')
