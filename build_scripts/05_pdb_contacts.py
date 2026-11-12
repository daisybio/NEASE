"""Step 5 (array job): residue contacts between protein chains of PDB entries.

For every entry of the shard, the first biological assembly is built (contacts of the crystal lattice are not
used) and residue pairs of different chains with heavy atoms within pdb.contact_distance are collected.
They give
    domain contacts     Pfam domain pairs in contact (SIFTS domain boundaries): the domain-domain catalogue
    interfaces          residue pairs in UniProt numbering, for chains mapped to UniProt: residue-level
                        interfaces and domain contacts of the human/mouse proteins

Writes <work_dir>/intermediate/pdb/contacts/shard_<n>.{domains,interfaces,log}.parquet

usage: python 05_pdb_contacts.py --config config.yaml [--task <shard>]   (Slurm: the array task ID)
"""
import argparse
import time
from collections import Counter, defaultdict

import gemmi
import pandas as pd

from common import download, load_config, log, paths, pdb_file, task_id
from contacts import build_assembly, residue_contacts


def position_index(rows, start, end, value):
    '''{chain: {label_seq: [value, ...]}} for residue ranges.'''
    index = defaultdict(lambda: defaultdict(list))
    for row in rows.itertuples(index=False):
        for s in range(getattr(row, start), getattr(row, end) + 1):
            index[row.chain][s].append(value(row, s))
    return index


def analyse(cfg, pdb, pfam, uniprot):
    settings = cfg['pdb']
    path = pdb_file(cfg, pdb)
    if not path.exists():
        if settings.get('mirror'):
            return None, None, dict(status='missing')
        download(cfg['sources']['pdb_entry'].format(pdb=pdb), path, retries=3, timeout=300)
    st = gemmi.read_structure(str(path))
    method = st.meta.experiments[0].method if st.meta.experiments else ''
    info = dict(method=method, resolution=st.resolution)
    if method not in settings['methods']:
        return None, None, dict(info, status='method')
    if method != 'SOLUTION NMR' and not 0 < st.resolution <= settings['max_resolution']:
        return None, None, dict(info, status='resolution')
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    model, kind = build_assembly(st, settings['max_assembly_chains'])
    info.update(assembly=kind, n_chains=len(model))
    if len(model) < 2:
        return None, None, dict(info, status='single chain')
    if model.count_atom_sites() > settings['max_atoms']:
        return None, None, dict(info, status='too large')

    # merge the copies of a chain pair (assembly copies and symmetric pairs)
    pairs = defaultdict(set)
    for (c1, c2), residues in residue_contacts(model, settings['contact_distance']).items():
        o1, o2 = c1.split('~')[0], c2.split('~')[0]
        if o1 <= o2:
            pairs[(o1, o2)] |= residues
        else:
            pairs[(o2, o1)] |= {(b, a) for a, b in residues}
    info['n_pairs'] = sum(len(p) for p in pairs.values())

    domains = position_index(pfam, 'label_start', 'label_end', lambda row, s: (row.pfam, row.uniprot, row.label_start))
    positions = position_index(uniprot, 'label_beg', 'label_end', lambda row, s: (row.uniprot, row.sp_beg + s - row.label_beg))
    domain_rows, interface_rows = [], []
    for (o1, o2), residues in pairs.items():
        counts = Counter((d1, d2) for s1, s2 in residues for d1 in domains[o1].get(s1, []) for d2 in domains[o2].get(s2, []))
        for ((pfam_1, uniprot_1, _), (pfam_2, uniprot_2, _)), n in counts.items():
            domain_rows.append((pdb, o1, o2, pfam_1, uniprot_1, pfam_2, uniprot_2, n))
        mapped = defaultdict(list)
        for s1, s2 in residues:
            for u1, p1 in positions[o1].get(s1, []):
                for u2, p2 in positions[o2].get(s2, []):
                    mapped[(u1, u2)].append((p1, p2))
        for (u1, u2), residue_pairs in mapped.items():
            residue_pairs = sorted(set(residue_pairs))
            interface_rows.append((pdb, o1, o2, u1, u2, [p for p, _ in residue_pairs], [p for _, p in residue_pairs], len(residue_pairs)))
    return domain_rows, interface_rows, dict(info, status='ok')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    parser.add_argument('--task', type=int, help='shard (default: $SLURM_ARRAY_TASK_ID)')
    args = parser.parse_args()
    cfg = load_config(args.config)
    shard = task_id(args)
    folder = paths(cfg)['intermediate'] / 'pdb'
    out = folder / 'contacts'
    out.mkdir(exist_ok=True)

    entries = pd.read_csv(folder / 'entries.tsv', sep='\t', dtype={'pdb': str})
    entries = entries.loc[entries['shard'] == shard, 'pdb'].tolist()
    pfam = pd.read_parquet(folder / 'pfam_map.parquet', filters=[('pdb', 'in', entries)]) if entries else None
    uniprot = pd.read_parquet(folder / 'uniprot_map.parquet', filters=[('pdb', 'in', entries)]) if entries else None
    log(f'shard {shard}: {len(entries)} entries')

    domain_rows, interface_rows, logs = [], [], []
    pfam_by_entry = dict(tuple(pfam.groupby('pdb'))) if entries else {}
    uniprot_by_entry = dict(tuple(uniprot.groupby('pdb'))) if entries else {}
    empty = pd.DataFrame(columns=['chain', 'label_start', 'label_end', 'pfam', 'uniprot', 'label_beg', 'sp_beg'])
    for i, pdb in enumerate(entries):
        start = time.time()
        try:
            d, f, info = analyse(cfg, pdb, pfam_by_entry.get(pdb, empty), uniprot_by_entry.get(pdb, empty))
            domain_rows += d or []
            interface_rows += f or []
        except Exception as e:  # damaged or unusual files: log and continue
            info = dict(status=f'error: {type(e).__name__}: {e}'[:200])
        logs.append(dict(pdb=pdb, seconds=round(time.time() - start, 2), **info))
        if (i + 1) % 100 == 0:
            log(f'  {i + 1}/{len(entries)}')

    pd.DataFrame(domain_rows, columns=['pdb', 'chain_1', 'chain_2', 'pfam_1', 'uniprot_1', 'pfam_2', 'uniprot_2', 'n_contacts']
                 ).to_parquet(out / f'shard_{shard:05d}.domains.parquet')
    pd.DataFrame(interface_rows, columns=['pdb', 'chain_1', 'chain_2', 'uniprot_1', 'uniprot_2', 'pos_1', 'pos_2', 'n_pairs']
                 ).to_parquet(out / f'shard_{shard:05d}.interfaces.parquet')
    logs = pd.DataFrame(logs)
    logs.to_parquet(out / f'shard_{shard:05d}.log.parquet')
    log(f'shard {shard}: {logs["status"].value_counts().to_dict()}')
