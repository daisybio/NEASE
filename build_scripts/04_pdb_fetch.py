"""Step 4: download the selected PDB entries (mmCIF, asymmetric unit) from RCSB.

Not needed when config pdb.mirror points to a local copy of the wwPDB archive (then this step only reports
missing entries). Needs internet access. Files that are already there are skipped, so the step can be rerun.

usage: python 04_pdb_fetch.py --config config.yaml [--threads 8] [--task 0 --ntasks 10]
"""
import argparse
import os
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from common import download, load_config, log, paths, pdb_file


def fetch(cfg, pdb):
    path = pdb_file(cfg, pdb)
    if path.exists():
        return 'present'
    if cfg['pdb'].get('mirror'):
        return 'missing'
    try:
        download(cfg['sources']['pdb_entry'].format(pdb=pdb), path, retries=3, timeout=300)
        return 'downloaded'
    except RuntimeError:
        return 'failed'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--task', type=int, help='part of the entries to download (with --ntasks; default: $SLURM_ARRAY_TASK_ID or 0)')
    parser.add_argument('--ntasks', type=int, default=1)
    args = parser.parse_args()
    cfg = load_config(args.config)
    task = args.task if args.task is not None else int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))

    entries = pd.read_csv(paths(cfg)['intermediate'] / 'pdb' / 'entries.tsv', sep='\t', dtype=str)['pdb'].tolist()
    entries = entries[task::args.ntasks]
    log(f'{len(entries)} entries')
    with ThreadPoolExecutor(args.threads) as pool:
        status = pd.Series(list(pool.map(lambda pdb: fetch(cfg, pdb), entries)), index=entries)
    log(status.value_counts().to_dict())
    if (status == 'failed').any():
        log('failed:', ' '.join(status.index[status == 'failed'][:50]))
