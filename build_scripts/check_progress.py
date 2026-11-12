"""Progress of the downloads of steps 4 (PDB entries) and 8 (AlphaFold models): files present / selected.

usage: python check_progress.py config.yaml
"""
import sys

import pandas as pd

from common import load_config, paths

cfg = load_config(sys.argv[1])
p = paths(cfg)
entries = pd.read_csv(p['intermediate'] / 'pdb' / 'entries.tsv', sep='\t', dtype=str)['pdb']
models = pd.read_parquet(p['intermediate'] / 'alphafold' / 'models.parquet', columns=['model', 'tar'])
pdb_done = sum((p['pdb'] / e[1:3] / f'{e}.cif.gz').exists() for e in entries)
models['done'] = [(p['alphafold'] / 'models' / f'{m}.cif.zst').exists() for m in models['model']]
tars = models.groupby('tar')['done'].all()
print(f'PDB entries:       {pdb_done} / {len(entries)}  ({len(entries) - pdb_done} left)')
print(f'AlphaFold models:  {models["done"].sum()} / {len(models)}  ({(~models["done"]).sum()} left)')
print(f'AlphaFold tars:    {tars.sum()} / {len(tars)}  ({(~tars).sum()} left)')
