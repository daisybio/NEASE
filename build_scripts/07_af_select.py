"""Step 7: select confident AlphaFold DB complex models of human / mouse protein pairs.

The AlphaFold DB complexes (EMBL-EBI / NVIDIA collaboration) come with interface confidence scores. Models of
two proteins of the same species with ipSAE >= alphafold.min_ipsae and pDockQ2 >= alphafold.min_pdockq2
(AFDB's high-confidence definition) are kept; homodimers only with alphafold.homodimers: true.

Writes <work_dir>/intermediate/alphafold/models.parquet: model, kind, species, uniprot_1, uniprot_2, scores, tar

usage: python 07_af_select.py --config config.yaml [--max-tars 2]
"""
import argparse

import pandas as pd

from common import SPECIES, input_path, load_config, log, paths

HETERODIMER_COLUMNS = ['modelEntityId', 'uniprot_ac_1', 'uniprot_ac_2', 'tax_id_1', 'tax_id_2', 'ipTM', 'pDockQ',
                       'max_ipSAE', 'max_pDockQ2_AB', 'local_tar_name']
HOMODIMER_COLUMNS = ['modelEntityId', 'uniprotAccession', 'taxId', 'ipTM', 'pDockQ', 'ipSAE_AB', 'ipSAE_BA',
                     'pDockQ2_AB', 'pDockQ2_BA', 'local_tar_name']


def heterodimers(path, taxa, settings):
    parts = []
    for chunk in pd.read_csv(path, usecols=HETERODIMER_COLUMNS, chunksize=1_000_000, dtype={'local_tar_name': str}):
        chunk = chunk[chunk['tax_id_1'].isin(taxa) & (chunk['tax_id_1'] == chunk['tax_id_2'])
                      & (chunk['max_ipSAE'] >= settings['min_ipsae']) & (chunk['max_pDockQ2_AB'] >= settings['min_pdockq2'])]
        parts.append(chunk.rename(columns={'modelEntityId': 'model', 'uniprot_ac_1': 'uniprot_1', 'uniprot_ac_2': 'uniprot_2',
                                           'tax_id_1': 'taxid', 'max_ipSAE': 'ipSAE', 'max_pDockQ2_AB': 'pDockQ2',
                                           'local_tar_name': 'tar'}).drop(columns='tax_id_2'))
    return pd.concat(parts, ignore_index=True).assign(kind='heterodimers')


def homodimers(path, taxa, settings):
    parts = []
    for chunk in pd.read_csv(path, usecols=HOMODIMER_COLUMNS, chunksize=1_000_000, dtype={'local_tar_name': str}):
        chunk = chunk.assign(ipSAE=chunk[['ipSAE_AB', 'ipSAE_BA']].max(axis=1), pDockQ2=chunk[['pDockQ2_AB', 'pDockQ2_BA']].max(axis=1))
        chunk = chunk[chunk['taxId'].isin(taxa) & (chunk['ipSAE'] >= settings['min_ipsae']) & (chunk['pDockQ2'] >= settings['min_pdockq2'])]
        parts.append(pd.DataFrame({'model': chunk['modelEntityId'], 'uniprot_1': chunk['uniprotAccession'],
                                   'uniprot_2': chunk['uniprotAccession'], 'taxid': chunk['taxId'], 'ipTM': chunk['ipTM'],
                                   'pDockQ': chunk['pDockQ'], 'ipSAE': chunk['ipSAE'], 'pDockQ2': chunk['pDockQ2'],
                                   'tar': chunk['local_tar_name']}))
    return pd.concat(parts, ignore_index=True).assign(kind='homodimers')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    parser.add_argument('--max-tars', type=int, help='only the models of a few tar files (for testing)')
    args = parser.parse_args()
    cfg = load_config(args.config)
    settings = cfg['alphafold']
    out = paths(cfg)['intermediate'] / 'alphafold'
    out.mkdir(exist_ok=True)

    taxa = {SPECIES[s]['taxid']: s for s in cfg['species']}
    models = [heterodimers(input_path(cfg, 'afdb_heterodimers'), taxa, settings)]
    if settings['homodimers']:
        models.append(homodimers(input_path(cfg, 'afdb_homodimers'), taxa, settings))
    models = pd.concat(models, ignore_index=True)
    models['species'] = models['taxid'].map(taxa)
    models['tar'] = models['tar'].where(models['tar'].str.endswith('.tar'), models['tar'] + '.tar')
    if args.max_tars:
        keep = models.groupby(['kind', 'tar']).size().sort_values(ascending=False).index[:args.max_tars]
        models = models.set_index(['kind', 'tar']).loc[keep].reset_index()
    models.to_parquet(out / 'models.parquet')
    log(f'{len(models)} models:', models.groupby(['species', 'kind']).size().to_dict(),
        f'in {models[["kind", "tar"]].drop_duplicates().shape[0]} tar files')
