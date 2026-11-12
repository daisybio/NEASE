"""Step 1: download the input tables (needs internet access: run it on a login/transfer node if compute nodes have none).

PDB entries and AlphaFold models are fetched by 04_pdb_fetch.py and 08_af_fetch.py.
Files that are already there are kept (--force downloads them again). A manifest with the download
dates is written to <work_dir>/downloads/manifest.tsv.

usage: python 01_download.py --config config.yaml [--only biogrid intact ...] [--force]
"""
import argparse
from datetime import datetime

import pandas as pd

from common import HUMAN_ONLY, download, input_path, load_config, log, paths, per_species, source_url

# Sources that are unreliable and not needed (HuRI is also part of IntAct)
OPTIONAL = {'huri'}


def needed_sources(cfg):
    names = ['sifts_pfam', 'sifts_uniprot', 'pfam_dead', 'ppidm', 'elm_instances', 'elm_interactions', 'elm_domains',
             'uniprot_idmapping', 'uniprot_proteome', 'uniprot_sprot', 'uniprot_varsplic', 'ensembl_gtf', 'ensembl_pep', 'ensembl_uniprot', 'ensembl_entrez']
    for source in cfg['ppi']['sources']:
        names += {'bioplex': ['bioplex_293t', 'bioplex_hct116'], 'string': ['string_physical', 'string_aliases']}.get(source, [source])
    if (cfg.get('boltz') or {}).get('enabled'):
        names.append('boltz')
    if 'Mouse' in cfg['species']:
        names.append('compara_homologies')
    if cfg['alphafold']['enabled']:
        names.append('afdb_heterodimers')
        if cfg['alphafold']['homodimers']:
            names.append('afdb_homodimers')
    return names


def jobs(cfg, names):
    for name in names:
        if not per_species(name):
            if name in HUMAN_ONLY and 'Human' not in cfg['species']:
                continue
            yield name, None
        else:
            for species in cfg['species']:
                yield name, species


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    parser.add_argument('--only', nargs='+', help='download only these sources')
    parser.add_argument('--force', action='store_true', help='download again even if the file exists')
    args = parser.parse_args()
    cfg = load_config(args.config)

    names = args.only or needed_sources(cfg)
    manifest_path = paths(cfg)['downloads'] / 'manifest.tsv'
    manifest = pd.read_csv(manifest_path, sep='\t') if manifest_path.exists() else pd.DataFrame(
        columns=['source', 'species', 'url', 'file', 'bytes', 'downloaded'])
    failed = []
    for name, species in jobs(cfg, names):
        url, dest = source_url(cfg, name, species), input_path(cfg, name, species)
        fresh = args.force or not dest.exists()
        log(f'{name} {species or ""}: {url}')
        try:
            download(url, dest, retries=3 if name in OPTIONAL else 5, force=args.force)
        except RuntimeError as e:
            log(f'  FAILED: {e}' + (' (optional source, skipped)' if name in OPTIONAL else ''))
            if name not in OPTIONAL:
                failed.append(f'{name} {species or ""}')
            continue
        if fresh:
            manifest = manifest[~((manifest['source'] == name) & (manifest['species'].fillna('') == (species or '')))]
            manifest.loc[len(manifest)] = [name, species, url, str(dest), dest.stat().st_size, datetime.now().isoformat(timespec='seconds')]
            manifest.to_csv(manifest_path, sep='\t', index=False)
    if failed:
        raise SystemExit(f'Downloads failed: {failed}')
    log('done')
