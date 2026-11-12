"""Step 8 (array job, needs internet access): download the selected AlphaFold models.

Two ways (alphafold.download in the configuration):
  models  single model files from the AlphaFold DB website (default; ~1 MB per model, ~10 GB for human and mouse)
  tars    the tar files of the FTP collection that contain the models. They hold thousands of models each and are
          large (median 1.3 GB; ~3 TB for the human and mouse models), so this takes days. Each tar file is deleted
          after the selected models are extracted.
Models are kept as <work_dir>/alphafold/models/<model>.cif.zst. Models that are already there are skipped, and a
model or tar file that fails is skipped and reported at the end, so failed tasks can simply be rerun.

usage: python 08_af_fetch.py --config config.yaml [--task 0] [--ntasks 200] [--threads 2]   (Slurm: the array task ID)
"""
import argparse
import tarfile
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import zstandard

from common import download, load_config, log, paths, task_id

MODEL_URL = 'https://alphafold.ebi.ac.uk/files/{model}-model_v1.cif'  # if sources.afdb_model is not configured


def fetch_model(cfg, folder, model):
    '''Download one model file and store it zstd-compressed, like the models in the tar files.'''
    raw = folder / 'models' / f'{model}.cif'
    target = folder / 'models' / f'{model}.cif.zst'
    try:
        download(cfg['sources'].get('afdb_model', MODEL_URL).format(model=model), raw, retries=3)
        target.with_suffix('.part').write_bytes(zstandard.ZstdCompressor().compress(raw.read_bytes()))
        target.with_suffix('.part').replace(target)
        return 'downloaded'
    except (RuntimeError, OSError) as e:
        log(f'  {model}: FAILED ({type(e).__name__}: {e})')
        return 'failed'
    finally:
        raw.unlink(missing_ok=True)


def fetch_models(cfg, folder, models, threads):
    missing = sorted(m for m in models['model'] if not (folder / 'models' / f'{m}.cif.zst').exists())
    log(f'{len(models)} models, {len(missing)} to download')
    if not missing:
        return []
    # first run: one model first; if the website refuses the downloads, stop instead of failing every model
    if len(missing) == len(models):
        if fetch_model(cfg, folder, missing[0]) == 'failed':
            raise SystemExit(f'Model files cannot be downloaded from {cfg["sources"].get("afdb_model", MODEL_URL)} (see above). '
                             'Set alphafold.download: tars in the configuration to use the tar files instead.')
        missing = missing[1:]
    with ThreadPoolExecutor(threads) as pool:
        status = pd.Series(list(pool.map(lambda m: fetch_model(cfg, folder, m), missing)), index=missing, dtype=object)
    log(status.value_counts().to_dict())
    return status.index[status == 'failed'].tolist()


def fetch_tars(cfg, folder, models):
    tars = sorted(models[['kind', 'tar']].drop_duplicates().itertuples(index=False, name=None))
    log(f'{len(tars)} tar files')
    failed = []
    for kind, tar in tars:
        wanted = set(models.loc[(models['kind'] == kind) & (models['tar'] == tar), 'model'])
        missing = {m for m in wanted if not (folder / 'models' / f'{m}.cif.zst').exists()}
        if not missing:
            continue
        path = folder / 'tars' / f'{kind}_{tar}'
        found = 0
        try:
            download(cfg['sources']['afdb_tar'].format(kind=kind, tar=tar), path)
            with tarfile.open(path) as archive:
                for member in archive:
                    # every model is there as mmCIF (.cif.zst) and PDB format (.pdb.zst)
                    name = member.name.split('/')[-1]
                    model = name.split('-model')[0]
                    if model in missing and name.endswith('.cif.zst'):
                        data = archive.extractfile(member).read()
                        target = folder / 'models' / f'{model}.cif.zst'
                        target.with_suffix('.part').write_bytes(data)
                        target.with_suffix('.part').replace(target)
                        found += 1
        except (RuntimeError, OSError, tarfile.TarError) as e:
            # skip this tar file and go on with the others; a rerun fetches it again
            log(f'  {kind}/{tar}: FAILED ({type(e).__name__}: {e})')
            failed.append(f'{kind}/{tar}')
            continue
        finally:
            path.unlink(missing_ok=True)
        log(f'  {kind}/{tar}: {found} of {len(missing)} models extracted')
    return failed


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    parser.add_argument('--task', type=int, help='default: $SLURM_ARRAY_TASK_ID')
    parser.add_argument('--ntasks', type=int, help='default: alphafold.tasks of the configuration')
    parser.add_argument('--threads', type=int, default=2, help='parallel downloads of model files (download: models)')
    args = parser.parse_args()
    cfg = load_config(args.config)
    task, ntasks = task_id(args), args.ntasks or cfg['alphafold']['tasks']
    mode = cfg['alphafold'].get('download', 'models')
    folder = paths(cfg)['alphafold']
    (folder / 'models').mkdir(exist_ok=True)
    (folder / 'tars').mkdir(exist_ok=True)

    models = pd.read_parquet(paths(cfg)['intermediate'] / 'alphafold' / 'models.parquet', columns=['model', 'kind', 'tar'])
    log(f'task {task}/{ntasks}, download: {mode}')
    if mode == 'models':
        failed = fetch_models(cfg, folder, models.sort_values('model').iloc[task::ntasks], args.threads)
    else:
        tars = sorted(models[['kind', 'tar']].drop_duplicates().itertuples(index=False, name=None))[task::ntasks]
        failed = fetch_tars(cfg, folder, models.merge(pd.DataFrame(tars, columns=['kind', 'tar'])))
    if failed:
        raise SystemExit(f'task {task}: {len(failed)} {mode} failed (rerun this task): {failed[:20]}')
    log(f'task {task}: done')
