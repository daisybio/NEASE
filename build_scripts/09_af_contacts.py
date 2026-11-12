"""Step 9 (array job): residue-level interfaces of the AlphaFold complex models.

Residue pairs of the two chains with heavy atoms within alphafold.contact_distance, both residues with
pLDDT >= alphafold.min_plddt. Residue numbers of AFDB models are positions in the UniProt sequence.

Writes <work_dir>/intermediate/alphafold/contacts/task_<n>.parquet:
    model, species, uniprot_1, uniprot_2, pos_1, pos_2 (residue pairs), n_pairs

usage: python 09_af_contacts.py --config config.yaml [--task 0] [--ntasks 200]   (Slurm: the array task ID)
"""
import argparse

import gemmi
import pandas as pd
import zstandard

from common import load_config, log, paths, task_id
from contacts import drop_low_confidence, residue_contacts


def model_interface(path, settings):
    with open(path, 'rb') as f:
        text = zstandard.ZstdDecompressor().stream_reader(f).read().decode()
    block = gemmi.cif.read_string(text).sole_block()
    entity_uniprot = {row[0]: (row[1], int(row[2])) for row in
                      block.find('_ma_target_ref_db_details.', ['target_entity_id', 'db_accession', 'seq_db_align_begin'])}
    chain_entity = {row[0]: row[1] for row in block.find('_struct_asym.', ['id', 'entity_id'])}
    st = gemmi.make_structure_from_block(block)
    st.setup_entities()
    st.remove_hydrogens()
    model = st[0]
    # uniprot and offset of every chain (via its label_asym_id)
    chains = {chain.name: entity_uniprot[chain_entity[chain[0].subchain]] for chain in model}
    drop_low_confidence(model, settings['min_plddt'])
    rows = []
    for (c1, c2), residues in residue_contacts(model, settings['contact_distance']).items():
        (u1, o1), (u2, o2) = chains[c1], chains[c2]
        pairs = sorted((s1 + o1 - 1, s2 + o2 - 1) for s1, s2 in residues)
        if u1 > u2:
            u1, u2, pairs = u2, u1, sorted((b, a) for a, b in pairs)
        rows.append((u1, u2, [a for a, _ in pairs], [b for _, b in pairs], len(pairs)))
    return rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    parser.add_argument('--task', type=int, help='default: $SLURM_ARRAY_TASK_ID')
    parser.add_argument('--ntasks', type=int, help='default: alphafold.tasks of the configuration')
    args = parser.parse_args()
    cfg = load_config(args.config)
    settings = cfg['alphafold']
    task, ntasks = task_id(args), args.ntasks or settings['tasks']
    folder = paths(cfg)['intermediate'] / 'alphafold'
    (folder / 'contacts').mkdir(exist_ok=True)

    models = pd.read_parquet(folder / 'models.parquet').sort_values('model').iloc[task::ntasks]
    log(f'task {task}/{ntasks}: {len(models)} models')
    rows, missing, failed = [], 0, 0
    for model in models.itertuples(index=False):
        path = paths(cfg)['alphafold'] / 'models' / f'{model.model}.cif.zst'
        if not path.exists():
            missing += 1
            continue
        try:
            for u1, u2, pos_1, pos_2, n in model_interface(path, settings):
                if n >= settings['min_interface_contacts']:
                    rows.append((model.model, model.species, u1, u2, pos_1, pos_2, n, model.ipSAE, model.pDockQ2))
        except Exception as e:  # damaged files: report and continue
            failed += 1
            log(f'  {model.model}: {type(e).__name__}: {e}')
    table = pd.DataFrame(rows, columns=['model', 'species', 'uniprot_1', 'uniprot_2', 'pos_1', 'pos_2', 'n_pairs', 'ipSAE', 'pDockQ2'])
    table.to_parquet(folder / 'contacts' / f'task_{task:05d}.parquet')
    log(f'task {task}: {len(table)} interfaces, {missing} models missing (rerun 08_af_fetch.py), {failed} failed')
