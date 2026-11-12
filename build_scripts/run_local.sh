#!/bin/bash
# Run the NEASE network build on one machine (no Slurm), array steps in parallel processes.
#
#   bash run_local.sh config.yaml [steps]        e.g. bash run_local.sh config.yaml 11 12 13 14
#
# Activate the build environment first (conda activate nease-network-build). JOBS: parallel processes.
set -euo pipefail

CONFIG=$(realpath "${1:?usage: run_local.sh config.yaml [steps]}")
shift
STEPS=${*:-"1 2 3 4 5 6 7 8 9 10 11 12 15 17 18 19 13 14 16"}
JOBS=${JOBS:-8}
cd "$(dirname "$0")"

read_config() { python -c "import yaml,sys; c=yaml.safe_load(open('$CONFIG')); print(eval(sys.argv[1]))" "$1"; }
PDB_SHARDS=$(read_config "c['pdb']['shards']")
AF_TASKS=$(read_config "c['alphafold']['tasks']")
AF_ENABLED=$(read_config "c['alphafold']['enabled']")
I3D_ENABLED=$(read_config "(c.get('interactome3d') or {}).get('enabled', False)")
HOMOLOGS_ENABLED=$(read_config "(c.get('homologs') or {}).get('enabled', False)")
HUMANPPI_ENABLED=$(read_config "(c.get('humanppi') or {}).get('enabled', False)")
BOLTZ_ENABLED=$(read_config "(c.get('boltz') or {}).get('enabled', False)")

array() {  # array <script> <n tasks>
    seq 0 $(($2 - 1)) | xargs -P "$JOBS" -I{} python "$1" --config "$CONFIG" --task {}
}

for step in $STEPS; do
    echo "=== step $step"
    case $step in
        1) python 01_download.py --config "$CONFIG" ;;
        2) python 02_ensembl.py --config "$CONFIG" ;;
        3) python 03_pdb_select.py --config "$CONFIG" ;;
        4) python 04_pdb_fetch.py --config "$CONFIG" --threads "$JOBS" ;;
        5) array 05_pdb_contacts.py "$PDB_SHARDS" ;;
        6) python 06_pdb_merge.py --config "$CONFIG" ;;
        7) [[ "$AF_ENABLED" == True ]] && python 07_af_select.py --config "$CONFIG" ;;
        8) [[ "$AF_ENABLED" == True ]] && array 08_af_fetch.py "$AF_TASKS" ;;
        9) [[ "$AF_ENABLED" == True ]] && array 09_af_contacts.py "$AF_TASKS" ;;
        10) python 10_orthologs.py --config "$CONFIG" ;;
        11) python 11_ppi.py --config "$CONFIG" ;;
        12) python 12_elm.py --config "$CONFIG" ;;
        15) [[ "$I3D_ENABLED" == True ]] && python 15_interactome3d.py --config "$CONFIG" --threads "$JOBS" ;;
        17) [[ "$HOMOLOGS_ENABLED" == True ]] && python 17_pdb_homologs.py --config "$CONFIG" ;;
        18) [[ "$HUMANPPI_ENABLED" == True ]] && python 18_humanppi.py --config "$CONFIG" ;;
        19) [[ "$BOLTZ_ENABLED" == True ]] && python 19_boltz.py --config "$CONFIG" ;;
        13) python 13_assemble.py --config "$CONFIG" ;;
        14) python 14_validate.py --config "$CONFIG" ;;
        16) python 16_diagnose_losses.py --config "$CONFIG" ;;
    esac
done
