#!/bin/bash
#SBATCH --job-name=nease_submit
#SBATCH --partition=day
#SBATCH --time=00:15:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=nease_submit_%j.out
# Submit the whole NEASE network build to Slurm, with the dependencies between the steps.
# From the build_scripts folder, either as a job or directly on a login node:
#
#   sbatch slurm/submit_all.sh config.yaml            all steps
#   sbatch slurm/submit_all.sh config.yaml 5 6 13     only some steps (no dependencies on steps that are not submitted)
#   bash slurm/submit_all.sh config.yaml              the same without a job
#
# This script only submits the jobs of the steps and ends within a minute; the steps then run as their own jobs.
# Its own log: nease_submit_<job ID>.out in the folder it was submitted from (sbatch).
# Logs of the steps: <work_dir>/logs/<job name>_<job ID>[_<array task>].out
#
# Adapt the settings below to your cluster. Steps 1, 4 and 8 download data: if compute nodes have no
# internet access, set DOWNLOAD_PARTITION to a partition that has (e.g. a transfer partition), or run these
# steps on a login/transfer node and submit the rest with the step numbers.
set -euo pipefail

CONFIG=$(realpath "${1:?usage: submit_all.sh config.yaml [steps]}")
shift
STEPS=${*:-"1 2 3 4 5 6 7 8 9 10 11 12 15 17 18 19 13 14 16"}

PARTITION=${PARTITION:-day}                       # compute partition
DOWNLOAD_PARTITION=${DOWNLOAD_PARTITION:-$PARTITION}   # steps 1, 4, 8: needs internet access
export NEASE_BUILD_ACTIVATE=${NEASE_BUILD_ACTIVATE:-"module load miniconda; conda activate nease-network-build"}
if [[ -z "${NEASE_BUILD_SCRIPTS:-}" && -n "${SLURM_JOB_ID:-}" ]]; then
    # under sbatch this script runs from a copy in Slurm's spool folder: find the scripts from the submission folder
    for dir in "$SLURM_SUBMIT_DIR" "$SLURM_SUBMIT_DIR/build_scripts" "$SLURM_SUBMIT_DIR/.."; do
        [[ -f "$dir/common.py" ]] && NEASE_BUILD_SCRIPTS=$(cd "$dir" && pwd) && break
    done
    : "${NEASE_BUILD_SCRIPTS:?build scripts not found: submit from the build_scripts folder (sbatch slurm/submit_all.sh config.yaml)}"
fi
export NEASE_BUILD_SCRIPTS=${NEASE_BUILD_SCRIPTS:-$(cd "$(dirname "$0")/.." && pwd)}

for p in "$PARTITION" "$DOWNLOAD_PARTITION"; do
    if [[ -z "$(sinfo -h -p "$p" -o %P 2>/dev/null)" ]]; then
        echo "Partition '$p' does not exist; available: $(sinfo -h -o %P | tr '\n' ' ')" >&2
        echo "Set PARTITION / DOWNLOAD_PARTITION at the top of slurm/submit_all.sh" >&2
        exit 1
    fi
done

cd "$NEASE_BUILD_SCRIPTS"
eval "$NEASE_BUILD_ACTIVATE"
# <work_dir>/logs (created with the other work folders; fails here if work_dir cannot be created)
LOGS=$(python3 -c "import sys; from common import load_config, paths; print(paths(load_config(sys.argv[1]))['logs'])" "$CONFIG") ||
    { echo "Cannot create the folders in work_dir: set work_dir in $CONFIG to a folder you can write to" >&2; exit 1; }
read_config() { python3 -c "import yaml,sys; c=yaml.safe_load(open('$CONFIG')); print(eval(sys.argv[1]))" "$1"; }
PDB_SHARDS=$(read_config "c['pdb']['shards']")
AF_TASKS=$(read_config "c['alphafold']['tasks']")
AF_ENABLED=$(read_config "c['alphafold']['enabled']")
I3D_ENABLED=$(read_config "(c.get('interactome3d') or {}).get('enabled', False)")

declare -A JOB
# submit <step> <name> <partition> <time> <mem> <cpus> <array or -> <dependency steps (comma separated) or -> command...
submit() {
    local step=$1 name=$2 partition=$3 time=$4 mem=$5 cpus=$6 array=$7 deps=$8
    shift 8
    [[ " $STEPS " == *" $step "* ]] || return 0
    local options=(--parsable --job-name "nease_$name" --partition "$partition" --time "$time" --mem "$mem" --cpus-per-task "$cpus")
    if [[ "$array" != - ]]; then
        options+=(--array "$array" --output "$LOGS/%x_%A_%a.out")
    else
        options+=(--output "$LOGS/%x_%j.out")
    fi
    local after=()
    if [[ "$deps" != - ]]; then
        for d in ${deps//,/ }; do [[ -n "${JOB[$d]:-}" ]] && after+=("${JOB[$d]}"); done
    fi
    (( ${#after[@]} )) && options+=(--dependency "afterok:$(IFS=:; echo "${after[*]}")")
    JOB[$step]=$(sbatch "${options[@]}" slurm/job.sbatch "$@")
    echo "step $step ($name): job ${JOB[$step]}"
}

submit 1  download      "$DOWNLOAD_PARTITION" 12:00:00 8G  2  -  -       python 01_download.py --config "$CONFIG"
submit 2  ensembl       "$PARTITION"          2:00:00  32G 8  -  1       python 02_ensembl.py --config "$CONFIG"
submit 3  pdb_select    "$PARTITION"          1:00:00  16G 1  -  1       python 03_pdb_select.py --config "$CONFIG"
submit 4  pdb_fetch     "$DOWNLOAD_PARTITION" 24:00:00 8G  4  0-9 3     python 04_pdb_fetch.py --config "$CONFIG" --threads 4 --ntasks 10
submit 5  pdb_contacts  "$PARTITION"          4:00:00  16G 1  "0-$((PDB_SHARDS - 1))%100" 4  python 05_pdb_contacts.py --config "$CONFIG"
submit 6  pdb_merge     "$PARTITION"          2:00:00  64G 1  -  5       python 06_pdb_merge.py --config "$CONFIG"
if [[ "$AF_ENABLED" == True ]]; then
submit 7  af_select     "$PARTITION"          2:00:00  16G 1  -  1       python 07_af_select.py --config "$CONFIG"
submit 8  af_fetch      "$DOWNLOAD_PARTITION" 12:00:00 8G  1  "0-$((AF_TASKS - 1))%10" 7  python 08_af_fetch.py --config "$CONFIG"
submit 9  af_contacts   "$PARTITION"          4:00:00  8G  1  "0-$((AF_TASKS - 1))%100" 8  python 09_af_contacts.py --config "$CONFIG"
fi
submit 10 orthologs     "$PARTITION"          1:00:00  16G 16 -  2       python 10_orthologs.py --config "$CONFIG"
submit 11 ppi           "$PARTITION"          4:00:00  64G 1  -  2,10    python 11_ppi.py --config "$CONFIG"
submit 12 elm           "$PARTITION"          2:00:00  32G 1  -  10,11   python 12_elm.py --config "$CONFIG"
if [[ "$I3D_ENABLED" == True ]]; then
# downloads the Interactome3D files unless they are there (python 15_interactome3d.py --config ... --download-only)
submit 15 interactome3d "$DOWNLOAD_PARTITION" 4:00:00  16G 8 -  2       python 15_interactome3d.py --config "$CONFIG" --threads 8
fi
HOMOLOGS_ENABLED=$(read_config "(c.get('homologs') or {}).get('enabled', False)")
if [[ "$HOMOLOGS_ENABLED" == True ]]; then
submit 17 homologs      "$PARTITION"          8:00:00  96G 16 -  2,5     python 17_pdb_homologs.py --config "$CONFIG"
fi
HUMANPPI_ENABLED=$(read_config "(c.get('humanppi') or {}).get('enabled', False)")
if [[ "$HUMANPPI_ENABLED" == True ]]; then
# reads the humanPPI model archive (67 GB) once; or run it elsewhere and copy intermediate/humanppi/ (see README)
submit 18 humanppi      "$PARTITION"          8:00:00  64G 16 -  1       python 18_humanppi.py --config "$CONFIG"
fi
BOLTZ_ENABLED=$(read_config "(c.get('boltz') or {}).get('enabled', False)")
if [[ "$BOLTZ_ENABLED" == True ]]; then
submit 19 boltz         "$PARTITION"          1:00:00  16G 8  -  1,2     python 19_boltz.py --config "$CONFIG"
fi
submit 13 assemble      "$PARTITION"          8:00:00  128G 1 -  6,9,12,15,17,18,19 python 13_assemble.py --config "$CONFIG"
submit 14 validate      "$PARTITION"          2:00:00  32G 1  -  13      python 14_validate.py --config "$CONFIG"
submit 16 losses        "$PARTITION"          2:00:00  32G 1  -  13      python 16_diagnose_losses.py --config "$CONFIG"
