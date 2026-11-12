#!/bin/bash
# Generic job of the NEASE network build: runs the command given as arguments in the build environment.
# Submitted by submit_all.sh (which sets the resources, the array, the dependencies and the log file).

set -euo pipefail
# activate the environment (set NEASE_BUILD_ACTIVATE in submit_all.sh, e.g. "module load miniconda; conda activate nease-network-build")
eval "${NEASE_BUILD_ACTIVATE:-true}"
cd "${NEASE_BUILD_SCRIPTS:?}"
echo "$(date) $(hostname) task ${SLURM_ARRAY_TASK_ID:-}: $*"
"$@"
echo "$(date) done"
