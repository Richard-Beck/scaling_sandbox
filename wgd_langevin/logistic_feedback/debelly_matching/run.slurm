#!/usr/bin/env bash
#SBATCH --job-name=mca-reduction
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:45:00
#SBATCH --qos=small
#SBATCH --output=wgd_langevin/logistic_feedback/debelly_matching/results_self_activation/slurm-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
conda run --no-capture-output -n scaling_sandbox python wgd_langevin/logistic_feedback/debelly_matching/campaign.py "$@"
