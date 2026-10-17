#!/usr/bin/env bash
#SBATCH --job-name=debelly-local-rac
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --qos=small
#SBATCH --array=0-2
#SBATCH --output=polarity1d/results/campaigns/debelly_fixed_substrate_L21_local_rac/logs/%A_%a.out
#SBATCH --error=polarity1d/results/campaigns/debelly_fixed_substrate_L21_local_rac/logs/%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
root=polarity1d/results/campaigns/debelly_fixed_substrate_L21_local_rac
case "${SLURM_ARRAY_TASK_ID}" in
  0) args=(--density 1 --output "$root/density_1");;
  1) args=(--density 10 --output "$root/density_10");;
  2) args=(--density 50 --output "$root/density_50");;
  *) exit 2;;
esac
conda run --no-capture-output -n scaling_sandbox python polarity1d/scripts/run_debelly_fixed_substrate.py --local-rac --step 0.0625 "${args[@]}"
