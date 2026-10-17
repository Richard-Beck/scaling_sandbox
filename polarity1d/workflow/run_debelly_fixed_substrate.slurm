#!/usr/bin/env bash
#SBATCH --job-name=debelly-clutch
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --qos=small
#SBATCH --array=0-7
#SBATCH --output=polarity1d/results/campaigns/debelly_fixed_substrate_L21/logs/%A_%a.out
#SBATCH --error=polarity1d/results/campaigns/debelly_fixed_substrate_L21/logs/%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
root=polarity1d/results/campaigns/debelly_fixed_substrate_L21
case "${SLURM_ARRAY_TASK_ID}" in
  0) args=(--density 1 --output "$root/density_1");;
  1) args=(--density 10 --output "$root/density_10");;
  2) args=(--density 50 --output "$root/density_50");;
  3) args=(--end 1000 --output "$root/check_base");;
  4) args=(--end 1000 --step 0.0625 --output "$root/check_timestep");;
  5) args=(--end 1000 --spacing 0.0125 --output "$root/check_substrate");;
  6) args=(--end 1000 --cells 210 --spacing 0.0125 --output "$root/check_cell_grid");;
  7) args=(--step 0.0625 --output "$root/check_long_timestep");;
  *) exit 2;;
esac
conda run --no-capture-output -n scaling_sandbox python polarity1d/scripts/run_debelly_fixed_substrate.py "${args[@]}"
