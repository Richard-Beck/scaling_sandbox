#!/bin/bash
#SBATCH --job-name=maree06-area
#SBATCH --array=0-5
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --qos=small
#SBATCH --output=maree06-area-%A_%a.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
make maree2006
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
lambdas=(1e-9 1e-7 1e-6 1e-5 1e-4 1e-3)
lambda="${lambdas[$SLURM_ARRAY_TASK_ID]}"
./maree2006 --steps 200 --seed 601 --no-snapshots \
  --area-mode cpm_lattice_area --lambda-cpm "$lambda" \
  --output "runs/area_mode_sweep/lambda_${lambda}"
