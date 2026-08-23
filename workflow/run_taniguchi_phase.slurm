#!/bin/bash
#SBATCH --job-name=taniguchi-phase
#SBATCH --array=0-3
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --qos=normal
#SBATCH --output=reports/output/taniguchi/logs/%x-%A_%a.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p reports/output/taniguchi/logs
conditions=(4F 4G 4H 4I)
conda run -n scaling_sandbox python scripts/run_taniguchi.py \
  --model=phase --condition="${conditions[SLURM_ARRAY_TASK_ID]}" --end-time=20
