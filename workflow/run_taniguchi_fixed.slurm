#!/bin/bash
#SBATCH --job-name=taniguchi-fixed
#SBATCH --array=0-3
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --qos=small
#SBATCH --output=reports/output/taniguchi/logs/%x-%A_%a.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p reports/output/taniguchi/logs
conditions=(6.0 5.7 4.2 2.5)
conda run -n scaling_sandbox python scripts/run_taniguchi.py \
  --model=fixed --condition="${conditions[SLURM_ARRAY_TASK_ID]}" --end-time=30
