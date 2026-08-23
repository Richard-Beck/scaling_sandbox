#!/bin/bash
#SBATCH --job-name=giese-wp
#SBATCH --array=0-13
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --qos=small
#SBATCH --output=reports/output/giese_wp/logs/%x-%A_%a.out

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
mkdir -p reports/output/giese_wp/logs
conda run -n scaling_sandbox python scripts/run_giese_native.py \
  --task-index="${SLURM_ARRAY_TASK_ID}"

