#!/bin/bash
#SBATCH --job-name=debelly-2026-xval
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=polarity1d/results/cross_validation/debelly_2026_v1/logs/task-%A_%a.out

set -euo pipefail

repository="${SLURM_SUBMIT_DIR}"
cases=(competition_l5 reversal_l15 tracking_l39)
mkdir -p "${repository}/polarity1d/results/cross_validation/debelly_2026_v1/logs"
cd "${repository}"
conda run -n scaling_sandbox python polarity1d/scripts/cross_validate_debelly_2026.py \
  "${cases[${SLURM_ARRAY_TASK_ID}]}"
