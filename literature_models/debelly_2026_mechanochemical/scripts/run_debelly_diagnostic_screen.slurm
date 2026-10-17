#!/bin/bash
#SBATCH --job-name=debelly-diagnostic
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:30:00
#SBATCH --output=logs/debelly-diagnostic-%A_%a.out

set -euo pipefail
mkdir -p logs
conda run -n scaling_sandbox python scripts/run_debelly_regime_screen.py \
  --diagnostic-task "$SLURM_ARRAY_TASK_ID"
