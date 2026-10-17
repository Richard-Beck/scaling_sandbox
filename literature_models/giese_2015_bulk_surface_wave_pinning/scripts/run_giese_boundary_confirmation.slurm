#!/bin/bash
#SBATCH --job-name=giese-confirm
#SBATCH --qos=small
#SBATCH --time=00:08:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --array=12-13
#SBATCH --output=reports/output/giese_wp/boundary_confirmation/logs/%A_%a.out
#SBATCH --error=reports/output/giese_wp/boundary_confirmation/logs/%A_%a.err

set -euo pipefail
cd /home/4473331/projects/scaling_sandbox
mkdir -p reports/output/giese_wp/boundary_confirmation/logs
conda run -n scaling_sandbox python scripts/run_giese_boundary_pilot.py \
  --task-index="${SLURM_ARRAY_TASK_ID}" --end-time=2000 \
  --output=reports/output/giese_wp/boundary_confirmation
