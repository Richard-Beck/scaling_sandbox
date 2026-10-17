#!/usr/bin/env bash
#SBATCH --job-name=figure_revision
#SBATCH --array=0-121
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:40:00
#SBATCH --qos=small
#SBATCH --output=reports/figure_revision/logs/%A_%a.out
#SBATCH --error=reports/figure_revision/logs/%A_%a.err
set -euo pipefail
cd /home/4473331/projects/scaling_sandbox
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
/home/4473331/.conda/envs/scaling_sandbox/bin/python reports/figure_revision/run_assays.py --index "$SLURM_ARRAY_TASK_ID"
