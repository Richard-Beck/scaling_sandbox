#!/usr/bin/env bash
#SBATCH --job-name=figure_fixed_peaks
#SBATCH --array=0-59
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --qos=small
#SBATCH --output=reports/figure_revision/logs/%A_%a.out
#SBATCH --error=reports/figure_revision/logs/%A_%a.err
set -euo pipefail
cd /home/4473331/projects/scaling_sandbox
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
/home/4473331/.conda/envs/scaling_sandbox/bin/python reports/figure_revision/fixed_peaks.py --index "$SLURM_ARRAY_TASK_ID"
