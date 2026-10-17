#!/bin/bash
#SBATCH --job-name=wgd-ordered-times
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=04:00:00
#SBATCH --output=wgd_langevin/persistence_ordered/results/logs/fit-%A_%a.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
fit_index=$((SLURM_ARRAY_TASK_ID + ${ORDERED_OFFSET:-0}))
/home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/persistence_ordered/campaign.py "$fit_index"
