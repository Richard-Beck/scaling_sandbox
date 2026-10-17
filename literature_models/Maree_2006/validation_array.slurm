#!/bin/bash
#SBATCH --job-name=maree-validation
#SBATCH --array=0-20
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=02:00:00
#SBATCH --qos=small
#SBATCH --output=results/validation_01/slurm-%A_%a.log
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
/home/4473331/.conda/envs/cpose/bin/python validation_campaign.py run --index "$SLURM_ARRAY_TASK_ID"
