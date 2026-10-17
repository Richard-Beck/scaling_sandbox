#!/bin/bash
#SBATCH --job-name=maree-lambda
#SBATCH --array=0-20
#SBATCH --cpus-per-task=4
#SBATCH --mem=3G
#SBATCH --time=00:30:00
#SBATCH --qos=small
#SBATCH --output=results/lambda_sweep_01/slurm-%A_%a.log
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
/home/4473331/.conda/envs/cpose/bin/python lambda_sweep.py run --index "$SLURM_ARRAY_TASK_ID"
