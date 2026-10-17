#!/bin/bash
#SBATCH --job-name=maree-bias-report
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:10:00
#SBATCH --qos=small
#SBATCH --output=results/lambda_bias_01/report-%j.log
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
/home/4473331/.conda/envs/cpose/bin/python lambda_bias.py report
