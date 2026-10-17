#!/bin/bash
#SBATCH --job-name=maree-report
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --qos=small
#SBATCH --output=results/validation_01/report-%j.log
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
/home/4473331/.conda/envs/cpose/bin/python validation_campaign.py report
