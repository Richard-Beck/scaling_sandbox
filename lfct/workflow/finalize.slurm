#!/bin/bash
#SBATCH --job-name=lfct-finalize
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=lfct/results/logs/finalize_%j.out
#SBATCH --error=lfct/results/logs/finalize_%j.err
set -euo pipefail
cd /home/4473331/projects/scaling_sandbox
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=2
/home/4473331/.conda/envs/cpose/bin/python -u lfct/scripts/finalize.py
