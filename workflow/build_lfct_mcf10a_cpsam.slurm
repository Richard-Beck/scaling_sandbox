#!/bin/bash
#SBATCH --job-name=mcf102-viewer
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=reports/lfct_mcf10a_fld102_cpsam_assets/logs/build_%j.out
#SBATCH --error=reports/lfct_mcf10a_fld102_cpsam_assets/logs/build_%j.err
set -euo pipefail
cd /home/4473331/projects/scaling_sandbox
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=4
/home/4473331/.conda/envs/cpose/bin/python -u scripts/build_lfct_mcf10a_cpsam.py
