#!/bin/bash
#SBATCH --job-name=fig4-amplitude
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:20:00
#SBATCH --array=0-259
#SBATCH --output=reports/figure4_heading_amplitude_20261023_%A_%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /home/4473331/.conda/envs/scaling_sandbox/bin/python -B polarity1d/manuscript/code/mechanisms/figure4_round5_analysis_20261023.py heading-worker --index "$SLURM_ARRAY_TASK_ID"
