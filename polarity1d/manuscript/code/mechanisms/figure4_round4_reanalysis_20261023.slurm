#!/bin/bash
#SBATCH --job-name=fig4-event-rescoring
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/home/4473331/projects/scaling_sandbox/reports/figure4_event_scoring_20261023_%j.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /home/4473331/projects/scaling_sandbox
exec /home/4473331/.conda/envs/scaling_sandbox/bin/python -B polarity1d/manuscript/code/mechanisms/figure4_round4_reanalysis_20261023.py --reanalyze --analysis-only
