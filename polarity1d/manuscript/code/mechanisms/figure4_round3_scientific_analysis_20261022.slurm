#!/bin/bash
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=03:00:00
#SBATCH --output=/home/4473331/projects/scaling_sandbox/polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure4/supporting/metadata/scientific_analysis_%j.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /home/4473331/projects/scaling_sandbox
exec /home/4473331/.conda/envs/scaling_sandbox/bin/python -B polarity1d/manuscript/code/mechanisms/figure4_round3_analysis_20261022.py --source round2 --nulls all
