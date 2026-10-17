#!/bin/bash
#SBATCH --job-name=refined-ou
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=wgd_langevin/refined_ou_20261003/results/logs/fit-%A_%a.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/refined_ou_20261003/scripts/fit.py "$SLURM_ARRAY_TASK_ID"
