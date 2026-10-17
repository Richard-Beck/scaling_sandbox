#!/bin/bash
#SBATCH --job-name=ou-fit
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=wgd_langevin/ou_bridge/results/logs/fit-%A_%a.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/ou_bridge/fit.py "$SLURM_ARRAY_TASK_ID"
