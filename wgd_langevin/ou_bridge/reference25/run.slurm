#!/bin/bash
#SBATCH --job-name=ou-reference25
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=02:00:00
#SBATCH --output=wgd_langevin/ou_bridge/reference25/results/logs/%A_%a.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg
/home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/ou_bridge/reference25/campaign.py "$1" "${SLURM_ARRAY_TASK_ID:-0}"
