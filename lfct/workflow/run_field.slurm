#!/bin/bash
#SBATCH --job-name=lfct-field
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:a30:1
#SBATCH --exclude=hpctpa3pc0066,hpctpa3pc0067,hpctpa3pc0068,hpctpa3pc0069
#SBATCH --time=01:00:00
#SBATCH --output=lfct/results/logs/field_%A_%a.out
#SBATCH --error=lfct/results/logs/field_%A_%a.err
set -euo pipefail
cd /home/4473331/projects/scaling_sandbox
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=4
/home/4473331/.conda/envs/cpose/bin/python -u lfct/scripts/pipeline.py --index "${SLURM_ARRAY_TASK_ID:-0}" "$@"
