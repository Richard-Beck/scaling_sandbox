#!/bin/bash
#SBATCH --job-name=mcf102-cpsam
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:a30:1
#SBATCH --exclude=hpctpa3pc0066,hpctpa3pc0067,hpctpa3pc0068,hpctpa3pc0069
#SBATCH --time=02:00:00
#SBATCH --output=reports/lfct_mcf10a_fld102_cpsam_assets/logs/segment_%A_%a.out
#SBATCH --error=reports/lfct_mcf10a_fld102_cpsam_assets/logs/segment_%A_%a.err
set -euo pipefail
cd /home/4473331/projects/scaling_sandbox
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=4
/home/4473331/.conda/envs/cpose/bin/python -u scripts/segment_lfct_mcf10a.py --shard "${SLURM_ARRAY_TASK_ID:-0}" "$@"
