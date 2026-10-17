#!/bin/bash
#SBATCH --job-name=lfct-smoke
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:a30:1
#SBATCH --exclude=hpctpa3pc0066,hpctpa3pc0067,hpctpa3pc0068,hpctpa3pc0069
#SBATCH --time=00:30:00
#SBATCH --output=lfct/results/logs/smoke_%j.out
#SBATCH --error=lfct/results/logs/smoke_%j.err
set -euo pipefail
cd /home/4473331/projects/scaling_sandbox
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=4
# One GPU allocation for this tiny test; process four fields serially.
for lfct_index in 0 1 2 3; do
  /home/4473331/.conda/envs/cpose/bin/python -u lfct/scripts/pipeline.py --profile smoke --index "$lfct_index"
done
