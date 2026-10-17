#!/bin/bash
#SBATCH --job-name=xuan-memory-ordered
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=wgd_langevin/persistence_only/results/xuan_memory_ordered/fit-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/persistence_only/xuan_memory_ordered.py
