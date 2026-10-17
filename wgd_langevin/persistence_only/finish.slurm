#!/bin/bash
#SBATCH --job-name=wgd-persistence-report
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:15:00
#SBATCH --output=wgd_langevin/persistence_only/results/logs/report-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/persistence_only/figures.py
/home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/persistence_only/audit.py --require-complete
