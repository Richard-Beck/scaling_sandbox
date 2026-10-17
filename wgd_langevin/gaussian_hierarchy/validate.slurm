#!/bin/bash
#SBATCH --job-name=wgd-gauss-check
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=wgd_langevin/gaussian_hierarchy/results/logs/check-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/gaussian_hierarchy/validate.py
