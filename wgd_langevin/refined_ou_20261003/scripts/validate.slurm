#!/bin/bash
#SBATCH --job-name=refined-ou-check
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=wgd_langevin/refined_ou_20261003/results/logs/validate-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/refined_ou_20261003/scripts/validate.py
