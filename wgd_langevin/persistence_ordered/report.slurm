#!/bin/bash
#SBATCH --job-name=wgd-ordered-report
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=09:00:00
#SBATCH --signal=B:TERM@60
#SBATCH --output=wgd_langevin/persistence_ordered/results/logs/report-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec /home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/persistence_ordered/publish.py --watch
