#!/bin/bash
#SBATCH --job-name=hosny-poisson
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --output=wgd_langevin/hosny/poisson_velocity/logs/%j.out
#SBATCH --error=wgd_langevin/hosny/poisson_velocity/logs/%j.err
set -euo pipefail
cd /home/4473331/projects/scaling_sandbox
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/hosny/poisson_velocity/run.py --workers 4
/home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/hosny/poisson_velocity/plot.py
