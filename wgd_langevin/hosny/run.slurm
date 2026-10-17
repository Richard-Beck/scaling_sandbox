#!/bin/bash
#SBATCH --job-name=hosny-langevin
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=wgd_langevin/hosny/logs/%j.out
#SBATCH --error=wgd_langevin/hosny/logs/%j.err
set -euo pipefail
cd /home/4473331/projects/scaling_sandbox
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1
/home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/hosny/run.py --workers 4
