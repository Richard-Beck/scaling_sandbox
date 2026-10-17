#!/bin/bash
#SBATCH --job-name=refined_ou_hosny
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --output=wgd_langevin/refined_ou_20261003/results/logs/hosny_%j.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/home/4473331/.conda/envs/scaling_sandbox/bin/python
$PY wgd_langevin/refined_ou_20261003/scripts/hosny_sensitivity/analyze.py
$PY wgd_langevin/refined_ou_20261003/scripts/hosny_sensitivity/refine_surface.py
$PY wgd_langevin/refined_ou_20261003/scripts/hosny_sensitivity/build_figure.py
