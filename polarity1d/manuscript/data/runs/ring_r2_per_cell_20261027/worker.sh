#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /home/4473331/.conda/envs/scaling_sandbox/bin/python -B /home/4473331/projects/scaling_sandbox/polarity1d/manuscript/data/runs/ring_r2_per_cell_20261027/frozen_sources/campaign.py worker --slot "$SLURM_ARRAY_TASK_ID" --level "${R2_LEVEL:-0}"
