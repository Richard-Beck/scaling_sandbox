#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /home/4473331/.conda/envs/scaling_sandbox/bin/python -B /home/4473331/projects/scaling_sandbox/polarity1d/manuscript/data/runs/poisson_ring_modes_20261027/frozen_sources/production.py worker --run /home/4473331/projects/scaling_sandbox/polarity1d/manuscript/data/runs/poisson_ring_modes_20261027 --slot "$SLURM_ARRAY_TASK_ID" --level "${POISSON_LEVEL:-0}"
