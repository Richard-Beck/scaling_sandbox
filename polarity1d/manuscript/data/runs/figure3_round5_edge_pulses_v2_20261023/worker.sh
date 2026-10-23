#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /home/4473331/.conda/envs/scaling_sandbox/bin/python -B /home/4473331/projects/scaling_sandbox/polarity1d/manuscript/data/runs/figure3_round5_edge_pulses_v2_20261023/frozen_sources/figure3_round5_edge_20261023.py worker --slot "$SLURM_ARRAY_TASK_ID"
