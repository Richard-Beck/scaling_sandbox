#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /home/4473331/.conda/envs/scaling_sandbox/bin/python -B /home/4473331/projects/scaling_sandbox/polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round6_20261023/figure3/scripts/native_analysis.py worker --slot "$SLURM_ARRAY_TASK_ID"
