#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /home/4473331/.conda/envs/scaling_sandbox/bin/python -B /home/4473331/projects/scaling_sandbox/polarity1d/manuscript/data/runs/signature_holmes_retention_round2_20261022/worker.py worker --out /home/4473331/projects/scaling_sandbox/polarity1d/manuscript/data/runs/signature_holmes_retention_round2_20261022 --task "$SLURM_ARRAY_TASK_ID"
