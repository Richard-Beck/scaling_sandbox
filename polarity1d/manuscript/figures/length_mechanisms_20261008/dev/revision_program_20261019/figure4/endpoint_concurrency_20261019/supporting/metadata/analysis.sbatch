#!/bin/bash
#SBATCH --job-name=f4_endpoint_context
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:30:00
#SBATCH --output=/home/4473331/projects/scaling_sandbox/polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/endpoint_concurrency_20261019/supporting/metadata/analysis_%j.out
#SBATCH --error=/home/4473331/projects/scaling_sandbox/polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/endpoint_concurrency_20261019/supporting/metadata/analysis_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /home/4473331/projects/scaling_sandbox
test -f /home/4473331/projects/scaling_sandbox/polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/signature_assembly_20261019/supporting/metadata/FINAL_POSTPROCESS_RECEIPT.json
/home/4473331/.conda/envs/scaling_sandbox/bin/python /home/4473331/projects/scaling_sandbox/polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/endpoint_concurrency_20261019/supporting/metadata/frozen_endpoint_analysis.py
