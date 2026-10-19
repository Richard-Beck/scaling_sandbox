#!/bin/bash
#SBATCH --job-name=f4_portable_redraw
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:20:00
#SBATCH --output=/home/4473331/projects/scaling_sandbox/polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/portable_export_operation_20261019/%j.out
#SBATCH --error=/home/4473331/projects/scaling_sandbox/polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/portable_export_operation_20261019/%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /home/4473331/projects/scaling_sandbox
/home/4473331/.conda/envs/scaling_sandbox/bin/python /home/4473331/projects/scaling_sandbox/polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/portable_export_operation_20261019/frozen_export.py
/home/4473331/.conda/envs/scaling_sandbox/bin/python /home/4473331/projects/scaling_sandbox/polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/portable_reconstruction_20261019/helpers/redraw.py --out /home/4473331/projects/scaling_sandbox/polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/portable_redraw_check_20261019 --recompute-ci
