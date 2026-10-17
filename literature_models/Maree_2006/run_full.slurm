#!/bin/bash
#SBATCH --job-name=maree06
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --qos=small
#SBATCH --output=maree06-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
make maree2006
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
run_dir="${1:-runs/polarized_seed1}"
steps="${2:-12000}"
./maree2006 --steps "$steps" --seed 1 --output "$run_dir" --snapshot-seconds 10
