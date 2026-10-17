#!/bin/bash
#SBATCH --job-name=M7682_default
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=literature_models/M7682/results/default/slurm-%j.log
set -euo pipefail
module load gnuplot/5.4.2-GCCcore-10.3.0
export OMP_NUM_THREADS=1
ROOT="${SLURM_SUBMIT_DIR}/literature_models/M7682"
cd "$ROOT/results/default"
"$ROOT/software/morpheus-2.3.9" --num-threads 1 --file "$ROOT/model/model.xml" > simulation.log 2>&1
