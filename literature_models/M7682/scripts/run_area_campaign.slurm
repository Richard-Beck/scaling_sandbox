#!/bin/bash
#SBATCH --job-name=M7682-area
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --array=0-29
#SBATCH --output=literature_models/M7682/results/area_scaling_80k/slurm-%A_%a.log
set -euo pipefail
module load gnuplot/5.4.2-GCCcore-10.3.0
export OMP_NUM_THREADS=1
ROOT="${SLURM_SUBMIT_DIR}/literature_models/M7682"
BASE="$ROOT/results/area_scaling_80k"
TASK="$SLURM_ARRAY_TASK_ID"
GROUP=$(( TASK / 10 ))
SEED=$(( 19 + TASK % 10 ))
case "$GROUP" in
    0) SIZE=area_0.5x ;;
    1) SIZE=area_1x ;;
    2) SIZE=area_2x ;;
esac
cd "$BASE/$SIZE/seed_$SEED"
"$ROOT/software/morpheus-2.3.9" --num-threads 1 --file model.xml > simulation.log 2>&1
test -s celltracks.xml
test -s cells_displacement.dat
test -s plot_80000.png
printf 'complete\n' > COMPLETE
