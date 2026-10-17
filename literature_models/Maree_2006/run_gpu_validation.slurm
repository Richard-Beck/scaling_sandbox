#!/bin/bash
#SBATCH --job-name=maree06-gpuval
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --gres=gpu:a30:1
#SBATCH --exclude=hpctpa3pc0066,hpctpa3pc0067,hpctpa3pc0068,hpctpa3pc0069
#SBATCH --time=00:15:00
#SBATCH --qos=small
#SBATCH --output=maree06-gpuval-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
module load cuda12.2/toolkit/12.2.2
make gpu maree2006
./maree2006_gpu --steps 100 --seed 702 --frozen --snapshot \
  --output runs/gpu_validation/gpu_frozen_100
export OMP_NUM_THREADS=4
./maree2006 --steps 100 --seed 702 --frozen --snapshot-seconds 10 \
  --output runs/gpu_validation/cpu_frozen_100
python3 compare_runs.py runs/gpu_validation/cpu_frozen_100 \
  runs/gpu_validation/gpu_frozen_100
