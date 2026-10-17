#!/bin/bash
#SBATCH --job-name=maree06-gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --gres=gpu:a30:1
#SBATCH --exclude=hpctpa3pc0066,hpctpa3pc0067,hpctpa3pc0068,hpctpa3pc0069
#SBATCH --time=00:30:00
#SBATCH --qos=small
#SBATCH --output=maree06-gpu-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
module load cuda12.2/toolkit/12.2.2
make gpu maree2006
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
/usr/bin/time -f 'gpu_wall_s=%e gpu_user_s=%U gpu_rss_kb=%M' \
  ./maree2006_gpu --steps 100 --seed 701 --output runs/gpu_benchmark/gpu_100
export OMP_NUM_THREADS=4
/usr/bin/time -f 'cpu_wall_s=%e cpu_user_s=%U cpu_rss_kb=%M' \
  ./maree2006 --steps 100 --seed 701 --output runs/gpu_benchmark/cpu_100 --no-snapshots
