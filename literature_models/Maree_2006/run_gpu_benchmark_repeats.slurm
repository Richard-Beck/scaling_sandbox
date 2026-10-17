#!/bin/bash
#SBATCH --job-name=maree06-gpurep
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --gres=gpu:a30:1
#SBATCH --exclude=hpctpa3pc0066,hpctpa3pc0067,hpctpa3pc0068,hpctpa3pc0069
#SBATCH --time=00:15:00
#SBATCH --qos=small
#SBATCH --output=maree06-gpurep-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
module load cuda12.2/toolkit/12.2.2
make gpu maree2006
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
for replicate in 1 2 3; do
  /usr/bin/time -f "gpu_rep=${replicate} wall_s=%e user_s=%U rss_kb=%M" \
    ./maree2006_gpu --steps 100 --seed "$((800+replicate))" \
    --output "runs/gpu_benchmark/repeats/gpu_${replicate}"
done
export OMP_NUM_THREADS=4
for replicate in 1 2 3; do
  /usr/bin/time -f "cpu_rep=${replicate} wall_s=%e user_s=%U rss_kb=%M" \
    ./maree2006 --steps 100 --seed "$((800+replicate))" --no-snapshots \
    --output "runs/gpu_benchmark/repeats/cpu_${replicate}"
done
