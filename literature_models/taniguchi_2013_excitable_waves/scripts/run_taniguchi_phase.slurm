#!/bin/bash
#SBATCH --job-name=taniguchi-phase
#SBATCH --array=0-3
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --qos=short
#SBATCH --output=reports/output/taniguchi/logs/%x-%A_%a.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
# Accept submission from the repository root or this model bundle.
if [[ -d "literature_models/taniguchi_2013_excitable_waves" ]]; then
  cd "literature_models/taniguchi_2013_excitable_waves"
fi
mkdir -p reports/output/taniguchi/logs
conditions=(4F 4G 4H 4I)
seed="${TANIGUCHI_SEED:-1}"
output="${TANIGUCHI_OUTPUT:-reports/output/taniguchi}"
dt="${TANIGUCHI_DT:-0.0005}"
backend="${TANIGUCHI_BACKEND:-fast}"
output_interval="${TANIGUCHI_OUTPUT_INTERVAL:-0.05}"
conda run -n scaling_sandbox python scripts/run_taniguchi.py \
  --model=phase --condition="${conditions[SLURM_ARRAY_TASK_ID]}" --end-time=20 \
  --seed="${seed}" --dt="${dt}" --backend="${backend}" \
  --output-interval="${output_interval}" --output="${output}"
