#!/bin/bash
#SBATCH --job-name=taniguchi-render
#SBATCH --array=0-3
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --qos=short
#SBATCH --output=reports/output/taniguchi/logs/%x-%A_%a.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
# Accept submission from the repository root or this model bundle.
if [[ -d "literature_models/taniguchi_2013_excitable_waves" ]]; then
  cd "literature_models/taniguchi_2013_excitable_waves"
fi
mkdir -p reports/output/taniguchi/logs
module load FFmpeg/7.1.2-GCCcore-14.3.0
conditions=(4F 4G 4H 4I)
condition="${conditions[SLURM_ARRAY_TASK_ID]}"
input_dir="${TANIGUCHI_INPUT:-reports/output/taniguchi}"
output_dir="${TANIGUCHI_RENDER_OUTPUT:-reports/output/taniguchi/rendered}"
conda run -n scaling_sandbox python scripts/render_taniguchi.py \
  "${input_dir}/phase_${condition}.npz" --output-dir="${output_dir}" --mode=both
