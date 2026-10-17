#!/bin/bash
#SBATCH --job-name=polarity1d-db-titrate
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=polarity1d/results/campaigns/debelly_input_titration_v1/logs/task-%A_%a.out

set -euo pipefail

repository="${SLURM_SUBMIT_DIR}"
manifest="${POLARITY1D_MANIFEST:-${repository}/polarity1d/campaigns/debelly_input_titration_v1.json}"
output_root="${POLARITY1D_OUTPUT_ROOT:-${repository}/polarity1d/results/campaigns/debelly_input_titration_v1}"

mkdir -p "${output_root}/logs"
cd "${repository}"
force_arguments=()
if [[ "${POLARITY1D_FORCE:-0}" == "1" ]]; then
  force_arguments+=(--force)
fi
conda run -n scaling_sandbox python polarity1d/scripts/run_phase1_campaign.py \
  "${manifest}" "${output_root}" "${force_arguments[@]}"
