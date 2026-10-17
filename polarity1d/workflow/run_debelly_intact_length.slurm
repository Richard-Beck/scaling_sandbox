#!/bin/bash
#SBATCH --job-name=polarity1d-db-length
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:30:00
#SBATCH --output=polarity1d/results/campaigns/debelly_intact_length_v2/logs/task-%A_%a.out

set -euo pipefail

repository="${SLURM_SUBMIT_DIR}"
manifest="${POLARITY1D_MANIFEST:-${repository}/polarity1d/campaigns/debelly_intact_length_v2.json}"
output_root="${POLARITY1D_OUTPUT_ROOT:-${repository}/polarity1d/results/campaigns/debelly_intact_length_v2}"

mkdir -p "${output_root}/logs"
cd "${repository}"
conda run -n scaling_sandbox python polarity1d/scripts/run_phase1_campaign.py \
  "${manifest}" "${output_root}"
