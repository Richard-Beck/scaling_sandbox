#!/bin/bash
#SBATCH --job-name=jilkine-long-reversal
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=01:00:00
#SBATCH --output=polarity1d/results/campaigns/jilkine_reversal_long_conditioning_v1/logs/task-%A_%a.out

set -euo pipefail

repository="${SLURM_SUBMIT_DIR}"
manifest="${repository}/polarity1d/campaigns/jilkine_reversal_long_conditioning_v1.json"
output_root="${repository}/polarity1d/results/campaigns/jilkine_reversal_long_conditioning_v1"
mkdir -p "${output_root}/logs"
cd "${repository}"
conda run -n scaling_sandbox python polarity1d/scripts/run_phase1_campaign.py \
  "${manifest}" "${output_root}"
