#!/usr/bin/env bash
#SBATCH --job-name=debelly-rev-10ks
#SBATCH --array=0-31
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=01:00:00
#SBATCH --qos=short
#SBATCH --output=polarity1d/results/tmp/trial/debelly_2026_long_window_reversal_trial_v1/logs/task-%A_%a.out
#SBATCH --error=polarity1d/results/tmp/trial/debelly_2026_long_window_reversal_trial_v1/logs/task-%A_%a.err

set -euo pipefail
repository="${SLURM_SUBMIT_DIR}"
cd "${repository}"
manifest="${repository}/polarity1d/campaigns/debelly_2026_long_window_reversal_trial_v1.json"
output_root="${repository}/polarity1d/results/tmp/trial/debelly_2026_long_window_reversal_trial_v1"
mkdir -p "${output_root}/logs"

conda run -n scaling_sandbox python polarity1d/scripts/run_phase1_campaign.py \
  "${manifest}" "${output_root}"
