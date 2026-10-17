#!/bin/bash
#SBATCH --job-name=debelly-cue-trial
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:45:00
#SBATCH --output=polarity1d/results/tmp/trial/debelly_2026_single_cue_kinetics_trial_v1/logs/task-%A_%a.out

set -euo pipefail

repository="${SLURM_SUBMIT_DIR}"
manifest="${repository}/polarity1d/campaigns/debelly_2026_single_cue_kinetics_trial_v1.json"
output_root="${repository}/polarity1d/results/tmp/trial/debelly_2026_single_cue_kinetics_trial_v1"

mkdir -p "${output_root}/logs"
cd "${repository}"
conda run -n scaling_sandbox python polarity1d/scripts/run_phase1_campaign.py \
  "${manifest}" "${output_root}"
