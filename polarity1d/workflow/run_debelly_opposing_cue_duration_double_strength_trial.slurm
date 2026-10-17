#!/usr/bin/env bash
#SBATCH --job-name=debelly-opp-2x
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=01:00:00
#SBATCH --qos=small
#SBATCH --output=polarity1d/results/campaigns/debelly_2026_opposing_cue_duration_double_strength_trial/logs/task-%A_%a.out
#SBATCH --error=polarity1d/results/campaigns/debelly_2026_opposing_cue_duration_double_strength_trial/logs/task-%A_%a.err

set -euo pipefail
repository="${SLURM_SUBMIT_DIR}"
cd "${repository}"
manifest="${repository}/polarity1d/campaigns/debelly_2026_opposing_cue_duration_double_strength_trial.json"
output_root="${repository}/polarity1d/results/campaigns/debelly_2026_opposing_cue_duration_double_strength_trial"
mkdir -p "${output_root}/logs"
conda run -n scaling_sandbox python polarity1d/scripts/run_phase1_campaign.py "${manifest}" "${output_root}"
