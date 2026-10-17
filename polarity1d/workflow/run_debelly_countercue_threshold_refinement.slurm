#!/usr/bin/env bash
#SBATCH --job-name=debelly-f2-bisect
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=01:00:00
#SBATCH --qos=small
#SBATCH --output=polarity1d/results/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_threshold_refinement_logs/task-%A_%a.out
#SBATCH --error=polarity1d/results/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_threshold_refinement_logs/task-%A_%a.err
set -euo pipefail
repository="${SLURM_SUBMIT_DIR}"
round_number="${REFINEMENT_ROUND:?REFINEMENT_ROUND must be set}"
name="debelly_2026_standalone_figure2_julia_activation_defaults_threshold_refinement_r${round_number}"
cd "${repository}"
manifest="${repository}/polarity1d/campaigns/${name}.json"
output_root="${repository}/polarity1d/results/campaigns/${name}"
mkdir -p "${output_root}/logs"
conda run -n scaling_sandbox python polarity1d/scripts/run_phase1_campaign.py "${manifest}" "${output_root}"
