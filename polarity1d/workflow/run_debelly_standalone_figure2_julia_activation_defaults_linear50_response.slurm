#!/usr/bin/env bash
#SBATCH --job-name=debelly-f2-linear50
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=01:00:00
#SBATCH --qos=medium
#SBATCH --output=polarity1d/results/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_linear50_response/logs/task-%A_%a.out
#SBATCH --error=polarity1d/results/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_linear50_response/logs/task-%A_%a.err
set -euo pipefail
repository="${SLURM_SUBMIT_DIR}"
cd "${repository}"
manifest="${repository}/polarity1d/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_linear50_response.json"
output_root="${repository}/polarity1d/results/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_linear50_response"
mkdir -p "${output_root}/logs"
conda run -n scaling_sandbox python polarity1d/scripts/run_phase1_campaign.py "${manifest}" "${output_root}"
