#!/usr/bin/env bash
#SBATCH --job-name=wp1_track_source
#SBATCH --array=0-63
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:45:00
#SBATCH --qos=small
#SBATCH --output=polarity1d/results/campaigns/phase1_production_v1/logs/%A_%a.out
#SBATCH --error=polarity1d/results/campaigns/phase1_production_v1/logs/%A_%a.err

set -euo pipefail
repository="${SLURM_SUBMIT_DIR}"
cd "${repository}"
output_root="${repository}/polarity1d/results/campaigns/phase1_production_v1"
mkdir -p "${output_root}/logs"
conda run -n scaling_sandbox python polarity1d/scripts/run_phase1_campaign.py \
  polarity1d/campaigns/retained_phase1_tracking_v1.json "${output_root}"
