#!/usr/bin/env bash
#SBATCH --job-name=wp1_track_audit
#SBATCH --array=0-63
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:45:00
#SBATCH --qos=small
#SBATCH --output=polarity1d/results/campaigns/wp1_tracking_reliability_audit_v1/logs/%A_%a.out
#SBATCH --error=polarity1d/results/campaigns/wp1_tracking_reliability_audit_v1/logs/%A_%a.err

set -euo pipefail
cd /home/4473331/projects/scaling_sandbox
mkdir -p polarity1d/results/campaigns/wp1_tracking_reliability_audit_v1/logs
conda run -n scaling_sandbox python polarity1d/scripts/run_phase1_campaign.py \
  polarity1d/campaigns/wp1_tracking_reliability_audit_v1.json \
  polarity1d/results/campaigns/wp1_tracking_reliability_audit_v1
