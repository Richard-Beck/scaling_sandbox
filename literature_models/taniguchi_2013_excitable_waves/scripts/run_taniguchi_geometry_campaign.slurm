#!/bin/bash
#SBATCH --job-name=taniguchi-geom
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=02:00:00
#SBATCH --output=logs/taniguchi-geometry-%A_%a.out

set -euo pipefail
manifest=${1:?Usage: sbatch --array=... scripts/run_taniguchi_geometry_campaign.slurm CONFIG_MANIFEST}
config=$(awk -v target="$((SLURM_ARRAY_TASK_ID + 1))" \
  '!/^[[:space:]]*($|#)/ {seen++; if (seen == target) {print; exit}}' "$manifest")
if [[ -z "$config" ]]; then
  echo "No configuration at index $SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi
conda run -n scaling_sandbox python scripts/run_taniguchi_geometry_campaign.py "$config"
