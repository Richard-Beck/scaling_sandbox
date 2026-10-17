#!/bin/bash
#SBATCH --job-name=giese-diff
#SBATCH --array=0-104
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --qos=small
#SBATCH --output=reports/output/giese_wp/diffusion_screen/logs/%x-%A_%a.out

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
# Accept submission from the repository root or this model bundle.
if [[ -d "literature_models/giese_2015_bulk_surface_wave_pinning" ]]; then
  cd "literature_models/giese_2015_bulk_surface_wave_pinning"
fi
mkdir -p reports/output/giese_wp/diffusion_screen/logs
conda run -n scaling_sandbox python scripts/run_giese_diffusion_screen.py \
  --task-index="${SLURM_ARRAY_TASK_ID}"
