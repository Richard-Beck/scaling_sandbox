#!/usr/bin/env bash
#SBATCH --job-name=wp1_fig1_examples
#SBATCH --array=0-11
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=02:00:00
#SBATCH --qos=short
#SBATCH --output=polarity1d/results/campaigns/wp1_figure1_examples_v1/logs/%A_%a.out
#SBATCH --error=polarity1d/results/campaigns/wp1_figure1_examples_v1/logs/%A_%a.err

set -euo pipefail

cd /home/4473331/projects/scaling_sandbox
mkdir -p polarity1d/results/campaigns/wp1_figure1_examples_v1/logs

config=$(find polarity1d/results/campaigns/wp1_figure1_examples_v1/configs -name "$(printf '%02d' "$SLURM_ARRAY_TASK_ID")_*.json" -print -quit)
assay=$(basename "$config" .json | cut -d_ -f2)
model=$(basename "$config" .json | cut -d_ -f3-)
output="polarity1d/results/campaigns/wp1_figure1_examples_v1/tasks/${assay}_${model}"

conda run -n scaling_sandbox python polarity1d/run.py "$config" "$output"
