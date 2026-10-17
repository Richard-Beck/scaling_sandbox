#!/bin/bash
#SBATCH --job-name=taniguchi-plots
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/taniguchi-geometry-plots-%j.out

set -euo pipefail
conda run -n scaling_sandbox python scripts/plot_taniguchi_geometry_campaign.py
