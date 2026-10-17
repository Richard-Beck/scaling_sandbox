#!/bin/bash
#SBATCH --job-name=lfct20x-parse
#SBATCH --partition=red
#SBATCH --qos=normal
#SBATCH --array=1-25
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=reports/lfct_20x_track_summary_assets/logs/parse_%A_%a.out
#SBATCH --error=reports/lfct_20x_track_summary_assets/logs/parse_%A_%a.err

set -euo pipefail
cd /home/4473331/projects/scaling_sandbox
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null
export R_LIBS_USER=/tmp/lfct_empty_R_lib
/app/eb/software/Anaconda3/2024.02-1/envs/R-4.2.3/bin/Rscript scripts/summarize_lfct_trackmate_20x.R
