#!/bin/bash
#SBATCH --job-name=lfct20x-report
#SBATCH --partition=red
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=reports/lfct_20x_track_summary_assets/logs/aggregate_%j.out
#SBATCH --error=reports/lfct_20x_track_summary_assets/logs/aggregate_%j.err

set -euo pipefail
cd /home/4473331/projects/scaling_sandbox
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null
export R_LIBS_USER=/tmp/lfct_empty_R_lib
export PATH=/app/eb/software/Anaconda3/2024.02-1/envs/R-4.2.3/bin:$PATH
export LFCT_AGGREGATE=1
/app/eb/software/Anaconda3/2024.02-1/envs/R-4.2.3/bin/Rscript scripts/summarize_lfct_trackmate_20x.R
/app/eb/software/Anaconda3/2024.02-1/envs/R-4.2.3/bin/Rscript -e 'rmarkdown::render("reports/lfct_20x_track_summary.Rmd", output_file = "lfct_20x_track_summary.html", quiet = TRUE)'
