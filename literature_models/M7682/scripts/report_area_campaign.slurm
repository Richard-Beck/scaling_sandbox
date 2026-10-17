#!/bin/bash
#SBATCH --job-name=M7682-report
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=literature_models/M7682/results/area_scaling_80k/report-%j.log
set -euo pipefail
/home/4473331/.conda/envs/scaling_sandbox/bin/python literature_models/M7682/scripts/report_area_campaign.py
