#!/bin/bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"
mkdir -p reports/output/giese_wp/logs
sbatch workflow/run_giese_native.slurm
