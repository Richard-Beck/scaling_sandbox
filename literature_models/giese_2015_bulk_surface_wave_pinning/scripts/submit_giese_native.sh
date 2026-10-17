#!/bin/bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
mkdir -p reports/output/giese_wp/logs
sbatch scripts/run_giese_native.slurm
