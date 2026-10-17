#!/bin/bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
mkdir -p reports/output/giese_wp/diffusion_screen/logs
array_job=$(sbatch --parsable scripts/run_giese_diffusion_screen.slurm)
plot_job=$(sbatch --parsable --dependency="afterok:${array_job}" \
  scripts/collect_giese_diffusion_screen.slurm)
printf 'simulation_array_job=%s\nplot_job=%s\n' "${array_job}" "${plot_job}"
