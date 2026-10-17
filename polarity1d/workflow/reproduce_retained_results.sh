#!/usr/bin/env bash
# Submit, collect, analyze, and render every retained polarity1d result.

set -euo pipefail

repository="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${repository}"
run=(conda run -n scaling_sandbox python)

"${run[@]}" polarity1d/scripts/create_wp1_figure1_common_horizon_campaign.py
"${run[@]}" polarity1d/scripts/create_wp1_tracking_reliability_audit.py
"${run[@]}" polarity1d/scripts/create_wp1_tracking_supplement_campaign.py
"${run[@]}" polarity1d/scripts/create_otsuji_reversal_geometry_sweep.py
"${run[@]}" polarity1d/scripts/create_retained_phase1_tracking_campaign.py
"${run[@]}" polarity1d/scripts/create_debelly_standalone_figure2_julia_activation_defaults.py base
"${run[@]}" polarity1d/scripts/create_debelly_panel_a_10000.py
"${run[@]}" polarity1d/scripts/create_wp1_figure1_examples.py

job_ids=(
  "$(sbatch --parsable --array=0-199 polarity1d/workflow/run_debelly_standalone_figure2_julia_activation_defaults_base.slurm)"
  "$(sbatch --parsable --array=0-79 polarity1d/workflow/run_debelly_panel_a_10000.slurm)"
  "$(sbatch --parsable polarity1d/workflow/run_wp1_figure1_common_horizon.slurm)"
  "$(sbatch --parsable polarity1d/workflow/run_retained_phase1_tracking.slurm)"
  "$(sbatch --parsable polarity1d/workflow/run_wp1_tracking_reliability_audit.slurm)"
  "$(sbatch --parsable polarity1d/workflow/run_wp1_tracking_supplement.slurm)"
  "$(sbatch --parsable polarity1d/workflow/run_otsuji_reversal_geometry_sweep.slurm)"
  "$(sbatch --parsable polarity1d/workflow/run_wp1_figure1_examples.slurm)"
)
joined="$(IFS=,; echo "${job_ids[*]}")"
while squeue -h -j "${joined}" | grep -q .; do sleep 15; done
if sacct -n -X -j "${joined}" --format=State | awk 'NF && $1 != "COMPLETED" {bad=1} END {exit bad}'; then
  :
else
  echo "At least one first-stage SLURM task failed" >&2
  exit 1
fi

collect() {
  "${run[@]}" polarity1d/scripts/collect_phase1_campaign.py "$1" "$2"
}
collect polarity1d/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_base.json polarity1d/results/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_base
collect polarity1d/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_panel_a_10000.json polarity1d/results/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_panel_a_10000
collect polarity1d/campaigns/wp1_figure1_common_horizon_v1.json polarity1d/results/campaigns/wp1_figure1_common_horizon_v1
collect polarity1d/campaigns/retained_phase1_tracking_v1.json polarity1d/results/campaigns/phase1_production_v1
collect polarity1d/campaigns/wp1_tracking_reliability_audit_v1.json polarity1d/results/campaigns/wp1_tracking_reliability_audit_v1
collect polarity1d/campaigns/wp1_tracking_supplement_v1.json polarity1d/results/campaigns/wp1_tracking_supplement_v1
collect polarity1d/campaigns/otsuji_reversal_geometry_sweep_l15_v1.json polarity1d/results/campaigns/otsuji_reversal_geometry_sweep_l15_v1

"${run[@]}" polarity1d/scripts/export_retained_phase1_tracking.py polarity1d/results/campaigns/phase1_production_v1
"${run[@]}" polarity1d/scripts/analyze_wp1_tracking_reliability_audit.py polarity1d/campaigns/wp1_tracking_reliability_audit_v1.json polarity1d/results/campaigns/wp1_tracking_reliability_audit_v1
"${run[@]}" polarity1d/scripts/analyze_otsuji_reversal_geometry_sweep.py polarity1d/campaigns/otsuji_reversal_geometry_sweep_l15_v1.json polarity1d/results/campaigns/otsuji_reversal_geometry_sweep_l15_v1

"${run[@]}" polarity1d/scripts/create_debelly_standalone_figure2_julia_activation_defaults.py response
response_job="$(sbatch --parsable --wait --array=0-3059 polarity1d/workflow/run_debelly_standalone_figure2_julia_activation_defaults_linear50_response.slurm)"
collect polarity1d/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_linear50_response.json polarity1d/results/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_linear50_response

"${run[@]}" polarity1d/scripts/qualify_models.py
"${run[@]}" polarity1d/scripts/build_retained_wp1_figures.py
"${run[@]}" polarity1d/scripts/build_debelly_standalone_figure2_v2.py --julia-activation-defaults-linear50
"${run[@]}" polarity1d/scripts/build_debelly_figure2_kymograph_supplement.py
"${run[@]}" polarity1d/scripts/package_retained_results.py

echo "Rendered results are ready for the pixel checkpoint in polarity1d/results/figures"
