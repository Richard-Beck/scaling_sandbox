#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
model_root="$repo_root/src/debelly"

usage() {
  echo "Usage: $0 {instantiate|baseline|mechanochemical|local-inhibition|axisymmetric} [smoke|published]"
}

case "${1:-}" in
  instantiate) case_name=instantiate ;;
  baseline|mechanochemical|local-inhibition|axisymmetric) case_name="$1" ;;
  *) usage; exit 2 ;;
esac
profile="${2:-smoke}"
case "$profile" in smoke|published) ;; *) usage; exit 2 ;; esac

if [[ -n "${JULIA:-}" ]]; then
  julia_cmd="$JULIA"
elif command -v julia >/dev/null 2>&1; then
  julia_cmd="$(command -v julia)"
else
  echo "Julia is not on PATH. Load Julia 1.8+ or set JULIA=/path/to/julia." >&2
  exit 127
fi

if [[ "$case_name" == instantiate ]]; then
  exec "$julia_cmd" "$model_root/instantiate.jl"
fi

if [[ "$profile" == smoke ]]; then
  export DEBELLY_WRITE_OUTPUTS="${DEBELLY_WRITE_OUTPUTS:-0}"
  export DEBELLY_PARTITION="${DEBELLY_PARTITION:-20}"
  export DEBELLY_AXIS_N="${DEBELLY_AXIS_N:-40}"
  export DEBELLY_OPTO_TIME="${DEBELLY_OPTO_TIME:-2}"
  if [[ "$case_name" == axisymmetric ]]; then
    export DEBELLY_T="${DEBELLY_T:-20}"
    export DEBELLY_OPTO_START="${DEBELLY_OPTO_START:-2}"
    export DEBELLY_OPTO_STOP="${DEBELLY_OPTO_STOP:-12}"
  else
    export DEBELLY_T="${DEBELLY_T:-6}"
  fi
else
  export DEBELLY_WRITE_OUTPUTS="${DEBELLY_WRITE_OUTPUTS:-1}"
  if [[ "$case_name" == axisymmetric ]]; then
    export DEBELLY_AXIS_N="${DEBELLY_AXIS_N:-80}"
    export DEBELLY_T="${DEBELLY_T:-350}"
  fi
fi

case "$case_name" in
  baseline)
    project="$model_root/one_d"
    export DEBELLY_ALPHA_OPTO=0 DEBELLY_BETA_OPTO=0
    work_dir="$model_root/one_d"
    source_file="Mechanochemical_general_code.jl"
    ;;
  mechanochemical)
    project="$model_root/one_d"
    work_dir="$model_root/one_d/protocols"
    source_file="mechanochemical_opto_front2back.jl"
    ;;
  local-inhibition)
    project="$model_root/one_d"
    work_dir="$model_root/one_d/protocols"
    source_file="localinhibition_opto_front2back.jl"
    ;;
  axisymmetric)
    project="$model_root/axisymmetric"
    work_dir="$model_root/axisymmetric"
    source_file="examples/SurfaceViscousFlows.jl"
    ;;
esac

echo "Running De Belly case=$case_name profile=$profile"
echo "Julia: $julia_cmd"
start_seconds=$SECONDS
(
  cd "$work_dir"
  "$julia_cmd" --project="$project" "$source_file"
)
elapsed=$((SECONDS - start_seconds))
echo "DEBELLY_RUNNER case=$case_name profile=$profile wall_seconds=$elapsed"
