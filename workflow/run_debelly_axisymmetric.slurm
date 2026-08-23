#!/usr/bin/env bash
#SBATCH --job-name=debelly-axis
#SBATCH --qos=long
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/debelly-axis-%j.out

set -euo pipefail
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
mkdir -p logs

# Override JULIA before submission if Julia 1.8+ is not exposed by the module.
if [[ -z "${JULIA:-}" ]]; then
  module load Julia
fi

scripts/run_debelly.sh axisymmetric published
