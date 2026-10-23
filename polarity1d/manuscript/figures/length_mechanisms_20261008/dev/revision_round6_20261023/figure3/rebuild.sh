#!/bin/bash
set -euo pipefail
TASK_FIGURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TASK_PYTHON=/home/4473331/.conda/envs/scaling_sandbox/bin/python
if [[ "${1:-}" == --native ]]; then
  "$TASK_PYTHON" -B "$TASK_FIGURE_DIR/scripts/native_analysis.py" prepare
  "$TASK_PYTHON" -B "$TASK_FIGURE_DIR/scripts/native_analysis.py" submit
  echo 'Native workers submitted; run ./rebuild.sh after completion.'
  exit 0
fi
exec "$TASK_PYTHON" -B "$TASK_FIGURE_DIR/scripts/build.py"
