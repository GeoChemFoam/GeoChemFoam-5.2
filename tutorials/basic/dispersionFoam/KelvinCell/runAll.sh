#!/usr/bin/env bash

set -Eeuo pipefail
IFS=$'\n\t'

# Always run from the script's directory
cd "$(dirname "$0")"

echo "=== GeoChemFoam KelvinCell: run all tutorial steps ==="

run_step() {
  local script="$1"
  if [[ ! -f "$script" ]]; then
    echo "[ERROR] Missing script: $script" >&2
    exit 1
  fi
  if [[ ! -x "$script" ]]; then
    # Try to make it executable if it isn't already
    chmod +x "$script" || true
  fi
  echo "--- Running: $script ---"
  "./$script"
  echo "--- Done: $script ---"
}

# Order derived from the KelvinCell tutorial workflow:
# 1) Create mesh (and optional local refinement/smoothing in createMesh.sh)
# 2) Initialize and run flow; post-process flow
# 3) Initialize and run dispersion; post-process dispersion

run_step createMesh.sh

# If you want extra smoothing outside createMesh.sh's optional refinement, uncomment:
# run_step runSmoothSolidSurface.sh

run_step initCaseFlow.sh
run_step runCaseFlow.sh
run_step processFlow.sh

run_step initCaseDispersion.sh
run_step runCaseDispersion.sh
run_step processDispersion.sh

echo "=== All steps completed successfully ==="

