#!/usr/bin/env bash

set -Eeuo pipefail
cd "$(dirname "$0")"

run_task() {
  local name="$1"
  if [[ -f "$name" ]]; then
    [[ -x "$name" ]] || chmod +x "$name" || true
    echo "--- Running: $name ---"
    "./$name"
  else
    echo "--- Running command: $name ---"
    $name
  fi
  echo "--- Done: $name ---"
}

echo "=== Run: solidDisplacementElasticVOSFoam/sphereVOS ==="
run_task createMesh.sh
#run_task runSubsetMesh.sh
run_task initCaseStress.sh
run_task runCaseStress.sh
run_task processStress.sh
echo "=== Completed: solidDisplacementElasticVOSFoam/sphereVOS ==="

