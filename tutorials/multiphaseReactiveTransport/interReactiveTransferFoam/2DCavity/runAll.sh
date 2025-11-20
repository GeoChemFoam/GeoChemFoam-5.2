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

echo "=== Run: interReactiveTransferFoam/2DCacity ==="
run_task createMesh.sh
run_task initCase0.sh
run_task runCase0.sh
run_task initCaseRT.sh
run_task runCaseRT.sh
run_task processRT.sh
echo "=== Completed: interReactiveTransferFoam/2DCacity  ==="

