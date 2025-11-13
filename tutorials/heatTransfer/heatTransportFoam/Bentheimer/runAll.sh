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

echo "=== Run: heatTransportFoam/Bentheimer ==="
run_task createMesh.sh
run_task initCaseFlow.sh
run_task runCaseFlow.sh
run_task processFlow.sh
run_task initCaseHeat.sh
run_task runCaseHeat.sh
run_task processHeat.sh
echo "=== Completed: heatTransportFoam/Bentheimer ==="

