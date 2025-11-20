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

echo "=== Run: interReactiveTransportFoam/2DGrainstonesMultiphaseReactiveTransport ==="
run_task createMesh.sh
run_task runSnappyHexMesh.sh
run_task initCaseRT.sh
run_task runCaseRT.sh
run_task processRT.sh
echo "=== Completed: interReactiveTransportFoam/2DGrainstonesMultiphaseReactiveTransport   ==="

