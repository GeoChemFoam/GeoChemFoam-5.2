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

echo "=== Run: 2DGrainstoneBimolecular ==="
run_task createMesh.sh
run_task runSnappyHexMesh.sh
run_task initCaseFlow.sh
run_task runCaseFlow.sh
run_task processFlow.sh
run_task initCaseReactiveTransport.sh
run_task runCaseReactiveTransport.sh
run_task processReactiveTransport.sh
echo "=== Completed: 2DGrainstoneBimolecular ==="

