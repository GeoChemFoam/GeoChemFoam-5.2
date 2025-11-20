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

echo "=== Run: interTransferFoam/risingBubble2D ==="
run_task createMesh.sh
run_task initCase0.sh
run_task runCase0.sh
run_task initCaseTransfer.sh
run_task runCaseTransfer.sh
run_task processTransfer.sh
echo "=== Completed: interTransferFoam/risingBubble2D ==="

