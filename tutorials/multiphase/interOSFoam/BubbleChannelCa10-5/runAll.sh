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

echo "=== Run: interOSFoam/BubbleChannelCa10-5 ==="
run_task createMesh.sh
run_task initCase0.sh
run_task runCase0.sh
run_task initCaseTPFlow.sh
run_task runCaseTPFlow.sh
run_task processCase.sh
echo "=== Completed: interOSFoam/BubbleChannelCa10-5 ==="

