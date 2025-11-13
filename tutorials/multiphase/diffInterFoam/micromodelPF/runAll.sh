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

echo "=== Run: diffInterFoam/micromodelPF ==="
run_task createMesh.sh
run_task runSnappyHexMesh.sh
run_task initCaseSPFlow.sh
run_task runCaseSPFlow.sh
run_task processSPFlow.sh
run_task initCaseTPFlow.sh
run_task runCaseTPFlow.sh
run_task processTPFlow.sh
echo "=== Completed: diffInterFoam/micromodelPF ==="

