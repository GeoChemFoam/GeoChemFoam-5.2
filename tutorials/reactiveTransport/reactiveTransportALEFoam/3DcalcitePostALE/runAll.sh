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

echo "=== Run: 3DcalcitePostALE ==="
run_task createMesh.sh
run_task runSnappyHexMesh.sh
run_task initCase.sh
run_task runCase0.sh
run_task runCase.sh
run_task processCase.sh
echo "=== Completed: 3DcalcitePostALE ==="

