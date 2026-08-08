#!/usr/bin/env bash
set -uo pipefail
ROOT="/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/Blast.jl"
LOG="$ROOT/benchmark/reactant_full_reverse_long.log"
STATUS="$ROOT/benchmark/reactant_full_reverse_long.status"
START=$(date -u +%s)
{
  echo "start_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "branch=$(git -C "$ROOT" branch --show-current)"
  echo "head=$(git -C "$ROOT" rev-parse HEAD)"
  echo "timeout=6h"
  echo "driver=benchmark/run_staged_reactant_gradient.jl"
} > "$STATUS"
/usr/bin/time -v timeout --preserve-status --signal=TERM --kill-after=10m 6h \
  julia -t 8 --project="$ROOT/benchmark" "$ROOT/benchmark/run_staged_reactant_gradient.jl" \
  > "$LOG" 2>&1
CODE=$?
{
  echo "end_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "elapsed_seconds=$(($(date -u +%s)-START))"
  echo "exit_code=$CODE"
} >> "$STATUS"
exit "$CODE"
