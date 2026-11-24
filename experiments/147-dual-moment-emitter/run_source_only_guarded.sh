#!/bin/bash
set -uo pipefail

ROOT="~/Code/virtual-cell-challenge-2026"
WORK="$ROOT/experiments/147-dual-moment-emitter/work/source-feasibility-v2"
READY="$WORK/memory-ready.json"
STOP="$WORK/memory-stop"
REPORT="$WORK/memory-guard.json"

if [[ -e "$WORK" ]]; then
  echo "immutable source-feasibility-v2 directory already exists: $WORK" >&2
  exit 2
fi

mkdir -p "$WORK"
cd "$ROOT" || exit 2

PYTHONDONTWRITEBYTECODE=1 uv run python \
  experiments/147-dual-moment-emitter/run_feasibility.py &
ROOT_PID=$!

uv run python scripts/run_with_memory_guard.py \
  --monitor-pid "$ROOT_PID" \
  --ready-file "$READY" \
  --stop-file "$STOP" \
  --report-file "$REPORT" \
  --max-rss-gib 6 \
  --min-free-percent 20 \
  --poll-seconds 0.5 \
  --quiet &
MONITOR_PID=$!

for _ in $(seq 1 100); do
  [[ -f "$READY" ]] && break
  kill -0 "$ROOT_PID" 2>/dev/null || break
  sleep 0.05
done

wait "$ROOT_PID"
STATUS=$?
touch "$STOP"
wait "$MONITOR_PID"
MONITOR_STATUS=$?

if [[ "$STATUS" -ne 0 ]]; then
  exit "$STATUS"
fi
exit "$MONITOR_STATUS"
