#!/bin/bash
set -uo pipefail

ROOT="~/Code/virtual-cell-challenge-2026"
GUARD="$ROOT/experiments/147-dual-moment-emitter/work/h1-fast-build-guard-v2"
READY="$GUARD/memory-ready.json"
STOP="$GUARD/memory-stop"
REPORT="$GUARD/memory-guard.json"

if [[ -e "$GUARD" ]]; then
  echo "immutable build guard directory already exists: $GUARD" >&2
  exit 2
fi

mkdir -p "$GUARD"
cd "$ROOT" || exit 2

PYTHONDONTWRITEBYTECODE=1 uv run python \
  experiments/147-dual-moment-emitter/build_candidates.py &
ROOT_PID=$!

uv run python scripts/run_with_memory_guard.py \
  --monitor-pid "$ROOT_PID" \
  --ready-file "$READY" \
  --stop-file "$STOP" \
  --report-file "$REPORT" \
  --max-rss-gib 8 \
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
