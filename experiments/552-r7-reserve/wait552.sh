#!/bin/bash
# 552: build R7b's surface (build_r7b.py, construction only) after drive551's efp KA process exits and gate551 opens,
# so the two jobs never overlap. Capped at 2026-10-03 00:00Z; logs to work/build_r7b.log.
cd "$(dirname "$0")/../.." || exit 1
LOG=experiments/552-r7-reserve/work/build_r7b.log
KA_PID="$1"
echo "$(date -u +%FT%TZ) armed (waits for pid $KA_PID, then gate551)" >> "$LOG"
while kill -0 "$KA_PID" 2>/dev/null; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610030000 ] && { echo "$(date -u +%FT%TZ) cap reached waiting for the KA; END" >> "$LOG"; exit 1; }
  sleep 60; done
until uv run python experiments/551-r6-c10b/gate551.py > /dev/null 2>&1; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610030000 ] && { echo "$(date -u +%FT%TZ) cap reached waiting for the gate; END" >> "$LOG"; exit 1; }
  sleep 60; done
echo "$(date -u +%FT%TZ) gate open; building" >> "$LOG"
/usr/bin/time -l uv run python experiments/552-r7-reserve/build_r7b.py >> "$LOG" 2>&1
echo "$(date -u +%FT%TZ) rc=$? END" >> "$LOG"
