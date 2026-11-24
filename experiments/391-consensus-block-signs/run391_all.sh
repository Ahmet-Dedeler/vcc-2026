#!/usr/bin/env bash
# 391 driver: KA, three arms, report -- serial, memory-gated, stops on the first failure.
set -uo pipefail
cd "$(dirname "$0")/../.."
X=experiments/391-consensus-block-signs; LOG=$X/work/driver391.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
export NUMBA_CACHE_DIR=$X/work/numba-cache NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
for arm in cap15 k3c5b k3c5bshuf k3down; do
  n=0; until [ "$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')" -ge 35 ]; do n=$((n+1)); [ $n -ge 120 ] && { log "GATE never opened before $arm"; exit 3; }; sleep 30; done
  log "START $arm"
  if ! uv run python $X/run391.py $arm >> $X/work/arm-$arm.log 2>&1; then log "FAILED $arm (see work/arm-$arm.log): $(grep -E 'Error|assert|Traceback' $X/work/arm-$arm.log | tail -2 | tr '\n' ' ')"; exit 1; fi
  log "DONE $arm $(grep -E '^'$arm' ' $X/work/arm-$arm.log | tail -1 | cut -c1-300)"
done
if uv run python $X/report391.py k3c5b k3c5bshuf k3down >> $X/work/report391.log 2>&1; then log "REPORT OK $(tr '\n' ' ' < $X/work/report391.log | cut -c1-900)"; else log "REPORT FAILED $(tail -3 $X/work/report391.log | tr '\n' ' ')"; exit 1; fi
log "ALL DONE"
