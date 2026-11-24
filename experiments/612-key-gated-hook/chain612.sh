#!/bin/bash
# 612: the H1 real-emitter runs, serial, each through the council gate (one of ours at a time). Waits for the already-queued
# KA (c4kgall s0, gate pid passed as $1) to finish first. Terminal line "<ts> END" on every path.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/612-key-gated-hook; L=$X/work/logs/chain612.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
KA_PID=${1:-}
ARGS="base %s --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
log "chain612 start (KA gate pid ${KA_PID:-none})"
if [ -n "$KA_PID" ]; then while kill -0 "$KA_PID" 2>/dev/null; do sleep 30; done; fi
for mode in c4kg idkg; do
  for s in 0 1 2; do
    log "run $mode s$s"
    bash experiments/602-council-gate/gate602.sh "612-$mode-s$s" 8 -- uv run python $X/run612.py $mode $(printf "$ARGS" "$s") \
      > $X/work/logs/$mode-s$s.log 2>&1
    log "done $mode s$s rc=$? ($(tail -1 $X/work/logs/$mode-s$s.log | cut -c1-120))"
  done
done
log "END"
