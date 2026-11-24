#!/bin/bash
# 612 (reordered 12:15Z): the decisive idkg seeds first, then c4kg s2; each through gate602 (FIFO). c4kg s1 runs separately
# (gate pid 61796, launched by chain612). Terminal line "<ts> END" on every path.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/612-key-gated-hook; L=$X/work/logs/chain612b.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
log "chain612b start"
for job in "idkg 0" "idkg 1" "idkg 2" "c4kg 2"; do
  set -- $job; mode=$1; s=$2
  log "run $mode s$s"
  bash experiments/602-council-gate/gate602.sh "612-$mode-s$s" 8 -- uv run python $X/run612.py $mode base $s $ARGS \
    > $X/work/logs/$mode-s$s.log 2>&1
  log "done $mode s$s rc=$?"
done
log "END"
