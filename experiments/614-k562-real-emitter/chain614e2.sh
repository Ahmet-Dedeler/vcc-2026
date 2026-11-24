#!/bin/bash
# 614 A1 (relaunch after the NaN-safe metrics fix): the gate runs G1 / G2 / G3 (seeds 0-2) and the verdict.
# then the H1@11k control. One gate602 ticket per run.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/614-k562-real-emitter; L=$X/work/logs/chain614e2.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
MDE="--mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
log "chain614e2 start (data prep already done in chain614e)"
run(){ local mode=$1 nm=$2 tag=$3; shift 3
  log "run $tag"
  bash experiments/602-council-gate/gate602.sh "614-$tag" 8 -- uv run python $X/run614.py $mode $nm base "$@" > $X/work/logs/$tag.log 2>&1
  log "done $tag rc=$?"; }
for s in 0 1 2; do
  run k562 g1 g1-s$s $s --bulk base $MDE --zero-repair -5 --zs dn_k250_zs
  run k562 g2nomde g2nomde-s$s $s --bulk base
  run k562 g2mde g2mde-s$s $s --bulk base $MDE
  run k562cp g3cp g3cp-s$s $s --bulk base $MDE --zero-repair -5 --zs dn_k250_zs
done
uv run python $X/gate614.py >> $L 2>&1; log "gate614 rc=$?"
