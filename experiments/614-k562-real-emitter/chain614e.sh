#!/bin/bash
# 614 A1: the sub-gate reservoir -> reference -> G3 key, then the gate runs G1 / G2 / G3 (seeds 0-2) and the verdict,
# then the H1@11k control. One gate602 ticket per run.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/614-k562-real-emitter; L=$X/work/logs/chain614e.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
MDE="--mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
log "chain614e start"
bash experiments/602-council-gate/gate602.sh 614-sub 6 -- /usr/bin/time -l uv run python $X/sub614.py > $X/work/logs/sub614.log 2>&1; log "sub614 rc=$?"
bash experiments/602-council-gate/gate602.sh 614-ref-k562s 8 -- /usr/bin/time -l uv run python $X/ref614.py $X/work/k562-gwps-c250-sub.h5 $X/work/targets614.txt $X/work/ref_k562s > $X/work/logs/ref_k562s.log 2>&1; log "ref_k562s rc=$?"
uv run python $X/key614cp.py > $X/work/logs/key614cp.log 2>&1; rc=$?; log "key614cp rc=$rc"
[ $rc -eq 0 ] || exit 1
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
