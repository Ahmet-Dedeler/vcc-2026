#!/bin/bash
# 614 descriptive control H1@11k: after chain614e exits -> thin h1-test -> reference -> G1 / G2 runs (seeds 0-2) -> summary.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/614-k562-real-emitter; L=$X/work/logs/chain614d.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
MDE="--mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
P=$(cat $X/work/chain614e2.pid); log "chain614d start; waiting for chain614e pid $P"
while ps -p "$P" >/dev/null 2>&1; do sleep 30; done
bash experiments/602-council-gate/gate602.sh 614-thin 8 -- /usr/bin/time -l uv run python $X/thin614.py > $X/work/logs/thin614.log 2>&1
log "thin614 rc=$?"
bash experiments/602-council-gate/gate602.sh 614-ref-h1t 6 -- uv run python $X/ref614.py $X/work/h1-thin11k.h5 $X/work/targets_h1.txt $X/work/ref_h1t11k > $X/work/logs/ref_h1t11k.log 2>&1
log "ref_h1t11k rc=$?"
run(){ local nm=$1 tag=$2; shift 2
  log "run $tag"
  bash experiments/602-council-gate/gate602.sh "614-$tag" 8 -- uv run python $X/run614.py h1t11k $nm base "$@" > $X/work/logs/$tag.log 2>&1
  log "done $tag rc=$?"; }
for s in 0 1 2; do
  run g1t g1t-s$s $s --bulk ps $MDE --zs dn_k250_zs
  run g2t_nomde g2t_nomde-s$s $s --bulk ps
  run g2t_mde g2t_mde-s$s $s --bulk ps $MDE
done
uv run python $X/h1t614.py >> $L 2>&1
log "h1t614 rc=$?"
