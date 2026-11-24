#!/bin/bash
# 614 hook forms (PRED-614-levers.md): after chain614e (the gate) exits, c4h and idall3 on G1's construction, seeds 0-2.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/614-k562-real-emitter; L=$X/work/logs/chain614h.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
MDE="--mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
P=$(cat $X/work/chain614e2.pid); log "chain614h start; waiting for chain614e pid $P"
while ps -p "$P" >/dev/null 2>&1; do sleep 30; done
for s in 0 1 2; do
  for pair in "k562c4h hk_c4h" "k562idall3 hk_id"; do
    set -- $pair; mode=$1; nm=$2
    log "run $nm-s$s"
    bash experiments/602-council-gate/gate602.sh "614-$nm-s$s" 8 -- uv run python $X/run614.py $mode $nm base $s --bulk base $MDE --zero-repair -5 --zs dn_k250_zs > $X/work/logs/$nm-s$s.log 2>&1
    log "done $nm-s$s rc=$?"
  done
done
uv run python $X/hook614.py >> $L 2>&1; log "hook614 rc=$?"
