#!/bin/bash
# 614 G3: the call-propensity block (lcap, kept) on the K562 harness, seeds 0-2, args as G1 (the zs base); then the gate.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/614-k562-real-emitter; L=$X/work/logs/chain614c.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
MDE="--mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
log "chain614c start; waiting for the wrapper KA in chain614b"
until grep -q "ka614 rc=" $X/work/logs/chain614b.log 2>/dev/null; do sleep 20; done
if ! grep -q "ka614 rc=0" $X/work/logs/chain614b.log; then log "wrapper KA failed: no G3 runs"; exit 1; fi
for s in 0 1 2; do
  log "run g3cp-s$s"
  bash experiments/602-council-gate/gate602.sh "614-g3cp-s$s" 8 -- uv run python $X/run614.py k562cp g3cp base $s --bulk base $MDE --zero-repair -5 --zs dn_k250_zs > $X/work/logs/g3cp-s$s.log 2>&1
  log "done g3cp-s$s rc=$?"
done
P=$(cat $X/work/chain614b.pid); while ps -p "$P" >/dev/null 2>&1; do sleep 20; done
uv run python $X/gate614.py >> $L 2>&1
log "gate614 rc=$?"
