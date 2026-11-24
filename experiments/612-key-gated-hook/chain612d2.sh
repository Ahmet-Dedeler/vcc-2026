#!/bin/bash
# 612 (10-02 ~14:25Z): chain612d2 = chain612d WITHOUT the wait on chain612c (gate602 v3 has a 3rd slot; 4d0207),
# each through gate602, then score612. Pre-registered: 612/PREDICTION-idall3.md. Terminal line "<ts> END" on every path.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/612-key-gated-hook; L=$X/work/logs/chain612d2.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
log "chain612d2 start (no wait)"
for s in 0 1 2; do
  log "run idall3 s$s"
  bash experiments/602-council-gate/gate602.sh "612-idall3-s$s" 8 -- uv run python $X/run612.py idall3 base $s $ARGS \
    > $X/work/logs/idall3-s$s.log 2>&1
  log "done idall3 s$s rc=$?"
done
uv run python $X/score612.py > $X/work/logs/score612-idall3.log 2>&1
log "score612 rc=$?"
