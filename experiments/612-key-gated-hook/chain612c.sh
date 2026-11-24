#!/bin/bash
# 612 (10-02 ~14:00Z): KGID re-run in 604 v3's rule (idkg3: no identity-only extension), seeds 0-2, each through gate602
# (FIFO), then score612. d4cd5f flagged that the KGID prereg named the v1 artifact. Terminal line "<ts> END" on every path.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/612-key-gated-hook; L=$X/work/logs/chain612c.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
log "chain612c start"
for s in 0 1 2; do
  log "run idkg3 s$s"
  bash experiments/602-council-gate/gate602.sh "612-idkg3-s$s" 8 -- uv run python $X/run612.py idkg3 base $s $ARGS \
    > $X/work/logs/idkg3-s$s.log 2>&1
  log "done idkg3 s$s rc=$?"
done
uv run python $X/score612.py > $X/work/logs/score612-idkg3.log 2>&1
log "score612 rc=$?"
