#!/bin/bash
# 612 (10-02 ~14:12Z): block dumps for [d4cd5f]'s 630 stack on the hooks that ship: c4kgd and idkg3d, seed 0, --dump 85,
# each through gate602. KA: each dump run's stats == the scored mode's seed-0 stats (cmp). The blocks are then moved from
# h1_<mode>/work/blocks to h1_<mode>/blocks, the path rnup630 --src reads. Terminal line "<ts> END" on every path.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/612-key-gated-hook; L=$X/work/logs/chain612e.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
T=stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz
log "chain612e start"
for pair in "c4kgd c4kg" "idkg3d idkg3"; do
  set -- $pair; m=$1; ref=$2
  log "run $m s0 --dump 85"
  bash experiments/602-council-gate/gate602.sh "612-$m-s0" 8 -- uv run python $X/run612.py $m base 0 $ARGS --dump 85 \
    > $X/work/logs/$m-s0.log 2>&1
  rc=$?
  log "done $m s0 rc=$rc"
  if cmp -s "$X/work/h1_$m/$T" "$X/work/h1_$ref/$T"; then log "KA $m s0 stats == $ref s0: PASS"; else log "KA $m s0 stats == $ref s0: FAIL (or $ref s0 absent)"; fi
  if [ -d "$X/work/h1_$m/work/blocks" ] && [ ! -e "$X/work/h1_$m/blocks" ]; then
    mv "$X/work/h1_$m/work/blocks" "$X/work/h1_$m/blocks" && log "blocks -> $X/work/h1_$m/blocks ($(ls $X/work/h1_$m/blocks | wc -l | tr -d ' ') files)"
  else
    log "blocks not moved (source missing or target exists)"
  fi
done
