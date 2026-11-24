#!/bin/bash
# 612 (10-02 ~15:30Z): idall3d seeds 1-2 (--dump 85) for [d4cd5f]'s 3-seed 630 stack on the shipping hook, AFTER chain612f
# exits, each through gate602 (v3: one ticket per run). KA per seed: stats == idall3's scored seed. Blocks merged file by
# file into h1_idall3d/blocks (no overwrite). Terminal line "<ts> END" on every path.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/612-key-gated-hook; L=$X/work/logs/chain612g.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
ARGS="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
log "chain612g start; waiting for chain612f pid 91776"
while ps -p 91776 >/dev/null 2>&1; do sleep 20; done
for s in 1 2; do
  T=stats-base_s${s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz
  log "run idall3d s$s --dump 85"
  bash experiments/602-council-gate/gate602.sh "612-idall3d-s$s" 8 -- uv run python $X/run612.py idall3d base $s $ARGS --dump 85 \
    > $X/work/logs/idall3d-s$s.log 2>&1
  log "done idall3d s$s rc=$?"
  if cmp -s "$X/work/h1_idall3d/$T" "$X/work/h1_idall3/$T"; then log "KA idall3d s$s stats == idall3 s$s: PASS"; else log "KA idall3d s$s stats == idall3 s$s: FAIL"; fi
  n=0
  mkdir -p "$X/work/h1_idall3d/blocks"
  for f in "$X/work/h1_idall3d/work/blocks/"*_s${s}_*.npz; do
    [ -e "$f" ] || continue
    b=$(basename "$f"); [ -e "$X/work/h1_idall3d/blocks/$b" ] && continue
    mv "$f" "$X/work/h1_idall3d/blocks/$b" && n=$((n+1))
  done
  log "blocks s$s -> $X/work/h1_idall3d/blocks ($n files)"
done
