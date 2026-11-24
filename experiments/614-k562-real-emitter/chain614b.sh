#!/bin/bash
# 614 stage 5-6: the wrapper KA on H1 (vs 435), then the K562 gate runs G1 (zs) and G2 (MDE), seeds 0-2.
# One gate602 ticket per run. Terminal "<ts> END" on every path.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/614-k562-real-emitter; L=$X/work/logs/chain614b.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
MDE="--mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
run(){ local nm=$1 mode=$2 tag=$3; shift 3
  log "run $tag"
  bash experiments/602-council-gate/gate602.sh "614-$tag" 8 -- uv run python $X/run614.py $mode $nm base "$@" > $X/work/logs/$tag.log 2>&1
  log "done $tag rc=$?"; }
log "chain614b start"
run kah1 h1ka kah1-s0 0 --bulk ps $MDE --zs dn_k250_zs
uv run python $X/ka614.py >> $L 2>&1; KA=$?
log "ka614 rc=$KA"
if [ $KA -ne 0 ]; then log "wrapper KA FAILED: no K562 runs"; exit 1; fi
for s in 0 1 2; do
  run g1 k562 g1-s$s $s --bulk base $MDE --zero-repair -5 --zs dn_k250_zs
  run g2nomde k562 g2nomde-s$s $s --bulk base
  run g2mde k562 g2mde-s$s $s --bulk base $MDE
done
uv run python $X/gate614.py >> $L 2>&1
log "gate614 rc=$?"
