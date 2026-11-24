#!/bin/bash
# 661 (4) (PRED-661-625.md, committed before this runs): KA (flip618 h1 none == 612's c4kg, bit for bit), then 618's
# c4kg base and the 8 single-octile prior-only flips, seeds 0-2, zero repair -5; then score661oct.py.
# One gate602 ticket per run; END on every exit path.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/661-618-board-preds; L=$X/work/logs/chain661a.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
MDE="--mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
log "chain661a start"
run(){ local surf=$1 oct=$2 s=$3; shift 3
  log "run ${surf}_oct$oct-s$s"
  bash experiments/602-council-gate/gate602.sh "661-${surf}-oct$oct-s$s" 8 -- uv run python $X/flip618.py $surf $oct base $s --bulk ps $MDE --zs dn_k250_zs "$@" > $X/work/logs/${surf}_oct$oct-s$s.log 2>&1
  local rc=$?; log "done ${surf}_oct$oct-s$s rc=$rc"; return $rc; }
run h1 none 0 || exit 1
uv run python $X/ka661.py >> $L 2>&1; rc=$?; log "ka661 rc=$rc"; [ $rc -eq 0 ] || exit 1
for s in 0 1 2; do run tv none $s --zero-repair -5 || exit 1; done
for k in 0 1 2 3 4 5 6 7; do
  for s in 0 1 2; do run tv $k $s --zero-repair -5; done
done
bash experiments/602-council-gate/gate602.sh 661-score-oct 6 -- uv run python $X/score661oct.py >> $L 2>&1; log "score661oct rc=$?"
