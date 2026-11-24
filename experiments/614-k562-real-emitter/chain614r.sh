#!/bin/bash
# 614 A1-R: the reservoir x0.5 / x2 truths and references (for capacity-bound levers before they are read).
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/614-k562-real-emitter; L=$X/work/logs/chain614r.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
for sc in 0.5 2; do
  SUB614_SCALE=$sc bash experiments/602-council-gate/gate602.sh 614-sub$sc 6 -- uv run python $X/sub614.py > $X/work/logs/sub614-x$sc.log 2>&1; log "sub614 x$sc rc=$?"
  bash experiments/602-council-gate/gate602.sh 614-ref-s$sc 8 -- uv run python $X/ref614.py $X/work/k562-gwps-c250-sub$sc.h5 $X/work/targets614.txt $X/work/ref_k562s$sc > $X/work/logs/ref_k562s$sc.log 2>&1; log "ref_k562s$sc rc=$?"
done
