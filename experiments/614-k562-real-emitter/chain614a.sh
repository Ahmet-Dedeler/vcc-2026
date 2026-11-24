#!/bin/bash
# 614 stage 3: the K562 reference, then ref614's KA on H1 (must equal 323's reference). One gate602 ticket per run.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/614-k562-real-emitter; L=$X/work/logs/chain614a.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
log "chain614a start"
bash experiments/602-council-gate/gate602.sh 614-ref-k562 6 -- /usr/bin/time -l uv run python $X/ref614.py $X/work/k562-gwps-c250.h5 $X/work/targets614.txt $X/work/ref_k562 > $X/work/logs/ref_k562.log 2>&1
log "ref_k562 rc=$?"
bash experiments/602-council-gate/gate602.sh 614-ref-h1ka 6 -- /usr/bin/time -l uv run python $X/ref614.py data/harness/h1-test.h5ad $X/work/targets_h1.txt $X/work/ref_h1ka > $X/work/logs/ref_h1ka.log 2>&1
log "ref_h1ka rc=$?"
uv run python $X/kaRef614.py >> $L 2>&1
log "kaRef614 rc=$?"
