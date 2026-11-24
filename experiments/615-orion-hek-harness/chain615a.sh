#!/bin/bash
# 615 stage 1-2 (after the blob PASS, 395c186f): truth (P'), reference, keys, hook signs, then the low-power diagnostic.
# One gate602 ticket per heavy run. A failing step ends the chain with a terminal END line (AGENTS 6.2).
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/615-orion-hek-harness; L=$X/work/logs/chain615a.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
log "chain615a start"
bash experiments/602-council-gate/gate602.sh 615-truth 6 -- /usr/bin/time -l uv run python $X/truth615.py > $X/work/logs/truth615.log 2>&1
rc=$?; log "truth615 rc=$rc"; [ $rc -eq 0 ] || exit 1
bash experiments/602-council-gate/gate602.sh 615-ref 10 -- /usr/bin/time -l uv run python experiments/614-k562-real-emitter/ref614.py $X/work/hek-truth.h5 $X/work/targets615.txt $X/work/ref_hek > $X/work/logs/ref_hek.log 2>&1
rc=$?; log "ref_hek rc=$rc"; [ $rc -eq 0 ] || exit 1
uv run python $X/diag615.py >> $L 2>&1; log "diag615 rc=$?"
uv run python $X/key615.py > $X/work/logs/key615.log 2>&1; rc=$?; log "key615 rc=$rc"; [ $rc -eq 0 ] || exit 1
bash experiments/602-council-gate/gate602.sh 615-keycp 6 -- /usr/bin/time -l uv run python $X/key615cp.py > $X/work/logs/key615cp.log 2>&1
rc=$?; log "key615cp rc=$rc"; [ $rc -eq 0 ] || exit 1
uv run python $X/idsign615.py > $X/work/logs/idsign615.log 2>&1; log "idsign615 rc=$?"
