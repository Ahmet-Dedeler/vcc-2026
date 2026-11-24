#!/bin/bash
# 618 stage B (PLAN.md): the h1-trainval truth (trainval cells + h1-test's controls), its reference (614's ref614.py),
# then the power diagnostic. One gate602 ticket per heavy step; every exit path writes END.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/618-h1-trainval-surface; L=$X/work/logs/chain618b.log
mkdir -p $X/work/logs
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
log "chain618b start"
bash experiments/602-council-gate/gate602.sh 618-truth 6 -- /usr/bin/time -l uv run python $X/truth618.py > $X/work/logs/truth618.log 2>&1
rc=$?; log "truth618 rc=$rc"; [ $rc -eq 0 ] || exit 1
bash experiments/602-council-gate/gate602.sh 618-ref 12 -- /usr/bin/time -l uv run python experiments/614-k562-real-emitter/ref614.py $X/work/tv-truth.h5 $X/work/targets_tv.txt $X/work/ref_tv > $X/work/logs/ref_tv.log 2>&1
rc=$?; log "ref_tv rc=$rc"; [ $rc -eq 0 ] || exit 1
uv run python $X/diag618.py >> $L 2>&1; log "diag618 rc=$?"
