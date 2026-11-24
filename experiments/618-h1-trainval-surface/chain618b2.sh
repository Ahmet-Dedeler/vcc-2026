#!/bin/bash
# 618 stage B, RESUMED at the reference (19:52Z; truth618 already rc=0): the 12 GB ticket enters through polite618.sh so it
# never heads gate602's queue. Logs to the same chain618b.log. Every exit path writes END.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/618-h1-trainval-surface; L=$X/work/logs/chain618b.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
log "chain618b2 start (resume at ref_tv)"
bash $X/polite618.sh 618-ref 12 -- /usr/bin/time -l uv run python experiments/614-k562-real-emitter/ref614.py $X/work/tv-truth.h5 $X/work/targets_tv.txt $X/work/ref_tv > $X/work/logs/ref_tv.log 2>&1
rc=$?; log "ref_tv rc=$rc"; [ $rc -eq 0 ] || exit 1
uv run python $X/diag618.py >> $L 2>&1; log "diag618 rc=$?"
