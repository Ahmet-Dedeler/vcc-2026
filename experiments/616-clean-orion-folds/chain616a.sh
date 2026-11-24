#!/bin/bash
# 616 (PLAN.md): after 615's stage 1 (chain615a) exits: S0 known answers, the clean HEK truth + gates, HCT116's raw
# extraction (gate602) + clean truth + gates, then the clean-mode instrument chain. Every exit path writes END.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/616-clean-orion-folds; L=$X/work/logs/chain616a.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
P=$(cat experiments/615-orion-hek-harness/work/chain615a.pid); log "chain616a start; waiting for chain615a pid $P"
while ps -p "$P" >/dev/null 2>&1; do sleep 30; done
SB="$X/sandbox616.py"
chain(){ local tag=$1 mode=$2
  B=experiments/613-block-prior-validation
  for s in experiments/639-identity-axis-lolo/decomp639.py $B/block613.py $B/block613b.py $B/block613c.py \
           $B/block613d.py $B/block613d3.py $B/block613d4.py $B/block613d5.py $B/block613d6.py $B/block613e.py \
           $B/block613e2.py experiments/612-key-gated-hook/ds612b.py experiments/612-key-gated-hook/ds612c.py \
           experiments/633-updose-folds/dose633.py; do
    n=$(basename $s .py)
    /usr/bin/time -l uv run python $SB $tag $mode $s > $X/work/logs/$tag-$n.log 2>&1
    log "$tag $n rc=$?"
  done; }
# S0 (a) gates480 on the thinned HEK harness through the sandbox; (b) the identity chain
/usr/bin/time -l uv run python $SB gka identity experiments/480-nmae-instrument/gates480.py orion_hek293t data/harness/orion-hek293t-thinned.h5ad --blk 1024 > $X/work/logs/gka.log 2>&1
log "gka rc=$?"
chain ka identity
uv run python $X/ka616.py gka ka >> $L 2>&1; rc=$?; log "ka616 rc=$rc"; [ $rc -eq 0 ] || exit 1
# S1: clean truths and gate tables
uv run python $X/truth616.py HEK293T > $X/work/logs/truth-HEK293T.log 2>&1; rc=$?; log "truth HEK293T rc=$rc"; [ $rc -eq 0 ] || exit 1
/usr/bin/time -l uv run python $SB gclean identity experiments/480-nmae-instrument/gates480.py orion_hek293t_clean $X/work/truth-HEK293T.h5 --blk 1024 > $X/work/logs/gclean-hek.log 2>&1
rc=$?; log "gates HEK293T clean rc=$rc"; [ $rc -eq 0 ] || exit 1
cp -c $X/work/sandbox/gclean/experiments/480-nmae-instrument/work/gate-orion_hek293t_clean.npz $X/work/gate-orion_hek293t_clean.npz
if [ ! -f experiments/615-orion-hek-harness/work/hct-raw.h5 ]; then
  bash experiments/602-council-gate/gate602.sh 616-extract-hct 9 -- /usr/bin/time -l uv run python experiments/615-orion-hek-harness/extract615.py --line HCT116 --out experiments/615-orion-hek-harness/work/hct-raw.h5 > $X/work/logs/extract-hct.log 2>&1
  rc=$?; log "extract HCT116 rc=$rc"; [ $rc -eq 0 ] || exit 1
fi
uv run python $X/truth616.py HCT116 > $X/work/logs/truth-HCT116.log 2>&1; rc=$?; log "truth HCT116 rc=$rc"; [ $rc -eq 0 ] || exit 1
/usr/bin/time -l uv run python $SB gclean identity experiments/480-nmae-instrument/gates480.py orion_hct116_clean $X/work/truth-HCT116.h5 --blk 1024 > $X/work/logs/gclean-hct.log 2>&1
rc=$?; log "gates HCT116 clean rc=$rc"; [ $rc -eq 0 ] || exit 1
cp -c $X/work/sandbox/gclean/experiments/480-nmae-instrument/work/gate-orion_hct116_clean.npz $X/work/gate-orion_hct116_clean.npz
# S2: the clean-mode chain
chain clean clean
log "S2 done"
