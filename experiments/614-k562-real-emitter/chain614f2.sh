#!/bin/bash
# 614 faithful-key re-run, RESUMED (10-03 ~15:00Z) after fkey stage died on donor_weights() (difference (8), TRUSTGATE
# amendment): re-freeze -> commit -> stage -> operator (gate602 v6 quiet admission; the same code's measured peak 13.7 GB)
# -> compose -> the same 12 gate runs (run614 k562f / k562fcp) -> gate614f.py. END on every exit path.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/614-k562-real-emitter; L=$X/work/logs/chain614f.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
MDE="--mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
log "chain614f2 start (resume at freeze, after the donor-weights fix 46339de9)"
step(){ local s=$1 gb=$2 how=$3
  case $how in
    none) /usr/bin/time -l uv run python $X/fkey614.py $s > $X/work/logs/fkey-$s.log 2>&1 ;;
    gate) bash experiments/602-council-gate/gate602.sh "614f-$s" $gb -- /usr/bin/time -l uv run python $X/fkey614.py $s > $X/work/logs/fkey-$s.log 2>&1 ;;
    quiet) GATE602_QUIET=1 bash experiments/602-council-gate/gate602.sh "614f-$s" $gb -- /usr/bin/time -l uv run python $X/fkey614.py $s > $X/work/logs/fkey-$s.log 2>&1 ;;
  esac
  local rc=$?; log "fkey $s rc=$rc"; [ $rc -eq 0 ] || exit 1; }
step freeze 0 none
git add -f $X/work/fkey/inputs.json && git commit -q -m "614 fkey: inputs re-frozen after the donor-weights fix, before any measurement" -- $X/work/fkey/inputs.json
log "fkey inputs committed rc=$? $(git log -1 --format=%h -- $X/work/fkey/inputs.json)"
step stage 8 gate
step operator 14 quiet
step compose 6 gate
run(){ local mode=$1 nm=$2 tag=$3; shift 3
  log "run $tag"
  bash experiments/602-council-gate/gate602.sh "614-$tag" 8 -- uv run python $X/run614.py $mode $nm base "$@" > $X/work/logs/$tag.log 2>&1
  log "done $tag rc=$?"; }
for s in 0 1 2; do
  run k562f g1f g1f-s$s $s --bulk base $MDE --zero-repair -5 --zs dn_k250_zs
  run k562f g2nomdef g2nomdef-s$s $s --bulk base
  run k562f g2mdef g2mdef-s$s $s --bulk base $MDE
  run k562fcp g3cpf g3cpf-s$s $s --bulk base $MDE --zero-repair -5 --zs dn_k250_zs
done
uv run python $X/gate614f.py >> $L 2>&1; log "gate614f rc=$?"
