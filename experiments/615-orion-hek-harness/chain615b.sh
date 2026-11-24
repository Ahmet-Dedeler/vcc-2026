#!/bin/bash
# 615 gate (TRUSTGATE.md, 6e35ff09): after chain615a's PROCESS exits and its stages all logged rc=0, the wrapper KA
# (run615 h1ka == 435), then G1 / G2 / G3 x seeds 0-2 on the HEK harness, then gate615.py and the zero-total counts.
# One gate602 ticket per run; every exit path writes a terminal END line.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/615-orion-hek-harness; L=$X/work/logs/chain615b.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
MDE="--mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
P=$(cat $X/work/chain615a.pid); log "chain615b start; waiting for chain615a pid $P"
while ps -p "$P" >/dev/null 2>&1; do sleep 30; done
for step in truth615 ref_hek key615 key615cp idsign615; do
  /usr/bin/grep -q "$step rc=0" $X/work/logs/chain615a.log || { log "stage $step did not log rc=0: stopping"; exit 1; }
done
run(){ local mode=$1 nm=$2 tag=$3; shift 3
  log "run $tag"
  bash experiments/602-council-gate/gate602.sh "615-$tag" 8 -- uv run python $X/run615.py $mode $nm base "$@" > $X/work/logs/$tag.log 2>&1
  log "done $tag rc=$?"; }
run h1ka kah1 kah1-s0 0 --bulk ps $MDE --zs dn_k250_zs
uv run python $X/ka615.py >> $L 2>&1; rc=$?; log "ka615 rc=$rc"; [ $rc -eq 0 ] || exit 1
for s in 0 1 2; do
  run hek g1 g1-s$s $s --bulk base $MDE --zero-repair -5 --zs dn_k250_zs
  run hek g2nomde g2nomde-s$s $s --bulk base
  run hek g2mde g2mde-s$s $s --bulk base $MDE
  run hekcp g3cp g3cp-s$s $s --bulk base $MDE --zero-repair -5 --zs dn_k250_zs
done
uv run python $X/gate615.py >> $L 2>&1; log "gate615 rc=$?"
uv run python experiments/654-carrier-amplitude/zerocount654.py $X/work/g1/stats-*.npz $X/work/g2nomde/stats-*.npz $X/work/g2mde/stats-*.npz $X/work/g3cp/stats-*.npz >> $L 2>&1
log "zerocount rc=$?"
