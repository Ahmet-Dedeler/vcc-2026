#!/bin/bash
# 618 stage C (TRUSTGATE.md, 50f6708c): after chain618a (faithful key, ka618 PASS) and chain618b (truth / reference /
# power diagnostic) have EXITED with their success lines: the wrapper KA (run618e h1ka == 435), the gate runs G1 / G2 / G3
# x seeds 0-2 (--bulk ps, H1-like as 611), gate618.py, the zero-total counts; then, only on a delivery-only PASS, the hook
# forms (c4h / idall3 / idkg3 x 3 seeds) for the 613d re-check. One gate602 ticket per run; every exit path writes END.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/618-h1-trainval-surface; L=$X/work/logs/chain618c.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
MDE="--mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
for c in chain618a chain618b; do
  P=$(cat $X/work/$c.pid); log "waiting for $c pid $P"
  while ps -p "$P" >/dev/null 2>&1; do sleep 30; done
done
/usr/bin/grep -q "stage A done" $X/work/logs/chain618a.log || { log "chain618a did not finish (no 'stage A done'): stopping"; exit 1; }
/usr/bin/grep -q "ka618 rc=0" $X/work/logs/chain618a.log || { log "ka618 did not pass: stopping"; exit 1; }
/usr/bin/grep -q "ref_tv rc=0" $X/work/logs/chain618b.log || { log "ref_tv missing: stopping"; exit 1; }
run(){ local mode=$1 nm=$2 tag=$3; shift 3
  log "run $tag"
  bash experiments/602-council-gate/gate602.sh "618-$tag" 8 -- uv run python $X/run618e.py $mode $nm base "$@" > $X/work/logs/$tag.log 2>&1
  log "done $tag rc=$?"; }
run h1ka kah1 kah1-s0 0 --bulk ps $MDE --zs dn_k250_zs
uv run python $X/ka618e.py >> $L 2>&1; rc=$?; log "ka618e rc=$rc"; [ $rc -eq 0 ] || exit 1
for s in 0 1 2; do
  run tv g1 g1-s$s $s --bulk ps $MDE --zero-repair -5 --zs dn_k250_zs
  run tv g2nomde g2nomde-s$s $s --bulk ps
  run tv g2mde g2mde-s$s $s --bulk ps $MDE
  run tvcp g3cp g3cp-s$s $s --bulk ps $MDE --zero-repair -5 --zs dn_k250_zs
done
uv run python $X/gate618.py >> $L 2>&1; log "gate618 rc=$?"
uv run python experiments/654-carrier-amplitude/zerocount654.py $X/work/g1/stats-*.npz $X/work/g2nomde/stats-*.npz $X/work/g2mde/stats-*.npz $X/work/g3cp/stats-*.npz >> $L 2>&1
log "zerocount rc=$?"
/usr/bin/grep -q '"delivery_only": "PASS"' $X/work/gate618.json || { log "delivery-only verdict not PASS: hook runs skipped (TRUSTGATE)"; exit 0; }
for s in 0 1 2; do
  run tvc4h hk_c4h hk_c4h-s$s $s --bulk ps $MDE --zero-repair -5 --zs dn_k250_zs
  run tvidall3 hk_id hk_id-s$s $s --bulk ps $MDE --zero-repair -5 --zs dn_k250_zs
  run tvidkg3 hk_kg hk_kg-s$s $s --bulk ps $MDE --zero-repair -5 --zs dn_k250_zs
done
log "stage C done"
