#!/bin/bash
# 618 stage D (after chain618c EXITS): if the gate's delivery-only verdict is PASS, the registered 613d re-check
# (recheck618.py), then the two magnitude stacks (PRED-618-stacks.md): plain / x1.25 dump runs (3 seeds, zero repair -5),
# the KA stk_base == g1, rnup618 (a) / (b) / 630-alone, price618, zero-total counts, stacks618.py's decision.
# work/STACKS-DONE is written on EVERY exit path, so 614's faithful-key re-run (which waits for it) never blocks forever.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/618-h1-trainval-surface; L=$X/work/logs/chain618d.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"; touch $X/work/STACKS-DONE' EXIT
MDE="--mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
P=$(cat $X/work/chain618c.pid); log "chain618d start; waiting for chain618c pid $P"
while ps -p "$P" >/dev/null 2>&1; do sleep 60; done
/usr/bin/grep -q '"delivery_only": "PASS"' $X/work/gate618.json 2>/dev/null || { log "delivery-only verdict not PASS (or no gate618.json): re-check and stacks SKIPPED"; exit 0; }
uv run python $X/recheck618.py >> $L 2>&1; log "recheck618 rc=$?"
run(){ local mode=$1 nm=$2 design=$3 s=$4; shift 4
  log "run $nm-s$s"
  bash experiments/602-council-gate/gate602.sh "618-$nm-s$s" 8 -- uv run python $X/run618e.py $mode $nm $design $s "$@" > $X/work/logs/$nm-s$s.log 2>&1
  local rc=$?; log "done $nm-s$s rc=$rc"
  if [ -d "$X/work/$nm/work/blocks" ]; then
    mkdir -p "$X/work/$nm/blocks" && mv "$X/work/$nm/work/blocks/"*.npz "$X/work/$nm/blocks/" && log "blocks -> $X/work/$nm/blocks ($(ls $X/work/$nm/blocks | wc -l | tr -d ' ') files)"
  fi
  return $rc; }
for s in 0 1 2; do
  run tv stk_base base $s --bulk ps $MDE --zs dn_k250_zs --zero-repair -5 --dump 126 || exit 1
  run tv stk_x125 x125 $s --bulk ps $MDE --zs dn_k250_zs --zero-repair -5 --dump 126 || exit 1
done
uv run python $X/ka618stk.py >> $L 2>&1; rc=$?; log "ka618stk rc=$rc"; [ $rc -eq 0 ] || exit 1
TB="base_s{seed}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse"; TX="x125_s{seed}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse"
rn(){ local s=$1 v=$2 src=$3 nm=$4 tag=$5
  bash experiments/602-council-gate/gate602.sh "618-rnup-$nm-s$s" 6 -- uv run python $X/rnup618.py $s $v --src $X/work/$src --name $nm --tag "$tag" > $X/work/logs/rnup-$nm-s$s.log 2>&1
  log "rnup $nm s$s rc=$?"; }
for s in 0 1 2; do
  rn $s b_k100_q20r3p10 stk_x125 t618a "$TX"
  rn $s b_k100_q20r3p10_D200 stk_base t618b "$TB"
  rn $s b_k100_q20r3p10 stk_base t618u "$TB"
  for f in $X/work/stk_x125/stats-x125_s${s}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz; do
    [ -f "$f" ] && cp -c "$f" "$X/work/stacks/stats-rnup_t618x_s${s}_x125.npz"
  done
done
for spec in "t618a b_k100_q20r3p10" "t618b b_k100_q20r3p10_D200" "t618u b_k100_q20r3p10" "t618x x125"; do
  set -- $spec
  bash experiments/602-council-gate/gate602.sh "618-price-$1" 6 -- uv run python $X/price618.py $2 --base $X/work/stk_base --name $1 --tag "$TB" >> $X/work/logs/price618-$1.log 2>&1
  log "price618 $1 rc=$?"
done
uv run python experiments/654-carrier-amplitude/zerocount654.py $X/work/stk_base/stats-*.npz $X/work/stk_x125/stats-*.npz $X/work/stacks/stats-rnup_*.npz >> $L 2>&1
log "zerocount rc=$?"
uv run python $X/stacks618.py >> $L 2>&1; log "stacks618 rc=$?"
log "stage D done"
