#!/bin/bash
# 657 (PREREGISTRATION ed36f45f), Phase 1: emissions (c4h base with KA vs 618's hk_c4h; c4h_o4 = route A), the rnup copy's KA
# vs 618's t618u, and the 630-only stacks with / without _o4 (route B) on both emissions; seeds 0-2. One gate602 ticket per
# run. Phase 2 (DN variants) is a separate script. Terminal "<ts> END" on every path.
set -u
cd ~/Code/virtual-cell-challenge-2026 || exit 1
D=experiments/657-octile4-up-callneutral; LOG=$D/work/seq657.log; X618=experiments/618-h1-trainval-surface
MDE="--mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
TB="base_s{seed}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse"
ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
say() { echo "$(ts) $*" >> "$LOG"; }
emit() {  # mode name seed
  say "QUEUE emit $2 s$3"
  bash experiments/602-council-gate/gate602.sh "657-$2-s$3" 8 -- uv run python $D/run657e.py "$1" "$2" base "$3" --bulk ps $MDE \
    --zero-repair -5 --zs dn_k250_zs --dump 126 >> "$D/work/emit-$2-s$3.log" 2>&1
  local r=$?; say "DONE emit $2 s$3 rc=$r"
  if [ -d "$D/work/$2/work/blocks" ]; then
    mkdir -p "$D/work/$2/blocks" && mv "$D/work/$2/work/blocks/"*.npz "$D/work/$2/blocks/" && say "blocks -> $D/work/$2/blocks ($(ls $D/work/$2/blocks | wc -l | tr -d ' '))"
  fi
  return $r
}
rn() {  # seed src name variants...
  local s=$1 src=$2 nm=$3; shift 3
  say "QUEUE rnup $nm s$s: $*"
  bash experiments/602-council-gate/gate602.sh "657-rnup-$nm-s$s" 7 -- uv run python $D/rnup657.py "$s" "$@" --src "$src" --name "$nm" \
    --tag "$TB" >> "$D/work/rnup-$nm-s$s.log" 2>&1
  local r=$?; say "DONE rnup $nm s$s rc=$r"; return $r
}
say "START seq657 pid $$ (v2: KA-R = original rnup630 vs the copy on 657's own c4h s0 dumps; 618's stk_base dumps were cleaned)"
for s in 0 1 2; do
  emit tvc4h c4h "$s" || { say "emit c4h s$s failed: stop"; say "END"; exit 1; }
  uv run python $D/ka657.py E "$s" >> "$LOG" 2>&1 || { say "KA-E s$s FAILED: stop"; say "END"; exit 1; }
  rn "$s" "$D/work/c4h" c4h b_k100_q20r3p10 b_k100_q20r3p10_o4 || say "WARN rnup c4h s$s"
  if [ "$s" = 0 ]; then
    say "QUEUE KA-R (original rnup630 on c4h s0 dumps)"
    bash experiments/602-council-gate/gate602.sh "657-rnuporig-s0" 7 -- uv run python $D/rnup657orig.py 0 b_k100_q20r3p10 --src "$D/work/c4h" \
      --name c4h --tag "$TB" >> "$D/work/rnuporig-s0.log" 2>&1
    uv run python $D/ka657.py R 0 >> "$LOG" 2>&1 || { say "KA-R FAILED: stop"; say "END"; exit 1; }
  fi
  emit tvc4h_o4 c4h_o4 "$s" || say "WARN emit c4h_o4 s$s"
  rn "$s" "$D/work/c4h_o4" o4 b_k100_q20r3p10 b_k100_q20r3p10_o4 || say "WARN rnup o4 s$s"
done
say "PHASE 1 DONE"; say "END"
