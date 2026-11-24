#!/bin/bash
# 657 Phase 2 (PREREGISTRATION ed36f45f): the DN stacks (b_k100_q20r3p10_D200, d4cd5f's frozen form) with / without _o4 on
# BOTH emissions, seeds 0-2. Starts after Phase 1's process exits. Two chains in parallel (c4h dumps | c4h_o4 dumps), each
# serial, one gate602 ticket per run. Terminal "<ts> END" lines per chain and overall.
set -u
cd ~/Code/virtual-cell-challenge-2026 || exit 1
D=experiments/657-octile4-up-callneutral; LOG=$D/work/seq657b.log
TB="base_s{seed}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse"
ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
say() { echo "$(ts) $*" >> "$LOG"; }
say "START seq657b pid $$ (waits for Phase 1 pid $(cat $D/work/seq657.pid))"
P1=$(cat $D/work/seq657.pid)
while kill -0 "$P1" 2>/dev/null; do sleep 60; done
grep -q "PHASE 1 DONE" $D/work/seq657.log || { say "Phase 1 did not finish cleanly: no Phase 2"; say "END"; exit 1; }
chain() {  # src name
  local src=$1 nm=$2
  for s in 0 1 2; do
    say "QUEUE rnup DN $nm s$s"
    bash experiments/602-council-gate/gate602.sh "657dn-$nm-s$s" 7 -- uv run python $D/rnup657.py "$s" b_k100_q20r3p10_D200 \
      b_k100_q20r3p10_D200_o4 --src "$src" --name "$nm" --tag "$TB" >> "$D/work/rnupdn-$nm-s$s.log" 2>&1
    say "DONE rnup DN $nm s$s rc=$?"
  done
  say "CHAIN $nm END"
}
chain "$D/work/c4h" c4h &
chain "$D/work/c4h_o4" o4 &
wait
uv run python $D/analyze657.py >> "$LOG" 2>&1
say "PHASE 2 DONE"; say "END"
