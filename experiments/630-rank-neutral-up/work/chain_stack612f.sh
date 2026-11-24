#!/bin/bash
# 630 on idall3 (the ungated v3 identity hook, the likely October hook): for each seed k in 0,1,2, wait until 21dd4c's
# dump logs "KA idall3d s<k> ... PASS" (chain612f for s0, chain612g for s1-2; s0's line may lack the seed tag) and
# 85 base_s<k>_* blocks sit in h1_idall3d/blocks, then run the v3 UP pass (one gate602 ticket per run). Price after.
# Cap 12 h from start.
set -u
R=~/Code/virtual-cell-challenge-2026
W=$R/experiments/630-rank-neutral-up/work
D=$R/experiments/612-key-gated-hook/work
G=experiments/602-council-gate/gate602.sh
cd $R
END_TS=$(( $(date +%s) + 43200 ))
ka_ok() {   # seed -> 0 if a PASS line for it exists in either log
  local k=$1
  if [ "$k" = "0" ]; then /usr/bin/grep -Eqs "KA idall3d( s0)? .*PASS" $D/logs/chain612f.log $D/logs/chain612g.log
  else /usr/bin/grep -Eqs "KA idall3d s$k .*PASS" $D/logs/chain612g.log; fi
}
ka_fail() { /usr/bin/grep -Eqs "KA idall3d( s$1)? .*FAIL" $D/logs/chain612f.log $D/logs/chain612g.log; }
for k in 0 1 2; do
  until ka_ok $k && [ "$(ls $D/h1_idall3d/blocks/base_s${k}_* 2>/dev/null | wc -l | tr -d ' ')" = "85" ]; do
    ka_fail $k && { echo "$(date -u +%FT%TZ) KA FAIL s$k; stopping"; break 2; }
    [ "$(date +%s)" -gt "$END_TS" ] && { echo "$(date -u +%FT%TZ) cap reached at s$k"; break 2; }
    sleep 60
  done
  echo "$(date -u +%FT%TZ) idall3d s$k ready"
  bash $G 630-stack-idall3-s$k 6 -- uv run python experiments/630-rank-neutral-up/rnup630.py $k b_k100_q20r3p10 --src $D/h1_idall3d --name idall3 > $W/rnup-idall3-s$k.log 2>&1
  echo "$(date -u +%FT%TZ) idall3 s$k rc $?"
done
bash $G 630-stack-idall3-price 4 -- uv run python experiments/630-rank-neutral-up/price630.py --base $D/h1_idall3d --name idall3 b_k100_q20r3p10 > $W/price630-idall3.log 2>&1
echo "$(date -u +%FT%TZ) price rc $?"
echo "$(date -u +%FT%TZ) END"
