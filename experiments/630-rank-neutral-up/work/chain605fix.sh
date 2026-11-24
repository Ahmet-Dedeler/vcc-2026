#!/bin/bash
# 630 post-605-fix verification: (1) the fixed module reproduces every priced H1 run bit for bit (seeds 0-2 x the v2
# incumbent b_k100 and the v3 dose pair) => the pre-fix code never used a negative-gain donor or lowered a mean on H1;
# (2) A/B/C capacity with the fixed module (counters + per-target identity vs the pre-fix json). Every step logs; END last.
set -u
R=~/Code/virtual-cell-challenge-2026
W=$R/experiments/630-rank-neutral-up/work
cd $R
echo "$(date -u +%FT%TZ) start module sha $(shasum -a 256 experiments/630-rank-neutral-up/src_rank_up.py | cut -c1-16)"
for v in b_k100_q20r3p10 b_k075_q20r3p10 b_k100; do
  for s in 0 1 2; do
    uv run python experiments/630-rank-neutral-up/test_rank_up630.py $s $v 85 > $W/ka605-$v-s$s.log 2>&1
    echo "$(date -u +%FT%TZ) KA $v s$s rc $? : $(tail -1 $W/ka605-$v-s$s.log)"
  done
done
VARIANT=b_k100_q20r3p10 bash $W/capacity_chain.sh > $W/capacity_chain_v3fix.log 2>&1
echo "$(date -u +%FT%TZ) capacity rc $? : $(grep '^[ABC]:' $W/capacity_chain_v3fix.log | tr '\n' ' ')"
echo "$(date -u +%FT%TZ) END"
