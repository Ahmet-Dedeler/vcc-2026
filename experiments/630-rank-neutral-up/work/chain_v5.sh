#!/bin/bash
# 630 v5 DOWN-currency screen: b_k100_q20r3p10 with the DOWN-gene compensation currency (all / called-only DOWN genes;
# lfc floor x1.5 / x2.5), H1 seeds 0-2, then six-member pricing vs the v3 incumbent. One gate602 ticket PER RUN (v3 rule).
set -u
R=~/Code/virtual-cell-challenge-2026
W=$R/experiments/630-rank-neutral-up/work
G=experiments/602-council-gate/gate602.sh
cd $R
V=b_k100_q20r3p10
for s in 0 1 2; do
  bash $G 630-v5-s$s 6 -- /usr/bin/time -l uv run python experiments/630-rank-neutral-up/rnup630.py $s ${V}_d15 ${V}_d25 ${V}_d15c ${V}_d25c > $W/rnup-v5-s$s.log 2>&1
  echo "$(date -u +%FT%TZ) v5 s$s rc $? $(grep -o '"total_calls_changed": [0-9]*' $W/rnup-v5-s$s.log | tr '\n' ' ')"
done
bash $G 630-v5-price 4 -- uv run python experiments/630-rank-neutral-up/price630.py ${V} ${V}_d15 ${V}_d25 ${V}_d15c ${V}_d25c > $W/price630-v5.log 2>&1
echo "$(date -u +%FT%TZ) price rc $?"
echo "$(date -u +%FT%TZ) END"
