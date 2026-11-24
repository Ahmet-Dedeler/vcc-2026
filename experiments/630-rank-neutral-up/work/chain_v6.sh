#!/bin/bash
# 630 v6 PROBE (library drift as compensation): b_k100_q20r3p10_L on H1 seeds 0-2 + pricing. One gate602 ticket per run.
set -u
R=~/Code/virtual-cell-challenge-2026
W=$R/experiments/630-rank-neutral-up/work
G=experiments/602-council-gate/gate602.sh
cd $R
V=b_k100_q20r3p10
for s in 0 1 2; do
  bash $G 630-v6-s$s 6 -- /usr/bin/time -l uv run python experiments/630-rank-neutral-up/rnup630.py $s ${V}_L > $W/rnup-v6-s$s.log 2>&1
  echo "$(date -u +%FT%TZ) v6 s$s rc $? $(grep -oE '"(L_max_dAUC_all|L_n_over_1e4|total_calls_changed|frac_got|L_lib_change_max)": [-0-9.e]+' $W/rnup-v6-s$s.log | tr '\n' ' ')"
done
bash $G 630-v6-price 4 -- uv run python experiments/630-rank-neutral-up/price630.py ${V} ${V}_L > $W/price630-v6.log 2>&1
echo "$(date -u +%FT%TZ) price rc $?"
echo "$(date -u +%FT%TZ) END"
