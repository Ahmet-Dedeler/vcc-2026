#!/bin/bash
set -u
R=~/Code/virtual-cell-challenge-2026; W=$R/experiments/630-rank-neutral-up/work; X=$R/experiments/631-nmae-headroom/work
G=experiments/602-council-gate/gate602.sh; cd $R
for s in 0 1 2; do
  bash $G 631-v7-s$s 6 -- /usr/bin/time -l uv run python experiments/630-rank-neutral-up/rnup630.py $s b_k100_q20r5p10 b_k100_q20r3p10_P200 b_k100_q20r5p10_P200 > $X/rnup-v7-s$s.log 2>&1
  echo "$(date -u +%FT%TZ) v7 s$s rc $? $(grep -oE '"(dn_paid|sub_paid|total_calls_changed)": [-0-9.e+]+' $X/rnup-v7-s$s.log | tr '\n' ' ')"
done
bash $G 631-v7-price 4 -- uv run python experiments/630-rank-neutral-up/price630.py b_k100_q20r3p10 b_k100_q20r5p10 b_k100_q20r3p10_P200 b_k100_q20r5p10_P200 > $X/price-v7.log 2>&1
echo "$(date -u +%FT%TZ) price rc $?"
echo "$(date -u +%FT%TZ) END"
