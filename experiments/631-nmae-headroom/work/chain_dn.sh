#!/bin/bash
# 631 DOWN pass on H1: DOWN-only (b_k000..._D200) and UP+DOWN (b_k100_q20r3p10_D200), seeds 0-2, then pricing.
set -u
R=~/Code/virtual-cell-challenge-2026; W=$R/experiments/630-rank-neutral-up/work; X=$R/experiments/631-nmae-headroom/work
G=experiments/602-council-gate/gate602.sh; cd $R
for s in 0 1 2; do
  bash $G 631-dn-s$s 6 -- /usr/bin/time -l uv run python experiments/630-rank-neutral-up/rnup630.py $s b_k000_q20r3p10_D200 b_k100_q20r3p10_D200 > $X/rnup-dn-s$s.log 2>&1
  echo "$(date -u +%FT%TZ) dn s$s rc $? $(grep -oE '"(dn_frac|total_calls_changed)": [-0-9.e]+' $X/rnup-dn-s$s.log | tr '\n' ' ')"
done
bash $G 631-dn-price 4 -- uv run python experiments/630-rank-neutral-up/price630.py b_k100_q20r3p10 b_k000_q20r3p10_D200 b_k100_q20r3p10_D200 > $X/price-dn.log 2>&1
echo "$(date -u +%FT%TZ) price rc $?"
echo "$(date -u +%FT%TZ) END"
