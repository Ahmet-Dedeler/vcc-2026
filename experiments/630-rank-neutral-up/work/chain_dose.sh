#!/bin/bash
# 630 v3 dose bracket: K 0.5 and K 1.5 (q20 r3 p10) on seeds 0-2, so the board dose pair can be chosen where H1 separates.
set -u
R=~/Code/virtual-cell-challenge-2026
W=$R/experiments/630-rank-neutral-up/work
cd $R
for s in 0 1 2; do
  /usr/bin/time -l uv run python experiments/630-rank-neutral-up/rnup630.py $s b_k050_q20r3p10 b_k150_q20r3p10 > $W/rnup-dose-s$s.log 2>&1
  echo "$(date -u +%FT%TZ) dose s$s rc $?"
done
uv run python experiments/630-rank-neutral-up/price630.py b_k050_q20r3p10 b_k075_q20r3p10 b_k100_q20r3p10 b_k150_q20r3p10 > $W/price630-dose.log 2>&1
echo "$(date -u +%FT%TZ) price rc $?"
echo "$(date -u +%FT%TZ) END"
