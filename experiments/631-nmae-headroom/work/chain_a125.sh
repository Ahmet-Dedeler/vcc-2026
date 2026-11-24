#!/bin/bash
# 631 stacking test on 654's x1.25 + zero-repair construction: 630 v3, 630 + DN, DN alone, priced against the x1.25 zs base.
set -u
R=~/Code/virtual-cell-challenge-2026; X=$R/experiments/631-nmae-headroom/work; S=$X/a125
TAG="x125_s{seed}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse"
G=experiments/602-council-gate/gate602.sh; cd $R
for s in 0 1 2; do
  bash $G 631-a125-s$s 6 -- /usr/bin/time -l uv run python experiments/630-rank-neutral-up/rnup630.py $s b_k100_q20r3p10 b_k100_q20r3p10_D200 b_k000_q20r3p10_D200 --src $S --name a125 --tag "$TAG" > $X/rnup-a125-s$s.log 2>&1
  echo "$(date -u +%FT%TZ) a125 s$s rc $? $(grep -oE '"(dn_frac|total_calls_changed|zsKA_matches_run630_metrics)": [a-z0-9.e-]+' $X/rnup-a125-s$s.log | tr '\n' ' ')"
done
bash $G 631-a125-price 4 -- uv run python experiments/630-rank-neutral-up/price630.py --base $S --name a125 --tag "$TAG" b_k100_q20r3p10 b_k100_q20r3p10_D200 b_k000_q20r3p10_D200 > $X/price-a125.log 2>&1
echo "$(date -u +%FT%TZ) price rc $?"
echo "$(date -u +%FT%TZ) END"
