#!/bin/bash
# 631 DN-alone re-price with a TRUE K=0 (no UP pass), seeds 0-2, then the dn-mode module KA and pricing.
set -u
R=~/Code/virtual-cell-challenge-2026; X=$R/experiments/631-nmae-headroom/work
G=experiments/602-council-gate/gate602.sh; cd $R
for s in 0 1 2; do
  bash $G 631-dn3-s$s 6 -- uv run python experiments/630-rank-neutral-up/rnup630.py $s b_k000_q20r3p10_D200 > $X/rnup-dn3-s$s.log 2>&1
  echo "$(date -u +%FT%TZ) dn3 s$s rc $? $(grep -oE '"(dn_frac|n_moved|total_calls_changed)": [-0-9.e]+' $X/rnup-dn3-s$s.log | tr '\n' ' ')"
  bash $G 631-ka3-dn-s$s 5 -- uv run python experiments/631-nmae-headroom/test_dn_mag631.py $s dn 85 > $X/ka631-dn3-s$s.log 2>&1
  echo "$(date -u +%FT%TZ) KA dn s$s rc $? : $(grep -h 'PASS\|FAIL' $X/ka631-dn3-s$s.log)"
done
bash $G 631-dn3-price 4 -- uv run python experiments/630-rank-neutral-up/price630.py b_k100_q20r3p10 b_k000_q20r3p10_D200 b_k100_q20r3p10_D200 > $X/price-dn3.log 2>&1
echo "$(date -u +%FT%TZ) price rc $?"
echo "$(date -u +%FT%TZ) END"
