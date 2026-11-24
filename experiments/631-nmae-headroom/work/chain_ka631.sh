#!/bin/bash
set -u
R=~/Code/virtual-cell-challenge-2026; X=$R/experiments/631-nmae-headroom/work
G=experiments/602-council-gate/gate602.sh; cd $R
echo "$(date -u +%FT%TZ) start module sha $(shasum -a 256 experiments/631-nmae-headroom/src_dn_mag631.py | cut -c1-16)"
for mode in dn updn; do for s in 0 1 2; do
  bash $G 631-ka-$mode-s$s 5 -- uv run python experiments/631-nmae-headroom/test_dn_mag631.py $s $mode 85 > $X/ka631-$mode-s$s.log 2>&1
  echo "$(date -u +%FT%TZ) KA $mode s$s rc $? : $(tail -1 $X/ka631-$mode-s$s.log)"
done; done
echo "$(date -u +%FT%TZ) END"
