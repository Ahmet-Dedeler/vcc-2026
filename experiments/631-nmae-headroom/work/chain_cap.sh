#!/bin/bash
set -u
R=~/Code/virtual-cell-challenge-2026; X=$R/experiments/631-nmae-headroom/work; cd $R
for cfg in "1e-3 head" "3e-3 head" "1e-4 called" "1e-3 called"; do
  set -- $cfg
  bash experiments/602-council-gate/gate602.sh 631-cap-$1-$2 5 -- uv run python experiments/631-nmae-headroom/cap631.py 12 $1 $2 > $X/cap631_$1_$2.log 2>&1
  echo "$(date -u +%FT%TZ) $1 $2 rc $? $(grep '^MEAN' $X/cap631_$1_$2.log)"
done
echo "$(date -u +%FT%TZ) END"
