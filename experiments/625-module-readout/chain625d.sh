#!/bin/bash
# 625d chain: one job at a time through gate602 (FIFO), so this experiment never holds more than one of the 2 slots.
# Order (d4cd5f): the base `none` with dumps first (its UP stack is every octile's reference), then octiles 0 2 5 7; per seed 0, 1, 2.
cd ~/Code/virtual-cell-challenge-2026 || exit 1
ARGS="base SEED --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs --dump 85"
for SEED in 0 1 2; do
  for OCT in none 0 2 5 7; do
    A=${ARGS/SEED/$SEED}
    echo "$(date -u +%FT%TZ) start oct $OCT seed $SEED"
    bash experiments/602-council-gate/gate602.sh "625d-$OCT-s$SEED" 8 -- uv run python experiments/625-module-readout/flip625.py $OCT $A \
      > "experiments/625-module-readout/work/flip625_${OCT}_s${SEED}.log" 2>&1
    echo "$(date -u +%FT%TZ) end oct $OCT seed $SEED rc $?"
  done
done
echo "$(date -u +%FT%TZ) END"
