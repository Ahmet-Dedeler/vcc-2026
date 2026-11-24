#!/bin/bash
# 625f: oracle-direction calibration runs (oracle625.py), seed 0, octiles 0 (UP prior, block + hook) and 5 (DOWN prior,
# block only), one job at a time through gate602.
cd ~/Code/virtual-cell-challenge-2026 || exit 1
ARGS="base 0 --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
for MODE in o0 o5; do
  echo "$(date -u +%FT%TZ) start $MODE seed 0"
  bash experiments/602-council-gate/gate602.sh "625f-$MODE-s0" 8 -- uv run python experiments/625-module-readout/oracle625.py $MODE $ARGS \
    > "experiments/625-module-readout/work/oracle625_${MODE}_s0.log" 2>&1
  echo "$(date -u +%FT%TZ) end $MODE seed 0 rc $?"
done
echo "$(date -u +%FT%TZ) END"
