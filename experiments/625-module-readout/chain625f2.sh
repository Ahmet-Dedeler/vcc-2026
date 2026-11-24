#!/bin/bash
# 625f part 2: oracle-direction arms for octiles 2 (UP prior, block only) and 7 (strongest DOWN prior, call-set term),
# seed 0, one gate602 ticket per run. Started after chain625f.sh (o0, o5) has exited.
cd ~/Code/virtual-cell-challenge-2026 || exit 1
WAITPID=$1
end=$(( $(date +%s) + 4 * 3600 ))
while kill -0 "$WAITPID" 2>/dev/null; do
  [ "$(date +%s)" -ge "$end" ] && { echo "$(date -u +%FT%TZ) gave up waiting for $WAITPID"; echo "$(date -u +%FT%TZ) END"; exit 1; }
  sleep 30
done
ARGS="base 0 --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs"
for MODE in o2 o7; do
  echo "$(date -u +%FT%TZ) start $MODE seed 0"
  bash experiments/602-council-gate/gate602.sh "625f-$MODE-s0" 8 -- uv run python experiments/625-module-readout/oracle625.py $MODE $ARGS \
    > "experiments/625-module-readout/work/oracle625_${MODE}_s0.log" 2>&1
  echo "$(date -u +%FT%TZ) end $MODE seed 0 rc $?"
done
echo "$(date -u +%FT%TZ) END"
