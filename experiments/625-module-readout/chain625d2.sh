#!/bin/bash
# 625d chain, part 2 (2026-10-02 ~13:55Z): the gate602 queue is long (four sessions), so the replicate seeds are cut to the
# two octiles that carry the seed-noise estimate: one UP-prior octile (0, hook + block) and one DOWN-prior octile (5, block
# only, the call-count side effect). Part 1 (chain625d.sh) ran none/0/2/5 seed 0 and was stopped while its octile-7 seed-0
# job (gate602 pid $1) was queued. That job runs on by itself; this chain waits for it to exit first, so 625's jobs stay serial.
# usage: nohup bash experiments/625-module-readout/chain625d2.sh <pid of the queued 625d-7-s0 gate job> > .../chain625d2.log 2>&1 &
cd ~/Code/virtual-cell-challenge-2026 || exit 1
WAITPID=$1
ARGS="base SEED --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs --dump 85"
echo "$(date -u +%FT%TZ) waiting for pid $WAITPID (625d-7-s0) to exit"
end=$(( $(date +%s) + 6 * 3600 ))
while kill -0 "$WAITPID" 2>/dev/null; do
  [ "$(date +%s)" -ge "$end" ] && { echo "$(date -u +%FT%TZ) gave up waiting for $WAITPID after 6 h"; echo "$(date -u +%FT%TZ) END"; exit 1; }
  sleep 30
done
echo "$(date -u +%FT%TZ) pid $WAITPID gone; gate log: $(grep '625d-7-s0' experiments/602-council-gate/work/gate602.log | tail -1)"
for SEED in 1 2; do
  for OCT in none 0 5; do
    A=${ARGS/SEED/$SEED}
    echo "$(date -u +%FT%TZ) start oct $OCT seed $SEED"
    bash experiments/602-council-gate/gate602.sh "625d-$OCT-s$SEED" 8 -- uv run python experiments/625-module-readout/flip625.py $OCT $A \
      > "experiments/625-module-readout/work/flip625_${OCT}_s${SEED}.log" 2>&1
    echo "$(date -u +%FT%TZ) end oct $OCT seed $SEED rc $?"
  done
done
echo "$(date -u +%FT%TZ) END"
