#!/bin/bash
# 258 -- run the four PLOLO folds one process at a time (memory), each gated on the machine's state.
set -u
cd ~/Code/virtual-cell-challenge-2026
W=experiments/258-union-dose-on-cd4svd5/work
for held in k562_gwps orion_hct116 orion_hek293t kolf2; do
  # gate: no >= 15 GiB job on the machine, and free memory >= 40 %
  while true; do
    big=$(ps aux | awk '$6/1048576 >= 15 {n++} END {print n+0}')
    freepct=$(memory_pressure | awk '/free percentage/ {gsub("%","",$5); print $5}')
    if [ "$big" -eq 0 ] && [ "${freepct:-0}" -ge 40 ]; then break; fi
    echo "[$(date -u +%H:%M:%SZ)] waiting: big=$big free=$freepct%" >> $W/run258.log
    sleep 60
  done
  echo "[$(date -u +%H:%M:%SZ)] fold $held start" >> $W/run258.log
  /usr/bin/time -l uv run python experiments/258-union-dose-on-cd4svd5/price258.py --fold $held > $W/fold-$held.stdout 2> $W/fold-$held.stderr
  rc=$?
  echo "[$(date -u +%H:%M:%SZ)] fold $held exit $rc; $(grep 'maximum resident' $W/fold-$held.stderr | awk '{printf "peak %.2f GiB", $1/1073741824}')" >> $W/run258.log
done
echo "[$(date -u +%H:%M:%SZ)] all folds done" >> $W/run258.log
