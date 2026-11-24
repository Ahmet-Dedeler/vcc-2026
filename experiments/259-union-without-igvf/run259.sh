#!/bin/bash
# 259 -- the PLOLO folds one process at a time (memory), each gated on the machine's state (258's gate).
set -u
cd ~/Code/virtual-cell-challenge-2026
W=experiments/259-union-without-igvf/work
for held in "$@"; do
  while true; do
    big=$(ps aux | awk '$6/1048576 >= 15 {n++} END {print n+0}')
    freepct=$(memory_pressure | awk '/free percentage/ {gsub("%","",$5); print $5}')
    if [ "$big" -eq 0 ] && [ "${freepct:-0}" -ge 40 ]; then break; fi
    echo "[$(date -u +%H:%M:%SZ)] waiting: big=$big free=$freepct%" >> $W/run259.log
    sleep 60
  done
  echo "[$(date -u +%H:%M:%SZ)] fold $held start" >> $W/run259.log
  /usr/bin/time -l uv run python experiments/259-union-without-igvf/price259.py --fold $held > $W/fold-$held.stdout 2> $W/fold-$held.stderr
  rc=$?
  echo "[$(date -u +%H:%M:%SZ)] fold $held exit $rc; $(grep 'maximum resident' $W/fold-$held.stderr | awk '{printf "peak %.2f GiB", $1/1073741824}')" >> $W/run259.log
done
echo "[$(date -u +%H:%M:%SZ)] folds done: $*" >> $W/run259.log
