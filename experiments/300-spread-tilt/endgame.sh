#!/bin/sh
cd ~/Code/virtual-cell-challenge-2026/experiments/300-spread-tilt
# wait for the in-flight run300 fold processes to finish, then stop the scheduler and re-allocate the two slots
while pgrep -f 'run300[.]py' > /dev/null; do sleep 10; done
pkill -f 'sched300[.]py'
sleep 2
date -u +"%H:%M:%SZ scheduler stopped; launching the mechanism row and the random-weight control" >> work/sched300.log
nohup uv run python mech300.py orion_hct116 tilt:d2.0 spread:d1.0 > work/mech300.nohup 2>&1 &
nohup uv run python run300.py k562_gwps rand:d1.0:s300001 > work/run300-k562_gwps.nohup 2>&1 &
wait
