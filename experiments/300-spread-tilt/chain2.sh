#!/bin/sh
cd ~/Code/virtual-cell-challenge-2026/experiments/300-spread-tilt
while pgrep -f 'run300[.]py orion_hek293t' > /dev/null; do sleep 10; done
uv run python run300.py orion_hct116 rand:d1.0:s300001 >> work/run300-orion_hct116.nohup 2>&1
uv run python run300.py k562_gwps rand:d1.0:s300001 >> work/run300-k562_gwps.nohup 2>&1
