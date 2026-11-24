#!/bin/sh
cd ~/Code/virtual-cell-challenge-2026/experiments/300-spread-tilt
export PYTHONDONTWRITEBYTECODE=1
uv run python run300.py orion_hct116 rand:d1.0:s300002 >> work/run300-orion_hct116.nohup 2>&1 &
uv run python run300.py orion_hek293t rand:d1.0:s300002 >> work/run300-orion_hek293t.nohup 2>&1 &
wait
uv run python run300.py k562_gwps rand:d1.0:s300002 >> work/run300-k562_gwps.nohup 2>&1
