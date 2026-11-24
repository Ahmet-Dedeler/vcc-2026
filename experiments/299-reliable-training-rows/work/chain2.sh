#!/bin/zsh
fold=$1
cd ~/Code/virtual-cell-challenge-2026/experiments/299-reliable-training-rows
while pgrep -f "run299.py $fold" > /dev/null; do sleep 15; done
PYTHONDONTWRITEBYTECODE=1 uv run python run299.py $fold wgt:steeprev wgt:steep >> work/run-$fold.nohup 2>&1
