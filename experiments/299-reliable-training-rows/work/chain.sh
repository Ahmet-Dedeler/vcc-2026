#!/bin/zsh
fold=$1; pid=$2
while kill -0 $pid 2>/dev/null; do sleep 20; done
cd ~/Code/virtual-cell-challenge-2026/experiments/299-reliable-training-rows
PYTHONDONTWRITEBYTECODE=1 uv run python run299.py $fold wgt:mild wgt:mildrev wgt:R >> work/run-$fold.nohup 2>&1
