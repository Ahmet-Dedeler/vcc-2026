#!/bin/zsh
# 264: four PLOLO folds, ONE fold per process, two processes at a time (peak ~5.5 GiB each on a 64 GiB machine).
# Pass 1: the mandate's construction (equal three-donor mean) on all four folds.  Pass 2: the matched sensitivity on the cancer folds.
cd ~/Code/virtual-cell-challenge-2026/experiments/264-majority-rule-cancer-folds
export PYTHONDONTWRITEBYTECODE=1
for pair in "k562_gwps orion_hct116" "orion_hek293t kolf2"; do
  for f in ${=pair}; do
    uv run python sweep264.py --fold $f > work/sweep264-$f.log 2>&1 &
    sleep 15
  done
  wait
done
touch work/sweep264-folds.done
for pair in "k562_gwps orion_hct116" "orion_hek293t"; do
  for f in ${=pair}; do
    uv run python sweep264.py --fold $f --kolf-weight matched > work/sweep264m-$f.log 2>&1 &
    sleep 15
  done
  wait
done
touch work/sweep264m-folds.done
