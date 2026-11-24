#!/bin/zsh
# 263: four PLOLO folds, ONE fold per process, two processes at a time (peak ~5.5 GiB each on a 64 GiB machine).
cd ~/Code/virtual-cell-challenge-2026/experiments/263-donor-selection-arbiters
for pair in "k562_gwps orion_hct116" "orion_hek293t kolf2"; do
  for f in ${=pair}; do
    uv run python sweep263.py --fold $f > work/sweep263-$f.log 2>&1 &
    sleep 20
  done
  wait
done
touch work/sweep263-folds.done
