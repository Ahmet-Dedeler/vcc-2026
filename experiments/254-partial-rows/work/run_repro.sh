#!/bin/bash
cd ~/Code/virtual-cell-challenge-2026
for m in partial composed; do for f in k562_gwps orion_hct116 orion_hek293t kolf2; do
  /usr/bin/time -l uv run python experiments/254-partial-rows/repro254.py --mode $m --fold $f > experiments/254-partial-rows/work/repro254-$m-$f.log 2>&1
  echo "$m $f rc=$? $(grep 'maximum resident' experiments/254-partial-rows/work/repro254-$m-$f.log | awk '{printf "%.2f GiB", $1/2^30}')"
done; done
uv run python experiments/254-partial-rows/repro254.py --summary > experiments/254-partial-rows/work/repro254-summary.stdout 2>&1; echo "summary rc=$?"
