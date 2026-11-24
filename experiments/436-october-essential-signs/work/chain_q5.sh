#!/bin/bash
# waits for queue4's driver to exit, then runs queue5 (only if queue4 ended with QUEUE DONE)
while kill -0 3945 2>/dev/null; do sleep 20; done
grep -q 'QUEUE DONE' "~/Code/virtual-cell-challenge-2026/experiments/436-october-essential-signs/work/drive436-q4.log" || { echo "queue4 did not finish cleanly: not starting queue5"; exit 1; }
exec "~/Code/virtual-cell-challenge-2026/experiments/436-october-essential-signs/drive436.sh" "~/Code/virtual-cell-challenge-2026/experiments/436-october-essential-signs/work/queue5.txt"
