#!/bin/bash
# waits for queue2's chain to exit, then runs queue3 (only if queue2 ended with QUEUE DONE)
while kill -0 68415 2>/dev/null; do sleep 20; done
grep -q 'QUEUE DONE' "~/Code/virtual-cell-challenge-2026/experiments/436-october-essential-signs/work/drive436-q2.log" || { echo "queue2 did not finish cleanly: not starting queue3"; exit 1; }
exec "~/Code/virtual-cell-challenge-2026/experiments/436-october-essential-signs/drive436.sh" "~/Code/virtual-cell-challenge-2026/experiments/436-october-essential-signs/work/queue3.txt"
