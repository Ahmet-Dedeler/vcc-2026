#!/bin/bash
# waits for queue1's driver to exit, then runs queue2 (only if queue1 ended with QUEUE DONE)
while kill -0 67105 2>/dev/null; do sleep 20; done
grep -q 'QUEUE DONE' "~/Code/virtual-cell-challenge-2026/experiments/436-october-essential-signs/work/drive436-q1.log" || { echo "queue1 did not finish cleanly: not starting queue2"; exit 1; }
exec "~/Code/virtual-cell-challenge-2026/experiments/436-october-essential-signs/drive436.sh" "~/Code/virtual-cell-challenge-2026/experiments/436-october-essential-signs/work/queue2.txt"
