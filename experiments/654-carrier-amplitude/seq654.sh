#!/bin/bash
# 654's runs, serially, one gate602 ticket per run (PREREGISTRATION b0bbfa20). KA first; a failed KA ends the sequence.
# Every exit path writes a terminal "<ts> END" line.
set -u
cd ~/Code/virtual-cell-challenge-2026 || exit 1
D=experiments/654-carrier-amplitude; LOG=$D/work/seq654.log
ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
say() { echo "$(ts) $*" >> "$LOG"; }
run() {  # design seed
  say "QUEUE $1 s$2"
  bash experiments/602-council-gate/gate602.sh "654-$1-s$2" 8 -- uv run python $D/run654.py cpm "$1" "$2" \
    --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs >> "$D/work/run654-$1-s$2.log" 2>&1
  rc=$?; say "DONE $1 s$2 rc=$rc"; return $rc
}
say "START seq654 pid $$"
run base 0 || { say "base s0 failed"; say "END"; exit 1; }
uv run python $D/ka654.py 0 >> "$LOG" 2>&1 || { say "KA s0 FAILED: stop"; say "END"; exit 1; }
for st in "x080 0" "x125 0" "base 1" "x080 1" "x125 1" "x080 2" "x125 2"; do
  set -- $st
  run "$1" "$2" || say "WARN $1 s$2 rc!=0"
  if [ "$1" = base ]; then uv run python $D/ka654.py "$2" >> "$LOG" 2>&1 || { say "KA s$2 FAILED: stop"; say "END"; exit 1; }; fi
done
say "ALL RUNS DONE"
say "END"
