#!/bin/bash
# 654 Amendment A1 (b1590abd): production-faithful arms with --zero-repair -5, serially, one gate602 ticket per run.
# Structural KA on each base_zr seed. Every exit path writes a terminal "<ts> END" line.
set -u
cd ~/Code/virtual-cell-challenge-2026 || exit 1
D=experiments/654-carrier-amplitude; LOG=$D/work/seq654b.log
ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
say() { echo "$(ts) $*" >> "$LOG"; }
run() {  # design seed
  say "QUEUE $1 zr s$2"
  bash experiments/602-council-gate/gate602.sh "654zr-$1-s$2" 8 -- uv run python $D/run654.py cpm "$1" "$2" \
    --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs --zero-repair -5 >> "$D/work/run654zr-$1-s$2.log" 2>&1
  rc=$?; say "DONE $1 zr s$2 rc=$rc"; return $rc
}
say "START seq654b pid $$"
for s in 0 1 2; do
  run base "$s" || say "WARN base zr s$s rc!=0"
  uv run python $D/kazr654.py "$s" >> "$LOG" 2>&1 || { say "KA zr s$s FAILED: stop"; say "END"; exit 1; }
  run x125 "$s" || say "WARN x125 zr s$s rc!=0"
done
say "ALL ZR RUNS DONE"
say "END"
