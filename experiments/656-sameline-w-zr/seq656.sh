#!/bin/bash
# 656 (PREREGISTRATION 28d99181): after 655's chain exits, 18 zero-repair runs of 450's weaker-bank w arms via gate602,
# seed 0 first with the structural KA per arm, then seeds 1-2; then the reader. Terminal "<ts> END" on every path.
set -u
cd ~/Code/virtual-cell-challenge-2026 || exit 1
D=experiments/656-sameline-w-zr; LOG=$D/work/seq656.log
ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
say() { echo "$(ts) $*" >> "$LOG"; }
say "START seq656 pid $$"
P655=$(cat experiments/655-carrier-amplitude-dose/work/seq655.pid 2>/dev/null)
while [ -n "$P655" ] && kill -0 "$P655" 2>/dev/null; do sleep 60; done
A="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zero-repair -5"
run() {  # w p seed
  say "QUEUE w$1 p$2 s$3"
  RUN450_FLIP=$2 RUN450_NOISE_SEED=0 RUN426_EMITTER=prod bash experiments/602-council-gate/gate602.sh "656-w$1-p$2-s$3" 8 -- \
    uv run python $D/run656c.py "$1" clean4 base "$3" $A >> "$D/work/run656-w$1-p$2-s$3.log" 2>&1
  r=$?; say "DONE w$1 p$2 s$3 rc=$r"; return $r
}
for p in 0.2 0.3; do for w in 0.5 1 2; do
  run "$w" "$p" 0 || say "WARN w$w p$p s0 rc!=0"
  uv run python $D/ka656.py "$w" "$p" >> "$LOG" 2>&1 || { say "KA w$w p$p FAILED: stop"; say "END"; exit 1; }
done; done
for s in 1 2; do for p in 0.2 0.3; do for w in 0.5 1 2; do run "$w" "$p" "$s" || say "WARN w$w p$p s$s rc!=0"; done; done; done
uv run python $D/read656.py >> "$LOG" 2>&1
say "ALL 656 RUNS DONE"; say "END"
