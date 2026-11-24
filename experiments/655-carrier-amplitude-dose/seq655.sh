#!/bin/bash
# 655 (PREREGISTRATION 18f4d2cb): after 654's zr arms + dumps exit, check the run condition; then KA (x125 zr s0 bit for bit
# vs 654) and x150 / x200 with --zero-repair -5, seeds 0-2, one gate602 ticket per run. Terminal "<ts> END" on every path.
set -u
cd ~/Code/virtual-cell-challenge-2026 || exit 1
D=experiments/655-carrier-amplitude-dose; D4=experiments/654-carrier-amplitude; LOG=$D/work/seq655.log
ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
say() { echo "$(ts) $*" >> "$LOG"; }
say "START seq655 pid $$"
for p in "$(cat $D4/work/seq654.pid)" "$(cat $D4/work/seq654b.pid)" "$(cat $D4/work/seq654c.pid)"; do
  while kill -0 "$p" 2>/dev/null; do sleep 30; done
done
uv run python $D4/analyze654.py --zr > /dev/null 2>&1
uv run python $D/cond655.py >> "$LOG" 2>&1; rc=$?
if [ $rc -ne 0 ]; then say "WITHDRAWN (condition not met): no 655 run"; say "END"; exit 0; fi
run() {  # design seed
  say "QUEUE $1 zr s$2"
  bash experiments/602-council-gate/gate602.sh "655-$1-s$2" 8 -- uv run python $D/run655.py cpm "$1" "$2" \
    --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs --zero-repair -5 >> "$D/work/run655-$1-s$2.log" 2>&1
  r=$?; say "DONE $1 zr s$2 rc=$r"; return $r
}
run x125 0 || { say "KA run failed"; say "END"; exit 1; }
uv run python - >> "$LOG" 2>&1 <<'PY' || { say "KA FAILED: stop"; say "END"; exit 1; }
import glob, sys, numpy as np
a = glob.glob("experiments/655-carrier-amplitude-dose/work/h1_cpm/stats-x125_s0_bps_zr-5_*_dnk250zs.npz")
b = glob.glob("experiments/654-carrier-amplitude/work/h1_cpm/stats-x125_s0_bps_zr-5_*_dnk250zs.npz")
assert len(a) == 1 and len(b) == 1, (a, b)
za, zb = np.load(a[0]), np.load(b[0])
ok = all(np.array_equal(za[k], zb[k], equal_nan=True) if za[k].dtype.kind == "f" else np.array_equal(za[k], zb[k]) for k in ("q", "p", "lfc", "sums"))
print("  655 KA (x125 zr s0 vs 654):", "KA PASS" if ok else "KA FAIL"); sys.exit(0 if ok else 1)
PY
for st in "x150 0" "x200 0" "x150 1" "x200 1" "x150 2" "x200 2"; do
  set -- $st; run "$1" "$2" || say "WARN $1 s$2 rc!=0"
done
say "ALL 655 RUNS DONE"; say "END"
