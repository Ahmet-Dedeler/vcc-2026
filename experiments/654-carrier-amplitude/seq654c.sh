#!/bin/bash
# 654 dumps for [d4cd5f] (UPDN priced on top of x1.25): x125 with --zero-repair -5 --dump 85, seeds 0-2, MODE cpm-n2200
# (identical code path to cpm: N_CPM 2200 is the default; separate stats dir work/h1_cpm-n2200). Blocks (pre-zs) go to
# work/blocks/. Starts after seq654b's process exits. KA per seed: stats == work/h1_cpm's x125 zr stats, bit for bit.
set -u
cd ~/Code/virtual-cell-challenge-2026 || exit 1
D=experiments/654-carrier-amplitude; LOG=$D/work/seq654c.log
ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
say() { echo "$(ts) $*" >> "$LOG"; }
say "START seq654c pid $$ (waits for seq654b pid $(cat $D/work/seq654b.pid))"
while kill -0 "$(cat $D/work/seq654b.pid)" 2>/dev/null; do sleep 30; done
for s in 0 1 2; do
  say "QUEUE x125 zr dump s$s"
  bash experiments/602-council-gate/gate602.sh "654dump-x125-s$s" 8 -- uv run python $D/run654.py cpm-n2200 x125 "$s" \
    --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs --zero-repair -5 --dump 85 \
    >> "$D/work/run654dump-x125-s$s.log" 2>&1
  say "DONE x125 zr dump s$s rc=$?"
  uv run python - "$s" >> "$LOG" 2>&1 <<'PY' || say "KA dump s$s FAILED"
import sys, glob, numpy as np
s = sys.argv[1]; D = "experiments/654-carrier-amplitude/work"
a = glob.glob(f"{D}/h1_cpm-n2200/stats-x125_s{s}_bps_zr-5_*_dnk250zs.npz"); b = glob.glob(f"{D}/h1_cpm/stats-x125_s{s}_bps_zr-5_*_dnk250zs.npz")
assert len(a) == 1 and len(b) == 1, (a, b)
za, zb = np.load(a[0]), np.load(b[0])
ok = all(np.array_equal(za[k], zb[k], equal_nan=True) if za[k].dtype.kind == "f" else np.array_equal(za[k], zb[k]) for k in ("q", "p", "lfc", "sums"))
print(f"  dump run s{s} vs h1_cpm x125 zr: {'KA PASS' if ok else 'KA FAIL'}"); sys.exit(0 if ok else 1)
PY
done
say "ALL DUMPS DONE"; say "END"
