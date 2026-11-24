#!/usr/bin/env python3
"""599 known answer: mode ka (run599's swaps, with ps_bulk's own bulk and no sharpening) must reproduce 435's base run:
every array of the stats npz (pre- and post-zs) bit for bit, and the raw metrics to 1e-12.
AMENDED before any arm run (2026-10-02 17:35Z): the first KA run had every array bit-equal but the pre-zs raw nmae at
0.9427040517722616 vs 435's ...618. The scorer itself is not bit-deterministic on identical inputs: 435's own saved arrays
score ...617 in three fresh processes (hashseed599.py, work/hashseed599.out), 2.2e-16 apart. The construction is what
the KA tests, so the arrays stay exact and the metrics get the 1e-12 tolerance. Exit 0 = PASS.
    uv run python experiments/599-tercilec-h1-de/ka599.py [seed]"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
s = int(sys.argv[1]) if len(sys.argv) > 1 else 0
TAG = f"base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse"
A, B = HERE / "work/h1_ka", ROOT / "experiments/435-zs-on-coarse/work/h1"
ok = True
for t in (TAG, TAG + "_dnk250zs"):
    za, zb = np.load(A / f"stats-{t}.npz", allow_pickle=True), np.load(B / f"stats-{t}.npz", allow_pickle=True)
    keys_ok = sorted(za.files) == sorted(zb.files)
    eq = {k: bool(np.array_equal(za[k], zb[k], equal_nan=np.issubdtype(np.asarray(za[k]).dtype, np.floating)))
          for k in zb.files if k in za.files}
    ma = json.loads((A / f"metrics-{t}.json").read_text())["raw"]
    mb = json.loads((B / f"metrics-{t}.json").read_text())["raw"]
    mclose = set(ma) == set(mb) and all(abs(float(ma[k]) - float(mb[k])) <= 1e-12 for k in mb)
    good = keys_ok and all(eq.values()) and mclose
    ok &= good
    print(f"[ka599] {t}: arrays {'all equal' if all(eq.values()) else eq} keys {keys_ok} raw metrics "
          f"{'equal' if ma == mb else ('within 1e-12' if mclose else (ma, mb))} -> {'PASS' if good else 'FAIL'}")
print("KA599 PASS" if ok else "KA599 FAIL")
sys.exit(0 if ok else 1)
