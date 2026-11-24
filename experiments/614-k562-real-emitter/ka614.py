#!/usr/bin/env python3
"""614 wrapper KA: run614 in h1ka mode (no swaps) must reproduce 435's H1 stats (the same run496 path, production
emitter, base s0 bps prodem coarse MDE, and its zs pass). Pass = sums bit-equal and p / q / lfc equal to 1e-12."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
A = ROOT / "experiments/435-zs-on-coarse/work/h1"
B = Path(__file__).resolve().parent / "work/kah1"
ok = True
for tag in ("stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz", "stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz"):
    a, b = np.load(A / tag), np.load(B / tag)
    se = np.array_equal(a["sums"], b["sums"])
    dm = max(float(np.nanmax(np.abs(np.asarray(a[k], float) - np.asarray(b[k], float)))) for k in ("p", "q", "lfc"))
    nan_same = all(np.array_equal(np.isnan(np.asarray(a[k], float)), np.isnan(np.asarray(b[k], float))) for k in ("p", "q", "lfc"))
    good = se and dm <= 1e-12 and nan_same
    ok &= good
    print(f"{tag}: sums {'equal' if se else 'DIFFER'}; max |Δ| p/q/lfc {dm:.2e}; NaN pattern {'same' if nan_same else 'DIFFERS'}: {'PASS' if good else 'FAIL'}")
print(f"KA run614 (h1ka) == 435: {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
