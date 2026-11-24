"""811 mechanism read R_tgt (PROPOSAL-2 §4): unweighted mean over targets of the per-target rendered cosine x/sqrt(a d)
(frame-averaged), arm minus base, per fold, from saved STANDARD-FORMAT arrays (light; < 20 MB read).
    uv run python experiments/811-lolo-model/rtgt811.py <run-dir> <arm> <base>"""
import sys
from pathlib import Path

import numpy as np


def cos_t(z):
    x, a, d = np.asarray(z["x"], float), np.asarray(z["a"], float), np.asarray(z["d"], float)
    with np.errstate(all="ignore"):
        c = x / np.sqrt(a * d[None, :])
    return np.nanmean(c, 0)


run, arm, base = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
for d in sorted((run / "bulk").glob("*")):
    if (d / f"{arm}.npz").exists() and (d / f"{base}.npz").exists():
        ca, cb = cos_t(np.load(d / f"{arm}.npz")), cos_t(np.load(d / f"{base}.npz"))
        ok = np.isfinite(ca) & np.isfinite(cb)
        print(f"{d.name:22s} R_tgt {np.mean(ca[ok] - cb[ok]):+.5f}  (n {ok.sum()}, share of targets up {np.mean(ca[ok] > cb[ok]):.2f})")
