"""271 -- diagnostic for the NEGATIVE whitened split-half reliability of the donor rows: is it an IN-SAMPLE effect of
the operator being fitted on the very rows it whitens?

For each cancer fold: the fold's operator (fitted on the two DONOR banks, the held line absent) is applied to
  (i)  each donor's half rows A, B  (in-sample: the donor's full-depth rows are in the fit)      -> cos(WA, WB)
  (ii) the HELD line's half rows A, B (out-of-sample: none of its rows are in the fit)           -> cos(WA, WB)
and, for both, the null cos(WA_t, WB_t') over a fixed derangement (a common component would show here), the whitened
half-to-full cosine cos(WA_t, W real_t), and the same three quantities in the raw space.  The held line's half rows
are its DONOR-format rows (work/half/rows-<held>.npz) -- the same objects 267 built its truths from; the reserved
orion_hek293t half-B TRUTHS are not read (the row file is an input built by halfrows271.py from 267's sums, and this
diagnostic never scores anything against it).  Writes work/diagwht271.json.
    PYTHONDONTWRITEBYTECODE=1 uv run python diagwht271.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/271-donor-reliability-estimation"
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/246-direction-hypotheses"))
sys.path.insert(0, str(HERE))
import lib246 as L  # noqa: E402
import plolo  # noqa: E402
import est271 as E  # noqa: E402

plolo.WORK = WORK / "plolo"
L.WORK = WORK
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t")
OUT = {}


def stats(A, B, R, ok, seed=1):
    idx = np.flatnonzero(ok)
    pi = E.derangement(idx.size, seed)
    c = E.row_cos(A[idx], B[idx])
    null = E.row_cos(A[idx], B[idx[pi]])
    cf = E.row_cos(A[idx], R[idx])
    return dict(n=int(idx.size), median=float(np.nanmedian(c)), mean=float(np.nanmean(c)), share_le0=float(np.nanmean(c <= 0)),
                p10=float(np.nanpercentile(c, 10)), p90=float(np.nanpercentile(c, 90)),
                null_median=float(np.nanmedian(null)), null_mean=float(np.nanmean(null)), null_sd=float(np.nanstd(null)),
                half_to_full_median=float(np.nanmedian(cf)))


for held in FOLDS:
    fold = plolo.build_fold(held)
    wc = np.asarray(fold["wcols"], dtype=np.int64)
    op = plolo.Operator(L.op_path_for(held))
    rec = {}
    for line, role in [(s, "donor (in-sample)") for s in plolo.FOLDS[held]["sources"]] + [(held, "held (out-of-sample)")]:
        z = np.load(WORK / f"half/rows-{line}.npz", allow_pickle=True)
        ok = z["avail_A"] & z["avail_B"]
        A, B, R = z["A"].astype(np.float64), z["B"].astype(np.float64), z["real"].astype(np.float64)
        raw = stats(A[:, wc], B[:, wc], R[:, wc], ok)
        WA = op.apply(A * L.ALPHA, wc, L.GAMMA, L.LAM)[:, wc]
        WB = op.apply(B * L.ALPHA, wc, L.GAMMA, L.LAM)[:, wc]
        WR = op.apply(R * L.ALPHA, wc, L.GAMMA, L.LAM)[:, wc]
        wht = stats(WA, WB, WR, ok)
        rec[line] = dict(role=role, raw=raw, wht=wht)
        print(f"[{held}] {line:14s} {role:22s} raw: median {raw['median']:+.4f} share<=0 {raw['share_le0']:.3f} null {raw['null_median']:+.4f} half~full {raw['half_to_full_median']:.4f} | "
              f"wht: median {wht['median']:+.4f} share<=0 {wht['share_le0']:.3f} null {wht['null_median']:+.4f} (sd {wht['null_sd']:.4f}) half~full {wht['half_to_full_median']:.4f}", flush=True)
    OUT[held] = rec
    del op
(WORK / "diagwht271.json").write_text(json.dumps(OUT, indent=1))
print("written work/diagwht271.json")
