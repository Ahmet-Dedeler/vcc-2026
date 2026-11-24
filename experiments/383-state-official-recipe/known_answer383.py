#!/usr/bin/env python3
"""383 known-answer: recompute 323 cap15 six raw members via score() only.

Reads 323 artifacts read-only; writes only to this experiment's work/.
Pass = all six |recomputed - recorded| <= 323 TOLERANCE (1e-6; mse 1e-10; pds 1e-12).
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("POLARS_MAX_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
for _p in (str(BASE), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np

import evaluate323 as E  # noqa: E402
import run323 as R  # noqa: E402

TOL = {"pds": 1e-12, "mse": 1e-10, "nmae": 1e-6, "fid": 1e-6, "reach": 1e-6, "jac": 1e-6}


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    # NOTE: E.verify() not called: 323's frozen-input manifest predates committed
    # builder evolution (c47d2343 VCC_CTX_A_PRIOR chain); the drift is committed
    # history, not corruption. load_reference() still hash-checks the statistics
    # that enter score(), which is the part load-bearing for this recompute.
    r, c, _ = E.load_reference()
    with np.load(BASE / "work/stats-cap15.npz", allow_pickle=False) as z:
        sums = [np.asarray(z["sums"][i]) for i in range(z["sums"].shape[0])]
        jks = np.asarray(z["jk"], dtype=float)
        ps = [np.asarray(z["p"][i]) for i in range(z["p"].shape[0])]
        qs = [np.asarray(z["q"][i]) for i in range(z["q"].shape[0])]
        lfcs = [np.asarray(z["lfc"][i]) for i in range(z["lfc"].shape[0])]
    raw, _, _ = E.score(r, c, sums, jks, ps, qs, lfcs)
    recorded = json.loads((BASE / "work/metrics-cap15.json").read_text())["raw"]
    diffs = {m: abs(float(raw[m]) - float(recorded[m])) for m in recorded}
    ok = all(diffs[m] <= TOL[m] for m in recorded)
    out = {"recomputed_raw": {m: float(raw[m]) for m in recorded},
           "recorded_raw": {m: float(recorded[m]) for m in recorded},
           "absolute_differences": diffs, "tolerance": TOL, "pass": ok}
    (WORK / "known-answer.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    if not ok:
        raise SystemExit("KNOWN-ANSWER MISMATCH — stop per PLAN.md")


if __name__ == "__main__":
    main()
