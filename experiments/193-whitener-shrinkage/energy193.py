"""193 — where the prediction's energy sits in the donor covariance's spectrum.

Every regulariser in Phase 1 lost, and the losses grew monotonically with the amount of
shrinkage.  This script is the mechanism: it reports, per spectral decile, how much of the
predicted effect matrix's energy lives there BEFORE whitening and AFTER, and how much of the
between-row DISCRIMINATION (the quantity `pds` actually reads) each band carries.

`pds` is a cosine retrieval statistic, so the object that matters is not energy but energy that
DIFFERS between rows.  For each band `b` the script reports

    share_before  = ||P V_b||_F^2 / ||P V||_F^2
    share_after   = ||P V_b f_b||_F^2 / ||P V f||_F^2
    disc_b        = 1 - (mean off-diagonal cosine of the rows restricted to band b)

so a band with large `share_after` and high `disc_b` is a band the whitener is right to
amplify, whatever a covariance-estimation loss says about it.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

from plolo import build_fold                       # noqa: E402
from run_arms import CONFIGS                       # noqa: E402
from screen193 import FOLD_ORDER, sha256_file      # noqa: E402

A15 = CONFIGS["a15"]
EDGES_REL = (0.0, 1e-9, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0)


def main() -> None:
    out = {}
    for held in FOLD_ORDER:
        fold = build_fold(held)
        wcols = np.asarray(fold["wcols"])
        op = Path(fold["op_path"])
        with np.load(op) as z:
            ev = np.asarray(z["evals"], float)
            V = np.asarray(z["evecs"], float)
        pre = np.asarray(fold["pre"], dtype=np.float64) * A15["alpha"]
        truth = np.asarray(fold["truth"], dtype=np.float64)
        Y = pre[:, wcols] @ V                        # coefficients in the eigenbasis
        T = truth[:, wcols] @ V
        f = ev ** (-A15["gamma"] / 2.0)
        r = ev / ev.max()
        del V
        e_before = (Y * Y).sum(0)
        e_after = ((Y * f) ** 2).sum(0)
        tot_b, tot_a = e_before.sum(), e_after.sum()
        rows = []
        print(f"\n[{held}] operator {op.name} sha {sha256_file(op)[:8]}  "
              f"p={len(ev)}  {len(fold['targets'])} targets")
        print(f"  {'band (ev/evmax)':>22s} {'dims':>6s} {'E_before':>9s} {'E_after':>9s} "
              f"{'meanCos':>8s} {'signalfrac':>10s}")
        for lo, hi in zip(EDGES_REL[:-1], EDGES_REL[1:]):
            m = (r >= lo) & (r < hi) if hi < 1.0 else (r >= lo)
            if not m.any():
                continue
            A = Y[:, m]
            An = A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)
            G = An @ An.T
            n = G.shape[0]
            mc = float(G[~np.eye(n, dtype=bool)].mean())
            # how much of this band's predicted energy is aligned with the TRUTH's energy
            tb = T[:, m]
            corr = float(np.abs((A * tb).sum()) /
                         max(np.linalg.norm(A) * np.linalg.norm(tb), 1e-300))
            rec = dict(lo=lo, hi=hi, dims=int(m.sum()),
                       e_before=float(e_before[m].sum() / tot_b),
                       e_after=float(e_after[m].sum() / tot_a),
                       mean_cos=mc, truth_align=corr)
            rows.append(rec)
            print(f"  [{lo:8.0e},{hi:8.0e}) {int(m.sum()):6d} {rec['e_before']:9.5f} "
                  f"{rec['e_after']:9.5f} {mc:8.4f} {corr:10.5f}")
        out[held] = dict(op=str(op.relative_to(ROOT)), op_sha256=sha256_file(op), bands=rows)
        del Y, T
    (HERE / "work/energy.json").write_text(json.dumps(out, indent=1))
    print("\nwrote work/energy.json")


if __name__ == "__main__":
    main()
