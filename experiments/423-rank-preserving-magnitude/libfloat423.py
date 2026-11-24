"""Counterfactual: the same pass WITHOUT carrier-library compensation (sub-gate counts restored to the original, so
each moved carrier's library floats by its kept-gene change). Measures the rank-sum change on UNTOUCHED genes."""
import sys, numpy as np, scipy.sparse as sp
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, _ = E.load_reference()
sidx = np.flatnonzero(~cc.keep)
for i in map(int, sys.argv[2].split(",")):
    z = np.load(HERE / f"work/explore/block_s0_{i:02d}.npz")
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    n_car = int(z["n_car"])
    s_, j_, p, q, l = E.block_stats(X, cc)
    R0, T0 = rp423.r1_tie(X, cc)
    for v in sys.argv[1].split(","):
        X2, dg = rp423.postpass(X, n_car, cc, l, rp423.parse_variant(v))
        t = dg["_touched"]
        Xf = X2.copy(); Xf[:n_car][:, sidx] = X[:n_car][:, sidx]          # carriers' sub-gate back: libraries float
        Xf[:n_car, :] = Xf[:n_car, :]
        Xf2 = X2.copy(); Xf2[:n_car, sidx] = X[:n_car, sidx]
        dl = (Xf2.sum(1) - X.sum(1))[:n_car]
        R1, T1 = rp423.r1_tie(Xf2, cc)
        d = R1 - R0
        s2, j2, p2, q2, l2 = E.block_stats(Xf2, cc)
        print(i, v, f"carriers with library change {int((dl != 0).sum())}, |dlib|/lib max {np.max(np.abs(dl) / X.sum(1)[:n_car]):.4f}; "
              f"UNTOUCHED genes with U moved {int((d[~t] != 0).sum())} of {int((~t).sum())}, max |dU| {np.abs(d[~t]).max():.1f}, "
              f"sum |dU| {np.abs(d[~t]).sum():.0f}; call flips {int(((q < .05) != (q2 < .05)).sum())}", flush=True)
