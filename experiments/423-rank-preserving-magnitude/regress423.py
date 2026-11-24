"""Regression: the current rp423 on dumped seed-0 base blocks must reproduce the stored full-run lfc/p for those targets."""
import sys, numpy as np, scipy.sparse as sp
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, _ = E.load_reference()
bad = 0
for v in sys.argv[1].split(","):
    st = np.load(HERE / f"work/h1/stats-{v}_s0.npz")
    for i in map(int, sys.argv[2].split(",")):
        z = np.load(HERE / f"work/explore/block_s0_{i:02d}.npz")
        X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
        s_, j_, p, q, l = E.block_stats(X, cc)
        X2, dg = rp423.postpass(X, int(z["n_car"]), cc, l, rp423.parse_variant(v))
        s2, j2, p2, q2, l2 = E.block_stats(X2, cc)
        ok = np.array_equal(l2, st["lfc"][i], equal_nan=True) and np.array_equal(p2, st["p"][i]) and float(j2) == float(st["jk"][i])
        bad += not ok
        print(v, i, "REPRODUCED" if ok else "DIFFERS", flush=True)
print(bad, "differ")
