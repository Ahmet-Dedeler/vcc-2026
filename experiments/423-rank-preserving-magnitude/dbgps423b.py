import sys, numpy as np, scipy.sparse as sp, traceback
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, _ = E.load_reference()
orig_settle = rp423.settle_ballast
def settle(X, bal, sidx, L_target):
    cur = X[bal].sum(1); sub = X[np.ix_(bal, sidx)].sum(1)
    print("   settle: diff", (cur - L_target).tolist(), "sub-gate left", sub.tolist())
    return orig_settle(X, bal, sidx, L_target)
rp423.settle_ballast = settle
for i in range(3, 10):
    z = np.load(HERE / f"work/explore/block_s0_bps_{i:02d}.npz")
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    n_car = int(z["n_car"])
    s_, j_, p, q, l = E.block_stats(X, cc)
    try:
        X2, dg = rp423.postpass(X, n_car, cc, l, rp423.parse_variant("dn_k200"))
        print(i, "ok")
    except AssertionError as e:
        print(i, "FAIL", e)
