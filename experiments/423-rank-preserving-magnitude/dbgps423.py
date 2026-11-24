import sys, numpy as np, scipy.sparse as sp, traceback
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, _ = E.load_reference()
for i in (0, 1, 2):
    z = np.load(HERE / f"work/explore/block_s0_bps_{i:02d}.npz")
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    n_car = int(z["n_car"])
    print(i, "ballast libs", X[n_car:].sum(1), "ballast sub-gate", X[n_car:][:, ~cc.keep].sum(1), "carrier sub-gate/cell", X[:n_car][:, ~cc.keep].sum(1).mean().round(1))
    s_, j_, p, q, l = E.block_stats(X, cc)
    for v in ("dn_k250_zs", "all_k140_zs", "dn_k200", "dn_k200_zspc"):
        try:
            X2, dg = rp423.postpass(X, n_car, cc, l, rp423.parse_variant(v))
            print("  ", v, "ok", {k: dg[k] for k in ("n_moved", "lim_ballast_subgate", "ballast_sub_net", "sub_used")})
        except AssertionError as e:
            print("  ", v, "FAIL", e); traceback.print_exc(limit=3)
