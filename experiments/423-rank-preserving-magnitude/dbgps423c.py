import sys, numpy as np, scipy.sparse as sp
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, _ = E.load_reference()
fails = 0
for i in range(int(sys.argv[1]), int(sys.argv[2])):
    z = np.load(HERE / f"work/explore/block_s0_bps_{i:02d}.npz")
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    n_car = int(z["n_car"])
    s_, j_, p, q, l = E.block_stats(X, cc)
    for v in ("dn_k250_zs", "all_k140_zs", "dn_k200", "dn_k200_zspc"):
        try:
            X2, dg = rp423.postpass(X, n_car, cc, l, rp423.parse_variant(v))
        except AssertionError as e:
            fails += 1; print(i, v, "K", 400 - n_car, "FAIL", e, flush=True)
print("fails", fails)
