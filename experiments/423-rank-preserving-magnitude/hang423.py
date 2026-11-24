import sys, faulthandler, numpy as np, scipy.sparse as sp
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, emit_c = E.load_reference()
z = np.load(HERE / f"work/explore/block_s0_{int(sys.argv[2]):02d}.npz")
X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
n_car = int(z["n_car"])
s_, j_, p, q, l = E.block_stats(X, cc)
faulthandler.dump_traceback_later(25, exit=True)
X2, dg = rp423.postpass(X, n_car, cc, l, rp423.parse_variant(sys.argv[1]))
print("done", {k: dg[k] for k in ("n_moved", "frac_massw", "frac_massw_zs")})
