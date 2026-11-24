import sys, time, cProfile, pstats, numpy as np, scipy.sparse as sp
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, emit_c = E.load_reference()
z = np.load(HERE / "work/explore/block_s0_01.npz")
X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
n_car = int(z["n_car"])
s_, j_, p, q, l = E.block_stats(X, cc)
for v in sys.argv[1].split(","):
    t0 = time.time()
    pr = cProfile.Profile(); pr.enable()
    X2, dg = rp423.postpass(X, n_car, cc, l, rp423.parse_variant(v))
    pr.disable()
    print(v, f"{time.time()-t0:.1f}s", {k: dg[k] for k in ("n_moved", "frac_massw", "frac_massw_zs", "rolled", "offset_pairs")})
    pstats.Stats(pr).sort_stats("cumulative").print_stats(8)
