import sys, numpy as np, scipy.sparse as sp
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, emit_c = E.load_reference()
orig = rp423._offset
log = []
def wrapped(x, s, c, R, room, cw, sgn, T=64):
    res = orig(x, s, c, R, room, cw, sgn, T)
    log.append((R, res is not None))
    return res
rp423._offset = wrapped
for ti in (1, 2):
    z = np.load(HERE / f"work/explore/block_s0_{ti:02d}.npz")
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    s_, j_, p, q, l = E.block_stats(X, cc)
    X2, dg = rp423.postpass(X, int(z["n_car"]), cc, l, rp423.parse_variant("all_k140_ex"))
R = np.array([t[0] for t in log]); ok = np.array([t[1] for t in log])
print("offsets", len(R), "ok", ok.sum(), "failed |R| quantiles", np.percentile(np.abs(R[~ok]), [10, 25, 50, 75, 90]) if (~ok).any() else None)
print("ok |R| quantiles", np.percentile(np.abs(R[ok]), [10, 50, 90]))
