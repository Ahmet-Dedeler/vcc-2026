import sys, numpy as np, scipy.sparse as sp
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, emit_c = E.load_reference()
z = np.load(HERE / "work/explore/block_s0_01.npz")
X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
n_car = int(z["n_car"])
s_, j_, p, q, l = E.block_stats(X, cc)
# instrument _offset
orig = rp423._offset
log = []
def wrapped(x, s, c, R, room, cw, sgn, T=4):
    res = orig(x, s, c, R, room, cw, sgn, T)
    lt0, eq0 = rp423._pos(x, s, c)
    usable = (eq0 == 0)
    g1 = 2 * (rp423._pos(x + 1, s, c)[0] - lt0)
    l1 = 2 * (lt0 - rp423._pos(np.maximum(x - 1, 0), s, c)[0])
    log.append((R, sgn, res is not None, int(usable.sum()), int(((g1 > 0) & usable).sum()), int(((l1 > 0) & usable & (x >= 1)).sum()),
                int(np.sort(g1[usable])[-5:].sum()) if usable.any() else 0))
    return res
rp423._offset = wrapped
X2, dg = rp423.postpass(X, n_car, cc, l, rp423.parse_variant("all_k140_bal"))
ok = [t for t in log if t[2]]; bad = [t for t in log if not t[2]]
print(len(ok), len(bad))
for t in bad[:25]: print("FAIL R=%d sgn=%d usable=%d n_gain=%d n_loss=%d top5gain=%d" % (t[0], t[1], t[3], t[4], t[5], t[6]))
for t in ok[:10]: print("OK R=%d sgn=%d usable=%d n_gain=%d n_loss=%d" % (t[0], t[1], t[3], t[4], t[5]))
