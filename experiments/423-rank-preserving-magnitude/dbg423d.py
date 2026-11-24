import numpy as np, scipy.sparse as sp
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, emit_c = E.load_reference()
genes_k = r["genes"][cc.keep]
z = np.load(HERE / "work/explore/block_s0_01.npz")
X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
n_car = int(z["n_car"]); j = int(np.flatnonzero(genes_k == "MYO6")[0]); col = np.flatnonzero(cc.keep)[j]
L = X.sum(1); s = L / 1e6; c = cc.cs[j]
x = X[:n_car, col].copy(); xb = X[n_car:, col].copy()
cw = 1.0 / (400 * s[:n_car]); bw = 1.0 / (400 * s[n_car:])
groups, gw = rp423._levels(cw)
xl, xh, *_ = rp423.cell_bounds(x, s[:n_car], c)
room = np.full(n_car, 100)
A = (x / s[:n_car]).sum() / 400 + (xb / s[n_car:]).sum() / 400
target = (cc.ref_lin[j] + 1e-9) * 2 ** (2.5 * np.log2(A / cc.ref_lin[j])) - 1e-9
need = target - A
x2 = rp423._zs_step(x, xl, xh, room, cw, need, groups, gw)
rem = need - ((x2 - x) * cw).sum()
print("A", A, "target", target, "need", need, "after zs rem", rem, "holders", x[x > 0], "->", x2[x2 > 0])
xl2, xh2, *_ = rp423.cell_bounds(x2, s[:n_car], c)
hs = np.flatnonzero(x2 > 0)
for h in hs:
    ys = np.arange(1, xl2[h]); lt0 = np.searchsorted(c, np.log1p(x2[h] / s[h]))
    lty = np.searchsorted(c, np.log1p(ys / s[h])); print("h", h, "L", L[h], "x", x2[h], "xl", xl2[h], "lo(y)", list(zip(ys.tolist(), (2 * (lt0 - lty)).tolist())))
Ds = np.arange(0, 400)
v = np.log1p((xb[0] + Ds) / s[n_car]); db = 2 * np.searchsorted(c, v) + (np.searchsorted(c, v, 'right') - np.searchsorted(c, v)) - (2 * np.searchsorted(c, np.log1p(xb[0] / s[n_car])))
print("db(D) first 60:", db[:60].tolist())
print("db at 100..400 step 20:", db[100:400:20].tolist())
res = rp423._drop_dump(x2, xl2, xh2, room, xb, s[:n_car], s[n_car:], c, cw, bw, rem)
print("drop_dump result", None if res is None else (res[0][res[0] > 0], res[1]))
# instrument
orig_off = rp423._offset
calls = []
def off(x, s_, c_, R, room, cw_, sgn, T=64, frozen=None):
    out = orig_off(x, s_, c_, R, room, cw_, sgn, T, frozen)
    calls.append((R, frozen, out is not None, [(int(i), int(x[i])) for i in np.flatnonzero(x > 0)], [int(room[i]) for i in np.flatnonzero(x > 0)]))
    return out
rp423._offset = off
res = rp423._drop_dump(x2, xl2, xh2, room, xb, s[:n_car], s[n_car:], c, cw, bw, rem)
print(len(calls)); [print(cl) for cl in calls[:15]]
