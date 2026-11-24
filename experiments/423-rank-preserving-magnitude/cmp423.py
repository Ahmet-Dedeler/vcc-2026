"""Per-gene nmae contributions of two variants on one dumped target."""
import sys, numpy as np, scipy.sparse as sp
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, emit_c = E.load_reference()
genes_k = r["genes"][cc.keep]
ti = int(sys.argv[3])
z = np.load(HERE / f"work/explore/block_s0_{ti:02d}.npz")
X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
n_car = int(z["n_car"])
s_, j_, p, q, l = E.block_stats(X, cc)
ql = r["q"][ti]; ll = r["lfc"][ti]
S = (ql < .05) & np.isfinite(ll) & (ll != 0) & (genes_k != r["targets"][ti]); den = np.abs(ll[S]).sum()
out = {}
for v in (sys.argv[1], sys.argv[2]):
    X2, dg = rp423.postpass(X, n_car, cc, l, rp423.parse_variant(v))
    l2 = E.block_stats(X2, cc)[4]
    out[v] = (l2, np.where(S, (np.abs(l2 - ll) - np.abs(l - ll)) / den, 0))
a, b = out[sys.argv[1]][1], out[sys.argv[2]][1]
Xk = X[:n_car][:, cc.keep]; sc = X.sum(1)[:n_car] / 1e6
for j in np.argsort(a - b)[:12]:
    x = Xk[:, j]; nz = np.flatnonzero(x)
    above = (np.log1p(x / sc) > cc.cs[j][-1]).sum()
    print(f"{genes_k[j]:>10s} ref {ll[j]:+.2f} l {l[j]:+.3f} {sys.argv[1]} {out[sys.argv[1]][0][j]:+.3f} ({a[j]:+.5f}) {sys.argv[2]} {out[sys.argv[2]][0][j]:+.3f} ({b[j]:+.5f}) nonzero {nz.size} above {above} zeros_ctrl {(cc.cs[j]==0).sum()} ball {X[n_car:, np.flatnonzero(cc.keep)[j]]}")
