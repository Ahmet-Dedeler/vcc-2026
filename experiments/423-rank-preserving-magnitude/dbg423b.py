import sys, numpy as np, scipy.sparse as sp
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, emit_c = E.load_reference()
genes_k = r["genes"][cc.keep]
ti = int(sys.argv[1]); gname = sys.argv[2]
z = np.load(HERE / f"work/explore/block_s0_{ti:02d}.npz")
X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
n_car = int(z["n_car"])
j = int(np.flatnonzero(genes_k == gname)[0])
gcol = np.flatnonzero(cc.keep)[j]
L = X.sum(1); s = L / 1e6
x = X[:n_car, gcol]; xb = X[n_car:, gcol]
c = cc.cs[j]
nz = np.flatnonzero(x)
print("gene", gname, "ctrl max CPM", np.expm1(c[-1]), "ctrl zeros", (c == 0).sum(), "top ctrl CPMs", np.round(np.expm1(c[-8:]), 1))
print("nonzero carriers:", [(int(i), int(x[i]), round(float(x[i] / s[i]), 1), int(L[i])) for i in nz])
print("ballast count", xb, "CPM", xb / s[n_car:], "pos", np.searchsorted(c, np.log1p(xb / s[n_car:])))
xl, xh, lt, rt, tied = rp423.cell_bounds(x, s[:n_car], c)
print("holders xl", xl[nz], "xh", np.minimum(xh[nz], 10**9))
