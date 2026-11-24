"""Mass-holder census on the dumped base blocks: DE-universe genes whose nonzero carriers ALL sit above the control max."""
import numpy as np, scipy.sparse as sp
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN
E = RN.E
r, cc, _ = E.load_reference()
genes_k = r["genes"][cc.keep]; kcol = np.flatnonzero(cc.keep)
cmax = cc.cs[:, -1]
rows = []
for f in sorted((HERE / "work/explore").glob("block_s0_*.npz")):
    z = np.load(f); i = int(f.stem.split("_")[-1])
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    n_car = int(z["n_car"]); s = X.sum(1)[:n_car] / 1e6
    Xk = X[:n_car][:, kcol]
    V = np.log1p(Xk / s[:, None])
    nz = (Xk > 0).sum(0); above = (V > cmax[None, :]).sum(0)
    A = (X[:, kcol] / (X.sum(1) / 1e6)[:, None]).mean(0); l = np.log2((A + 1e-9) / (cc.ref_lin + 1e-9))
    holder = (nz > 0) & (above == nz)
    some = (above > 0) & ~holder
    S = (r["q"][i] < .05) & np.isfinite(r["lfc"][i]) & (r["lfc"][i] != 0) & (genes_k != r["targets"][i])
    ok = S & np.isfinite(l)
    acc = lambda m: float((np.sign(l[m & ok]) == np.sign(r["lfc"][i][m & ok])).mean()) if (m & ok).any() else np.nan
    rows.append((i, holder.sum(), (holder & (l < 0)).sum(), np.median(nz[holder]) if holder.any() else 0, some.sum(),
                 (holder & S).sum(), acc(holder), acc(some), acc(~holder & ~some), S.sum(),
                 float(np.median(np.abs(l[holder]))) if holder.any() else 0, float(np.median(np.abs(l[~holder & np.isfinite(l)])))))
    print("target %d: pure-holder genes %d (down %d, median holders %d), partial-above-max %d | in ref DE: %d; sign acc holder %.3f partial %.3f rest %.3f (|S|=%d); median |l| holder %.3f rest %.3f" % rows[-1])
