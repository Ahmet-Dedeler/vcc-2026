"""Robustness of exact-U to a float32 scorer path: recompute every DE-universe gene's rank sum with per-cell values
computed in float32 (x * 1e6 / lib in float32, log1p in float32, controls the same way) and count genes whose U moved."""
import sys, numpy as np, scipy.sparse as sp
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, _ = E.load_reference()
import h5py
from pathlib import Path as P
# controls in float32 (h1-test controls)
lab, genes, xr = RN.E.read_real()
xc = xr[lab == "non-targeting"]
lib_c = np.asarray(xc.sum(1)).ravel().astype(np.float32)
keep = cc.keep
ctrl32 = []
Xc = xc[:, np.flatnonzero(keep)].toarray().astype(np.float32)
Vc = np.log1p(Xc * (np.float32(1e6) / lib_c[:, None])).T      # (G, n) float32
Vc.sort(axis=1)
del Xc
def U32(block):
    x = block[:, keep].astype(np.float32); lib = block.sum(1).astype(np.float32)
    V = np.log1p(x * (np.float32(1e6) / lib[:, None])).T
    out = np.empty(V.shape[0])
    for g in range(V.shape[0]):
        lt = np.searchsorted(Vc[g], V[g], "left"); rt = np.searchsorted(Vc[g], V[g], "right")
        out[g] = (2 * lt + (rt - lt)).sum()
    return out
for i in map(int, sys.argv[2].split(",")):
    z = np.load(HERE / f"work/explore/block_s0_{i:02d}.npz")
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    s_, j_, p, q, l = E.block_stats(X, cc)
    u0 = U32(X)
    for v in sys.argv[1].split(","):
        X2, dg = rp423.postpass(X, int(z["n_car"]), cc, l, rp423.parse_variant(v))
        u1 = U32(X2); d = u1 - u0
        print(i, v, "float32 path: genes with U moved", int((d != 0).sum()), "of", int(dg["_touched"].sum()), "touched; max |dU| (half-ranks)", float(np.abs(d).max()), "sum |dU|", float(np.abs(d).sum()), flush=True)
