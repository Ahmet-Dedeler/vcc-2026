"""Which genes carry the pc nmae gain on a dumped target: carrier structure, sign agreement with the reference."""
import sys, numpy as np, scipy.sparse as sp
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, emit_c = E.load_reference()
genes_k = r["genes"][cc.keep]
for ti in map(int, sys.argv[2].split(",")):
    z = np.load(HERE / f"work/explore/block_s0_{ti:02d}.npz")
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    n_car = int(z["n_car"])
    s_, j_, p, q, l = E.block_stats(X, cc)
    X2, dg = rp423.postpass(X, n_car, cc, l, rp423.parse_variant(sys.argv[1]))
    s2, j2, p2, q2, l2 = E.block_stats(X2, cc)
    ql = r["q"][ti]; ll = r["lfc"][ti]
    S = (ql < .05) & np.isfinite(ll) & (ll != 0) & (genes_k != r["targets"][ti])
    den = np.abs(ll[S]).sum()
    contrib = (np.abs(l2 - ll) - np.abs(l - ll)) / den      # per-gene nmae change (only S genes count)
    contrib = np.where(S, contrib, 0)
    Xk = X[:n_car][:, cc.keep]
    nz = (Xk > 0).sum(0)
    L = X.sum(1); sc = L / 1e6
    above = np.array([(np.log1p(Xk[:, j] / sc[:n_car]) > cc.cs[j][-1]).sum() for j in range(Xk.shape[1])])
    print(f"target {ti} {r['targets'][ti]} |S|={S.sum()} nmae change {contrib.sum():+.4f}; touched&moved in S: {int((S & (l2 != l)).sum())}")
    moved = S & (l2 != l)
    sa = np.sign(l) == np.sign(ll)
    print(f"  moved-in-S sign acc {sa[moved].mean():.3f} (all S {sa[S].mean():.3f}); gain from right-sign {contrib[moved & sa].sum():+.4f}, wrong-sign {contrib[moved & ~sa].sum():+.4f}")
    print(f"  moved-in-S: median nonzero carriers {np.median(nz[moved]):.0f}, median above-max {np.median(above[moved]):.0f}; |l| med {np.median(np.abs(l[moved])):.3f}, |l2| med {np.median(np.abs(l2[moved])):.3f}, |ref| med {np.median(np.abs(ll[moved])):.3f}")
    o = np.argsort(contrib)[:8]
    for j in o:
        print(f"   {genes_k[j]:>10s} ref {ll[j]:+.3f} l {l[j]:+.3f} -> {l2[j]:+.3f} d_nmae {contrib[j]:+.5f} nonzero {nz[j]} above {above[j]} ctrl_cpm {cc.ref_lin[j]:.1f}")
