"""Quick check of a few variants on selected dumped targets (nmae on those targets + verification)."""
import sys, json, time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
HERE = Path(__file__).resolve().parent
variants = sys.argv[1].split(","); tis = [int(t) for t in sys.argv[2].split(",")]
import run423 as RN
import rp423
E = RN.E
r, cc, emit_c = E.load_reference()
genes_k = r["genes"][cc.keep]
def nmae_one(i, lfc):
    ql = r["q"][i]; ll = r["lfc"][i]
    S = (ql < .05) & np.isfinite(ll) & (ll != 0) & (genes_k != r["targets"][i])
    lp = np.where(np.isfinite(lfc), lfc, 0.0)
    return float(np.abs(lp[S] - ll[S]).sum() / np.abs(ll[S]).sum())
for i in tis:
    z = np.load(HERE / f"work/explore/block_s0_{i:02d}.npz")
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    n_car = int(z["n_car"])
    s, j, p, q, l = E.block_stats(X, cc)
    b = nmae_one(i, l)
    for v in variants:
        t0 = time.time()
        X2, dg = rp423.postpass(X, n_car, cc, l, rp423.parse_variant(v)); t1 = time.time()
        s2, j2, p2, q2, l2 = E.block_stats(X2, cc)
        dg.pop("_fr_by_gene"); tch = dg.pop("_touched")
        vf = rp423.verify(X, X2, n_car, cc, s, p, q, l, s2, p2, q2, l2, touched=tch)
        print(i, v, f"{t1-t0:.1f}s", f"dnmae {nmae_one(i, l2) - b:+.4f}",
              {k: dg[k] for k in ("n_moved", "frac_massw", "frac_massw_dn", "rolled", "drop_dump", "offset_pairs")},
              {k: vf[k] for k in ("n_genes_U_changed", "U_abs_max", "call_flips", "sign_flips", "p_logmax", "n_untouched_U_changed")}, flush=True)
