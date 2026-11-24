"""Quick loop on the dumped base blocks (first N targets): nmae on those targets + exact verification."""
import sys, json, time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
HERE = Path(__file__).resolve().parent
variants = sys.argv[1].split(",")
import run423 as RN
import rp423
E = RN.E
r, cc, emit_c = E.load_reference()
genes_k = r["genes"][cc.keep]
files = sorted((HERE / "work/explore").glob("block_s0_*.npz"))

def nmae_one(i, lfc):
    ql = r["q"][i]; ll = r["lfc"][i]
    S = (ql < .05) & np.isfinite(ll) & (ll != 0) & (genes_k != r["targets"][i])
    if S.sum() < 10:
        return np.nan
    lp = np.where(np.isfinite(lfc), lfc, 0.0)
    return float(np.abs(lp[S] - ll[S]).sum() / np.abs(ll[S]).sum())

res = {v: [] for v in ["base"] + variants}
for f in files:
    z = np.load(f); i = int(f.stem.split("_")[-1])
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    n_car = int(z["n_car"])
    s, j, p, q, l = E.block_stats(X, cc)
    res["base"].append(nmae_one(i, l))
    for v in variants:
        t0 = time.time()
        X2, dg = rp423.postpass(X, n_car, cc, l, rp423.parse_variant(v))
        t1 = time.time()
        s2, j2, p2, q2, l2 = E.block_stats(X2, cc)
        vf = rp423.verify(X, X2, n_car, cc, s, p, q, l, s2, p2, q2, l2)
        res[v].append(nmae_one(i, l2))
        dg.pop("_fr_by_gene"); dg.pop("_touched")
        print(i, v, f"pass {t1-t0:.1f}s", f"nmae {res['base'][-1]:.4f} -> {res[v][-1]:.4f}", "jk", round(j, 3), round(j2, 3),
              json.dumps({k: dg[k] for k in ("n_touch", "n_moved", "counts_moved", "frac_mean", "frac_massw", "frac_massw_up", "frac_massw_dn", "frac_massw_zs", "rolled", "lim_ballast_subgate", "ballast_sub_net", "dU_net_abs_max", "n_dU_ball", "dU_ball_abs_max", "offset_pairs", "sub_used", "sub_budget_exhausted_cells")}), json.dumps(vf), flush=True)
for v, a in res.items():
    a = np.array(a); print(v, "nmae mean over", np.isfinite(a).sum(), "targets:", np.nanmean(a))
