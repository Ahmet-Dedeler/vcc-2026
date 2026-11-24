"""812 descriptive — is there between-NTC-guide overdispersion in A/B/C's controls beyond cell sampling?
If the 46 NTC guides' 400-cell means of a gene vary more than 400-cell sampling allows (R_g = var_between / mean sampling
var >> 1), that gene's control mean carries a guide/batch term shared by every target's real lfc, and its real sign on
non-DE targets is less predictable. Per-cell arithmetic CPM (the scorer's DE unit), genes > 5 CPM in the context's controls.

usage: uv run python experiments/812-de-signs/ntcvar812.py   -> work/ntcvar812.json (read-only on data/challenge)
"""
from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
out = {}
for ctx in "ABC":
    p = ROOT / f"data/challenge/context_{ctx}.h5ad"
    with h5py.File(p, "r") as f:
        X = sp.csr_matrix((f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]), shape=tuple(f["X"].attrs["shape"]))
        ntc = f["obs/ntc_id"]
        if isinstance(ntc, h5py.Group):
            codes = ntc["codes"][:]
            cats = ntc["categories"][:].astype(str)
        else:
            raw = ntc[:].astype(str)
            cats, codes = np.unique(raw, return_inverse=True)
    lib = np.asarray(X.sum(1)).ravel()
    C = (sp.diags((1e6 / lib).astype(np.float32)) @ X).astype(np.float32)   # per-cell CPM
    del X
    C = C.tocsc()
    mean = np.asarray(C.mean(0)).ravel()
    keep = mean > 5
    C = C[:, keep]
    G = C.shape[1]
    k = len(cats)
    ind = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(k, len(codes)))
    n = np.asarray(ind.sum(1)).ravel()
    gm = (ind @ C).toarray() / n[:, None]               # guide means k x G
    gsq = (ind @ C.multiply(C)).toarray() / n[:, None]
    wvar = gsq - gm ** 2                                  # within-guide cell variance
    samp = (wvar / n[:, None]).mean(0)                    # expected variance of a guide mean
    between = gm.var(0, ddof=1)
    R = between / np.maximum(samp, 1e-30)
    # share of a 400-cell target mean's noise that a shared guide term would add (tau^2 = max(0, between - samp))
    tau2 = np.maximum(between - samp, 0)
    cellvar = (wvar * n[:, None]).sum(0) / n.sum()
    frac = tau2 / (tau2 + cellvar / 400)
    qs = [0.1, 0.25, 0.5, 0.75, 0.9, 0.99]
    out[ctx] = dict(genes=int(G), guides=int(k), cells_per_guide=[int(n.min()), int(n.max())],
                    R_quantiles={str(q): round(float(np.quantile(R, q)), 3) for q in qs},
                    share_R_gt_1_5=round(float((R > 1.5).mean()), 4), share_R_gt_2=round(float((R > 2).mean()), 4),
                    share_R_gt_4=round(float((R > 4).mean()), 4),
                    guide_term_share_of_target_noise_quantiles={str(q): round(float(np.quantile(frac, q)), 4) for q in qs},
                    R_by_cpm={f"{lo}-{hi:g}": round(float(np.median(R[(mean[keep] >= lo) & (mean[keep] < hi)])), 3)
                              for lo, hi in ((5, 20), (20, 100), (100, 1000), (1000, 1e9))})
    print(ctx, json.dumps(out[ctx]), flush=True)
json.dump(out, open(HERE / "work/ntcvar812.json", "w"), indent=1)
