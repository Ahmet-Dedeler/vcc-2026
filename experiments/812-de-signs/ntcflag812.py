"""812 H-NTCB footprint on A/B/C: how many genes would the guide-Huber control-bias rule flag (|z| >= 2 / 3) in each
challenge context, and where do they sit (control CPM rank: inside the 2,200-gene block or not)? Same estimator as
ntcbias812.py (H1: 5-6 flagged genes, M0 accuracy 0.66-0.69 vs CLEAN4 0.43-0.44, +0.00017 overall). Truth-free.

usage: uv run python experiments/812-de-signs/ntcflag812.py   -> work/ntcflag812.json (read-only on data/challenge)
"""
from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def wmedian(m, w):
    o = np.argsort(m, 0)
    ms = np.take_along_axis(m, o, 0)
    ws = w[o]
    cw = np.cumsum(ws, 0) / ws.sum(0)
    k = (cw >= 0.5).argmax(0)
    return ms[k, np.arange(m.shape[1])]


def huber(m, se, n, k=2.0, it=30):
    c = wmedian(m, n.astype(float))
    for _ in range(it):
        r = (m - c[None, :]) / np.maximum(se, 1e-12)
        w = np.where(np.abs(r) <= k, 1.0, k / np.maximum(np.abs(r), 1e-12)) * n[:, None]
        c = (w * m).sum(0) / w.sum(0)
    return c


out = {}
for ctx in "ABC":
    with h5py.File(ROOT / f"data/challenge/context_{ctx}.h5ad", "r") as f:
        X = sp.csr_matrix((f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]), shape=tuple(f["X"].attrs["shape"]))
        ntc = f["obs/ntc_id"]
        codes, cats = ntc["codes"][:], ntc["categories"][:].astype(str)
        vi = f["var"][f["var"].attrs["_index"]]
        gnames = (vi["values"][:] if isinstance(vi, h5py.Group) and "values" in vi else
                  vi["categories"][:][vi["codes"][:]] if isinstance(vi, h5py.Group) else vi[:]).astype(str)
    lib = np.asarray(X.sum(1)).ravel()
    C = (sp.diags((1e6 / lib).astype(np.float32)) @ X).astype(np.float32).tocsc()
    del X
    mean = np.asarray(C.mean(0)).ravel()
    keep = np.flatnonzero(mean > 5)
    C = C[:, keep]
    M = mean[keep].astype(np.float64)
    k = len(cats)
    ind = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(k, len(codes)))
    n = np.asarray(ind.sum(1)).ravel()
    gm = (ind @ C).toarray() / n[:, None]
    gsq = (ind @ C.multiply(C)).toarray() / n[:, None]
    v = np.maximum(gsq - gm ** 2, 0) * n[:, None] / (n[:, None] - 1)
    se = np.sqrt(v / n[:, None])
    cellvar = (v * (n[:, None] - 1)).sum(0) / (n.sum() - 1) + ((gm - M) ** 2 * n[:, None]).sum(0) / (n.sum() - 1)
    z = (M - huber(gm, se, n)) / np.sqrt(cellvar / n.sum())
    rank = np.empty(M.size, int)
    rank[np.argsort(-M, kind="stable")] = np.arange(M.size)
    rows = []
    for j in np.flatnonzero(np.abs(z) >= 2):
        dev = (gm[:, j] - M[j]) / se[:, j]
        g = int(np.argmax(np.abs(dev)))
        rows.append(dict(gene=str(gnames[keep[j]]), z=round(float(z[j]), 2), cpm=round(float(M[j]), 1), cpm_rank=int(rank[j]),
                         worst_guide=str(cats[g]), worst_guide_dev_se=round(float(dev[g]), 1),
                         worst_guide_fold=round(float(gm[g, j] / M[j]), 3)))
    out[ctx] = dict(genes=int(M.size), flagged_z2=int((np.abs(z) >= 2).sum()), flagged_z3=int((np.abs(z) >= 3).sum()),
                    flagged_z3_in_top2200=int(((np.abs(z) >= 3) & (rank < 2600)).sum()), rows=sorted(rows, key=lambda r: -abs(r["z"])))
    print(ctx, {kk: vv for kk, vv in out[ctx].items() if kk != "rows"}, out[ctx]["rows"][:8], flush=True)
json.dump(out, open(HERE / "work/ntcflag812.json", "w"), indent=1)
