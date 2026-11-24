"""615 stage 0, known-answer rows for the blob statistic (no raw Orion data is read here).

FINDINGS H.4 quotes the share of true effect energy sitting in the single target-independent mean
vector: orion-hek293t-thinned 44.1 %, orion-hct116-thinned 26.4 %, h1-test 3.6 % (mean pairwise cosine
+0.423 / +0.252 / +0.024). The script that produced those numbers is not in the repo, so this one
re-derives the statistic and records which definition reproduces them, BEFORE the statistic is
computed on the raw HEK293T batches (blob615.py).

    share = n * ||mean_t d_t||^2 / sum_t ||d_t||^2,   d_t = b_t - b_ctrl,
    b = log1p(TS * P / sum(P)) of the pooled counts (the scorer's pds space, cell-eval2 bulk TS 5e4),
    genes = the file's axis minus the surface's own target genes (pds drops every panel gene).

Variants (reported, the first is the a-priori definition): TS 5e4 / 1e4, all genes vs genes > 5 CPM in
control. Usage: uv run python blobka615.py   ->  work/blobka615.json
"""
from __future__ import annotations

import json
import os
import sys

import h5py
import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
SURF = {
    "orion-hek293t-thinned": ("data/harness/orion-hek293t-thinned.h5ad", 0.441, 0.423),
    "orion-hct116-thinned": ("data/harness/orion-hct116-thinned.h5ad", 0.264, 0.252),
    "h1-test": ("data/harness/h1-test.h5ad", 0.036, 0.024),
}


def pooled_sums(path: str):
    """Per-label summed counts (labels x genes, float64) from a CSR h5ad, plus labels and genes."""
    import anndata as ad
    A = ad.read_h5ad(path, backed="r")
    vals = A.obs["target_gene"].astype(str).to_numpy()
    genes = np.asarray(A.var_names.astype(str))
    A.file.close()
    cats, codes = np.unique(vals, return_inverse=True)
    with h5py.File(path, "r") as h:
        X = h["X"]
        n, g = int(X["indptr"].shape[0] - 1), len(genes)
        M = sp.csr_matrix((X["data"][:].astype(np.float64), X["indices"][:], X["indptr"][:]), shape=(n, g))
    G = sp.csr_matrix((np.ones(n), (codes, np.arange(n))), shape=(len(cats), n))
    S = np.asarray((G @ M).todense())
    ncell = np.bincount(codes, minlength=len(cats))
    return cats, S, ncell, genes


def share(D: np.ndarray):
    mu = D.mean(0)
    s = float(D.shape[0] * (mu @ mu) / (D * D).sum())
    U = D / np.linalg.norm(D, axis=1, keepdims=True)
    C = U @ U.T
    k = D.shape[0]
    mpc = float((C.sum() - np.trace(C)) / (k * (k - 1)))
    return s, mpc


def main() -> None:
    out = {}
    for name, (rel, ref_share, ref_cos) in SURF.items():
        cats, S, ncell, genes = pooled_sums(os.path.join(REPO, rel))
        ctrl = np.array([c in ("non-targeting", "control", "NTC") for c in cats])
        assert ctrl.sum() == 1, (name, cats[:5])
        tg = [c for c in cats[~ctrl]]
        keep = ~np.isin(genes, tg)
        P, c0 = S[~ctrl], S[ctrl][0]
        cpm_ctrl = 1e6 * c0 / c0.sum()
        rows = {}
        for ts in (5e4, 1e4):
            b = np.log1p(ts * P / P.sum(1, keepdims=True))
            bc = np.log1p(ts * c0 / c0.sum())
            D = b - bc
            for gname, gm in (("all", keep), ("cpm5", keep & (cpm_ctrl > 5))):
                s, mpc = share(D[:, gm])
                rows[f"ts{int(ts)}_{gname}"] = dict(share=round(s, 4), mean_pair_cos=round(mpc, 4),
                                                   genes=int(gm.sum()))
        out[name] = dict(targets=len(tg), cells_per_target=int(np.median(ncell[~ctrl])),
                         controls=int(ncell[ctrl][0]), ref_share=ref_share, ref_cos=ref_cos, rows=rows)
        print(name, json.dumps(out[name]), flush=True)
    os.makedirs(os.path.join(HERE, "work"), exist_ok=True)
    with open(os.path.join(HERE, "work", "blobka615.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    sys.exit(main())
