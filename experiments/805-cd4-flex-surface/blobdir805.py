"""805: is the CD4 truth's blob (split-half share 0.120) technical or biological? (lab lead's request, 10-16 03:10Z)

    uv run python experiments/805-cd4-flex-surface/blobdir805.py <truth.h5ad> <out.json>

In 615's space (pseudobulk log1p(5e4 · CPM), target genes excluded) on the truth's targets:
  mu    = mean effect over targets (the blob direction), pc1 = top PC of the mean-removed effects;
  hek   = the thinned Orion HEK293T truth's mean vector (615: the "long-gene axis", the known technical blob);
  depth = the controls' own depth axis: b(top-quartile-library controls) - b(bottom-quartile controls);
cosines among them, the share of mu's energy on the depth axis, the per-target Spearman of the library ratio
(target median library / control median) with the projection of d_t on the depth axis, and the split-half share with
the depth axis projected out of every effect. Descriptive: no rule attached (the gate's verdict stands).
"""
from __future__ import annotations

import json
import os
import sys

import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp
from scipy.stats import spearmanr

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "experiments/615-orion-hek-harness"))
import blob615  # noqa: E402

TS = 5e4


def b(x):
    x = np.atleast_2d(x)
    return np.log1p(TS * x / x.sum(1, keepdims=True))


def cos(a, c):
    return float(a @ c / (np.linalg.norm(a) * np.linalg.norm(c)))


def main() -> None:
    path, outp = sys.argv[1], sys.argv[2]
    A = ad.read_h5ad(path, backed="r")
    lab = A.obs["target_gene"].astype(str).to_numpy()
    genes = np.asarray(A.var_names.astype(str))
    A.file.close()
    with h5py.File(path, "r") as h:
        X = h["X"]
        M = sp.csr_matrix((X["data"][:].astype(np.float64), X["indices"][:], X["indptr"][:]), shape=(len(lab), len(genes)))
    lib = np.asarray(M.sum(1)).ravel()
    cats, codes = np.unique(lab, return_inverse=True)
    rng = np.random.default_rng(615)
    half = np.zeros(len(lab), dtype=np.int64)
    for k in range(len(cats)):
        idx = np.flatnonzero(codes == k)
        half[rng.permutation(idx)[len(idx) // 2:]] = 1
    grp = codes * 2 + half
    G = sp.csr_matrix((np.ones(len(lab)), (grp, np.arange(len(lab)))), shape=(2 * len(cats), len(lab)))
    S = np.asarray((G @ M).todense())
    c = int(np.flatnonzero(cats == "non-targeting")[0])
    keep = ~np.isin(genes, cats)
    tg = [k for k in range(len(cats)) if k != c]
    SA, SB = S[0::2][tg], S[1::2][tg]
    bc = b(S[2 * c] + S[2 * c + 1])[0]
    D = (b(SA + SB) - bc)[:, keep]
    DA, DB = (b(SA) - bc)[:, keep], (b(SB) - bc)[:, keep]
    mu = D.mean(0)
    Dc = D - mu
    _, sv, vt = np.linalg.svd(Dc, full_matrices=False)
    pc1 = vt[0]
    # the controls' own depth axis
    ci = np.flatnonzero(lab == "non-targeting")
    lo, hi = np.quantile(lib[ci], [0.25, 0.75])
    top = np.asarray(M[ci[lib[ci] >= hi]].sum(0)).ravel()
    bot = np.asarray(M[ci[lib[ci] <= lo]].sum(0)).ravel()
    depth = (b(top) - b(bot))[0][keep]
    # 615's technical blob direction (thinned HEK mean vector) on this gene axis
    _, hek_mu = blob615.old_surface("data/harness/orion-hek293t-thinned.h5ad", None, genes)
    hek = hek_mu[keep]
    u = depth / np.linalg.norm(depth)
    proj = D @ u
    ctrl_med = float(np.median(lib[ci]))
    libratio = np.array([np.median(lib[codes == k]) / ctrl_med for k in tg])

    def strip(Z):
        return Z - np.outer(Z @ u, u)

    n = D.shape[0]
    sh_strip = float(n * (strip(DA).mean(0) @ strip(DB).mean(0)) / np.einsum("ij,ij->", strip(DA), strip(DB)))
    sh = float(n * (DA.mean(0) @ DB.mean(0)) / np.einsum("ij,ij->", DA, DB))
    out = dict(truth=path, n_targets=n, share_sh=sh, share_sh_depth_axis_removed=sh_strip,
               cos_mu_hekblob=cos(mu, hek), cos_pc1_hekblob=abs(cos(pc1, hek)),
               cos_mu_depth=cos(mu, depth), cos_pc1_depth=abs(cos(pc1, depth)), cos_hekblob_depth=cos(hek, depth),
               mu_energy_share_on_depth=float((mu @ u) ** 2 / (mu @ mu)),
               pc1_share_centred=float(sv[0] ** 2 / (sv ** 2).sum()),
               spearman_libratio_vs_depth_projection=float(spearmanr(libratio, proj)[0]),
               libratio_quantiles=np.quantile(libratio, [0.1, 0.5, 0.9]).tolist(),
               control_library_q25_q75=[float(lo), float(hi)])
    json.dump(out, open(outp, "w"), indent=1)
    print(json.dumps(out))


if __name__ == "__main__":
    main()
