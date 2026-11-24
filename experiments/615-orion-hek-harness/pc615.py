#!/usr/bin/env python3
"""615 descriptive: is the long-gene axis (the thinned HEK truth's mean vector, FINDINGS H.4's blob direction) also the
top direction of BETWEEN-target variation on h1-test (Arc's chemistry) and on the raw HEK harness?
For each surface: the top-3 PCs of the mean-removed pds-space effects (TS 5e4, panel/target genes dropped), their
energy shares and |cos| with the blob direction, on the genes both axes share. Split-half: each PC is fitted on half A
and its energy share re-measured on half B, so noise does not choose the direction. -> work/pc615.json"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = [sys.argv[0]]
import blob615 as B  # noqa: E402  (read-only: the same profile / reader functions)


def halves(rel):
    import anndata as ad
    import h5py
    import scipy.sparse as sp
    if rel.endswith(".h5ad"):
        A = ad.read_h5ad(os.path.join(B.REPO, rel), backed="r")
        lab = A.obs["target_gene"].astype(str).to_numpy()
        genes = np.asarray(A.var_names.astype(str))
        A.file.close()
    else:                                              # 614/615 CSR-h5 truth
        with h5py.File(os.path.join(B.REPO, rel), "r") as h:
            lab = h["obs/target_gene"][:].astype(str)
            genes = h["var/_index"][:].astype(str)
    panel = set(open(os.path.join(B.REPO, "data/challenge/pert_counts.csv")).read().split()) if not rel.endswith(".h5ad") else set()
    cats, codes = np.unique(lab, return_inverse=True)
    rng = np.random.default_rng(B.SEED)
    half = np.zeros(len(lab), dtype=np.int64)
    for k in range(len(cats)):
        idx = np.flatnonzero(codes == k)
        half[rng.permutation(idx)[len(idx) // 2:]] = 1
    with h5py.File(os.path.join(B.REPO, rel), "r") as h:
        X = h["X"]
        M = sp.csr_matrix((X["data"][:].astype(np.float64), X["indices"][:], X["indptr"][:]), shape=(len(lab), len(genes)))
    G = sp.csr_matrix((np.ones(len(lab)), (codes * 2 + half, np.arange(len(lab)))), shape=(2 * len(cats), len(lab)))
    S = np.asarray((G @ M).todense())
    c = int(np.flatnonzero(np.isin(cats, ["non-targeting", "Non-Targeting", "control"]))[0])
    tg = [k for k in range(len(cats)) if k != c]
    bc = B.b(S[2 * c] + S[2 * c + 1])[0]
    keep = ~np.isin(genes, cats) & ~np.isin(genes, list(panel))
    DA, DB = (B.b(S[0::2][tg]) - bc)[:, keep], (B.b(S[1::2][tg]) - bc)[:, keep]
    full = (B.b(S[0::2][tg] + S[1::2][tg]) - bc)[:, keep]
    return genes[keep], DA, DB, full


def pcs(DA, DB, k=3):
    A = DA - DA.mean(0)
    Bc = DB - DB.mean(0)
    _, s, vt = np.linalg.svd(A, full_matrices=False)
    tot = float(np.einsum("ij,ij->", A, Bc))          # split-half total centred energy (noise-free)
    shares = [float(((A @ v) * (Bc @ v)).sum() / tot) for v in vt[:k]]
    return vt[:k], shares


def main():
    g_thin, _, _, full_thin = halves("data/harness/orion-hek293t-thinned.h5ad")
    blob = dict(zip(g_thin, full_thin.mean(0)))
    out = {}
    for name, rel in (("h1-test", "data/harness/h1-test.h5ad"),
                      ("orion-hek293t-thinned", "data/harness/orion-hek293t-thinned.h5ad"),
                      ("hek615-raw-truth", "experiments/615-orion-hek-harness/work/hek-truth.h5")):
        g, DA, DB, _ = halves(rel)
        bv = np.array([blob.get(x, 0.0) for x in g])
        V, sh = pcs(DA, DB)
        out[name] = {"pc_share_sh": sh, "pc_abs_cos_blob": [float(abs(v @ bv) / np.linalg.norm(bv)) for v in V]}
        print(name, json.dumps(out[name]), flush=True)
    (open(os.path.join(HERE, "work/pc615.json"), "w")).write(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
