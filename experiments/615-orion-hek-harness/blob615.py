"""615 stage 0: the technical-blob test on Orion RAW HEK293T (PREREGISTRATION-blob.md, b91a9b59).

Reads work/hek-raw.h5 (extract615.py) and, for the known-answer reference values, the three old surfaces.
Statistic: share = n * ||mean d||^2 / sum ||d_t||^2 in the scorer's pds space (TS 5e4, all genes minus the 300
panel genes); the DECISION statistic is its split-half form n * <mu_A, mu_B> / sum <d_tA, d_tB>.

Control configurations: P (all controls), P' (one pool, batch mixture of the target set's perturbed cells),
BM (each target's own batch mixture, per half), G (the old greedy builder's selection; descriptive).
Decision (harness subset): PASS iff share_sh(P') <= 0.10 and share_sh(BM) <= 0.10.

    uv run python blob615.py   ->  work/blob615.json, work/targets615.txt
"""
from __future__ import annotations

import json
import os
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
TS = 5e4
CAP, MIN_CELLS, N_HARNESS, SEED, THRESH = 250, 100, 100, 615, 0.10


def b(x: np.ndarray) -> np.ndarray:
    """pds-space profile of summed counts (rows = groups)."""
    x = np.atleast_2d(x)
    return np.log1p(TS * x / x.sum(1, keepdims=True))


def shares(D: np.ndarray, DA: np.ndarray | None = None, DB: np.ndarray | None = None) -> dict:
    n = D.shape[0]
    mu = D.mean(0)
    out = dict(n=n, share=float(n * (mu @ mu) / (D * D).sum()))
    U = D / np.linalg.norm(D, axis=1, keepdims=True)
    C = U @ U.T
    out["mean_pair_cos"] = float((C.sum() - np.trace(C)) / (n * (n - 1)))
    if DA is not None:
        out["share_sh"] = float(n * (DA.mean(0) @ DB.mean(0)) / np.einsum("ij,ij->", DA, DB))
    Dc = D - mu
    s = np.linalg.svd(Dc, compute_uv=False)
    out["top_pc_share_centred"] = float(s[0] ** 2 / (s ** 2).sum())
    return out


def top_pc(D: np.ndarray) -> np.ndarray:
    Dc = D - D.mean(0)
    _, _, vt = np.linalg.svd(Dc, full_matrices=False)
    return vt[0]


def cos(a: np.ndarray, c: np.ndarray) -> float:
    return float(a @ c / (np.linalg.norm(a) * np.linalg.norm(c)))


def old_surface(rel: str, panel_drop: np.ndarray | None, genes_ref: np.ndarray):
    """Plain and split-half shares of an old h5ad surface (labels split at random per label, seed 615)."""
    import anndata as ad
    A = ad.read_h5ad(os.path.join(REPO, rel), backed="r")
    lab = A.obs["target_gene"].astype(str).to_numpy()
    genes = np.asarray(A.var_names.astype(str))
    A.file.close()
    cats, codes = np.unique(lab, return_inverse=True)
    rng = np.random.default_rng(SEED)
    half = np.zeros(len(lab), dtype=np.int64)
    for k in range(len(cats)):
        idx = np.flatnonzero(codes == k)
        half[rng.permutation(idx)[len(idx) // 2:]] = 1
    with h5py.File(os.path.join(REPO, rel), "r") as h:
        X = h["X"]
        M = sp.csr_matrix((X["data"][:].astype(np.float64), X["indices"][:], X["indptr"][:]),
                          shape=(len(lab), len(genes)))
    grp = codes * 2 + half
    Gm = sp.csr_matrix((np.ones(len(lab)), (grp, np.arange(len(lab)))), shape=(2 * len(cats), len(lab)))
    S = np.asarray((Gm @ M).todense())
    ctrl = np.flatnonzero(np.isin(cats, ["non-targeting", "Non-Targeting", "control"]))
    assert len(ctrl) == 1
    c = int(ctrl[0])
    keep = ~np.isin(genes, cats)
    tg = [k for k in range(len(cats)) if k != c]
    SA, SB = S[0::2][tg], S[1::2][tg]
    C = S[2 * c] + S[2 * c + 1]
    bc = b(C)[0]
    D = (b(SA + SB) - bc)[:, keep]
    DA, DB = (b(SA) - bc)[:, keep], (b(SB) - bc)[:, keep]
    r = shares(D, DA, DB)
    # the surface's mean vector on the reference gene axis (for the long-gene-axis cosine)
    mu = pd.Series(D.mean(0), index=genes[keep]).reindex(genes_ref).fillna(0.0).to_numpy()
    return r, mu


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    panel = pd.read_csv(os.path.join(REPO, "data/challenge/pert_counts.csv"))["target_gene"].astype(str).to_numpy()
    with h5py.File(os.path.join(WORK, "hek-raw.h5"), "r") as h:
        lab = h["obs/target_gene"][:].astype(str)
        bat = h["obs/batch"][:].astype(np.int64)
        umi = h["obs/umi"][:].astype(np.float64)
        kind = h["obs/kind"][:]
        genes = h["var/gene_name"][:].astype(str)
        bnum = h["batches/number"][:]
        Cb = h["batches/ctrl_sums"][:]
        nC = h["batches/n_ctrl"][:].astype(np.float64)
        indptr = h["X/indptr"][:]
        nrow = len(lab)
    G = len(genes)
    keep = ~np.isin(genes, panel)
    bpos = {int(x): i for i, x in enumerate(bnum)}
    bi = np.array([bpos[int(x)] for x in bat])
    nB = len(bnum)
    Cbar = Cb / nC[:, None]                                    # per-batch mean control count vector
    umi_ctrl_b = Cb.sum(1) / nC                                # per-batch mean control UMI (official axis)

    pan = np.flatnonzero(kind == 1)
    vc = pd.Series(lab[pan]).value_counts()
    elig = sorted(vc.index[vc >= MIN_CELLS])
    harness = sorted(np.random.default_rng(SEED).choice(elig, size=N_HARNESS, replace=False).tolist())
    with open(os.path.join(WORK, "targets615.txt"), "w") as f:
        f.write("\n".join(harness) + "\n")
    pidx = {g: i for i, g in enumerate(panel)}

    # Groups: per eligible target, 2 random halves of a seeded <= 250-cell subset (rng per target = [615, panel
    # index], so a target's cells are the same in every set); G: the first 150 cells in file (= batch) order, halved.
    gid_h = np.full(nrow, -1, dtype=np.int64)
    gid_g = np.full(nrow, -1, dtype=np.int64)
    for k, t in enumerate(elig):
        rows = pan[lab[pan] == t]
        rng = np.random.default_rng([SEED, pidx[t]])
        sel = np.sort(rng.choice(rows, size=min(CAP, len(rows)), replace=False))
        perm = rng.permutation(sel)
        gid_h[perm[: len(perm) // 2]] = 2 * k
        gid_h[perm[len(perm) // 2:]] = 2 * k + 1
        first = rows[:150]
        pg = rng.permutation(first)
        gid_g[pg[: len(pg) // 2]] = 2 * k
        gid_g[pg[len(pg) // 2:]] = 2 * k + 1
    K = len(elig)
    # the harness's exact cells (rows of hek-raw.h5, per half), so truth615.py builds the object this test read
    np.savez(os.path.join(WORK, "harness615_rows.npz"), targets=np.asarray(harness),
             **{f"{t}__{hf}": np.flatnonzero(gid_h == 2 * elig.index(t) + hf) for t in harness for hf in (0, 1)})
    S_h = np.zeros((2 * K, G))
    S_g = np.zeros((2 * K, G))
    # batch mixtures (cells per batch) per group
    N_h = np.zeros((2 * K, nB))
    N_g = np.zeros((2 * K, nB))
    m = gid_h >= 0
    np.add.at(N_h, (gid_h[m], bi[m]), 1.0)
    m = gid_g >= 0
    np.add.at(N_g, (gid_g[m], bi[m]), 1.0)
    U_h = np.zeros(2 * K)
    np.add.at(U_h, gid_h[gid_h >= 0], umi[gid_h >= 0])

    need = np.flatnonzero((gid_h >= 0) | (gid_g >= 0))
    step = 4000
    with h5py.File(os.path.join(WORK, "hek-raw.h5"), "r") as h:
        dd, di = h["X/data"], h["X/indices"]
        for s0 in range(0, len(need), step):
            rows = need[s0:s0 + step]
            r0, r1 = int(rows[0]), int(rows[-1]) + 1
            a0, a1 = int(indptr[r0]), int(indptr[r1])
            data = dd[a0:a1].astype(np.float64)
            ind = di[a0:a1]
            ip = indptr[r0:r1 + 1] - a0
            Xb = sp.csr_matrix((data, ind, ip), shape=(r1 - r0, G))
            loc = np.arange(r0, r1)
            for gid, S in ((gid_h, S_h), (gid_g, S_g)):
                g = gid[loc]
                mm = g >= 0
                if mm.any():
                    Gm = sp.csr_matrix((np.ones(int(mm.sum())), (g[mm], np.flatnonzero(mm))),
                                       shape=(2 * K, r1 - r0))
                    S += np.asarray((Gm @ Xb).todense())
            print(f"  rows {s0 + len(rows):,}/{len(need):,}", flush=True)

    def read(set_targets: list[str], S: np.ndarray, N: np.ndarray, tag: str) -> dict:
        ks = [elig.index(t) for t in set_targets]
        A, B = S[[2 * k for k in ks]], S[[2 * k + 1 for k in ks]]
        NA, NB = N[[2 * k for k in ks]], N[[2 * k + 1 for k in ks]]
        full = b(A + B)
        out = {}
        # P: every control in every batch
        bc = b(Cb.sum(0))[0]
        out["P"] = shares((full - bc)[:, keep], (b(A) - bc)[:, keep], (b(B) - bc)[:, keep])
        # P': one pool with the set's perturbed batch mixture
        w = (NA + NB).sum(0)
        bp = b(w @ Cbar)[0]
        out["Pprime"] = shares((full - bp)[:, keep], (b(A) - bp)[:, keep], (b(B) - bp)[:, keep])
        out["Pprime"]["umi_ratio"] = float(((A + B).sum() / (NA + NB).sum()) / ((w @ umi_ctrl_b) / w.sum()))
        # BM: each target's (and half's) own batch mixture
        bmA, bmB, bmF = b(NA @ Cbar), b(NB @ Cbar), b((NA + NB) @ Cbar)
        out["BM"] = shares((full - bmF)[:, keep], (b(A) - bmA)[:, keep], (b(B) - bmB)[:, keep])
        ut = (A + B).sum(1) / (NA + NB).sum(1)
        uc = ((NA + NB) @ umi_ctrl_b) / (NA + NB).sum(1)
        out["BM"]["umi_ratio_median"] = float(np.median(ut / uc))
        out["_mu"] = {k: (full - x)[:, keep].mean(0) for k, x in (("P", bc), ("Pprime", bp), ("BM", bmF))}
        out["_pc"] = {"Pprime": top_pc((full - bp)[:, keep]), "BM": top_pc((full - bmF)[:, keep])}
        return out

    res = {"harness": read(harness, S_h, N_h, "h"), "eligible": read(elig, S_h, N_h, "e")}
    # G: first 150 cells per target in batch order vs controls from the first batches until 18,400
    kG = int(np.searchsorted(np.cumsum(nC), 18400) + 1)
    bG = b(Cb[:kG].sum(0))[0]
    ks = [elig.index(t) for t in harness]
    A, B = S_g[[2 * k for k in ks]], S_g[[2 * k + 1 for k in ks]]
    res["harness"]["G"] = shares((b(A + B) - bG)[:, keep], (b(A) - bG)[:, keep], (b(B) - bG)[:, keep])
    res["harness"]["G"].update(ctrl_batches=kG, ctrl_cells=int(nC[:kG].sum()),
                               pert_batch_span_median=float(np.median(
                                   [np.flatnonzero(N_g[2 * k] + N_g[2 * k + 1]).max() + 1 for k in ks])))

    # known-answer surfaces: split-half reference values and the thinned HEK mean vector (the blob direction)
    ref = {}
    g_ref = genes[keep]
    for name, rel in (("orion-hek293t-thinned", "data/harness/orion-hek293t-thinned.h5ad"),
                      ("orion-hct116-thinned", "data/harness/orion-hct116-thinned.h5ad"),
                      ("h1-test", "data/harness/h1-test.h5ad")):
        r, mu = old_surface(rel, None, g_ref)
        ref[name] = r
        if name == "orion-hek293t-thinned":
            mu_blob = mu
        print(name, json.dumps(r), flush=True)

    for setname in ("harness", "eligible"):
        R = res[setname]
        R["blob_axis_cos"] = {k: cos(v, mu_blob) for k, v in R.pop("_mu").items()}
        R["top_pc_blob_axis_cos"] = {k: abs(cos(v, mu_blob)) for k, v in R.pop("_pc").items()}
    H = res["harness"]
    verdict = ("PASS" if H["Pprime"]["share_sh"] <= THRESH and H["BM"]["share_sh"] <= THRESH else
               "STOP: per-target batch composition" if H["BM"]["share_sh"] <= THRESH else
               "STOP: not batches")
    out = dict(prereg="PREREGISTRATION-blob.md (b91a9b59)", threshold=THRESH, verdict=verdict,
               n_eligible=len(elig), harness_targets=len(harness), cells_cap=CAP,
               harness_cells_median=float(np.median((N_h[[2 * elig.index(t) for t in harness]]
                                                     + N_h[[2 * elig.index(t) + 1 for t in harness]]).sum(1))),
               results=res, reference=ref)
    with open(os.path.join(WORK, "blob615.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ("results", "reference")}))
    for setname in ("harness", "eligible"):
        for cfg in ("P", "Pprime", "BM", "G"):
            if cfg in res[setname]:
                r = res[setname][cfg]
                print(f"{setname:9s} {cfg:7s} share_sh {r['share_sh']:.4f}  share {r['share']:.4f}  "
                      f"pair_cos {r['mean_pair_cos']:+.4f}  topPC {r['top_pc_share_centred']:.4f}")
    print("blob-axis cos (harness):", json.dumps(res["harness"]["blob_axis_cos"]),
          "top-PC:", json.dumps(res["harness"]["top_pc_blob_axis_cos"]))
    print("VERDICT:", verdict)


if __name__ == "__main__":
    sys.exit(main())
