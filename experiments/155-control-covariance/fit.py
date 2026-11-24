"""155 — fit the control-covariance direction estimator.

HYPOTHESIS.  In the 18,400 non-targeting control cells of a context, a target gene T's
expression varies naturally.  Cells low in T are a weak natural knockdown of T.  Regressing
every gene on T across control cells, after residualising the dominant confounds (library
size and global cell state), estimates a CONTEXT-SPECIFIC downstream direction for T that
needs no cross-line transfer.

    knockdown of T moves T down, so the predicted effect is  -beta_{T,g}

CONTROLS ONLY.  This script reads exclusively the cells labelled `non-targeting` in
orion-hct116-thinned.h5ad.  The 132 perturbation cohorts are never touched here; they are
opened only by the scorer, through 095's cached reference DE table.

Covariates: intercept, log library size, and the top-k PCs of the log-normalised control
matrix (k swept).  The whole estimator is one residualisation plus one cross-product:
    R = X - Z (Z^+ X);   beta[t,g] = <R_t, R_g> / <R_t, R_t>

Writes work/coef_k{K}.npz with `beta`, `corr`, `genes`, `targets`, plus diagnostics.
"""
from __future__ import annotations

import argparse
import os
import sys

import anndata as ad
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
HARNESS = os.path.join(REPO, "data/harness/orion-hct116-thinned.h5ad")
SURFACE = os.path.join(REPO, "experiments/095-clean-harness/de_thinned.npz")
CTRL, PERT_COL = "non-targeting", "target_gene"


def norm_log(x: np.ndarray) -> np.ndarray:
    """095's transform, verbatim: CPM to 1e6 per cell, then log1p."""
    s = x.sum(1)
    s[s == 0] = 1.0
    return np.log1p(x * (1e6 / s[:, None]))


def build_control_matrix():
    """Stream the control cells out of the harness onto the surface's gated gene axis."""
    z = np.load(SURFACE, allow_pickle=True)
    genes = np.asarray([str(g) for g in z["genes"]])
    targets = np.asarray([str(p) for p in z["perts"]])
    del z

    A = ad.read_h5ad(HARNESS)                          # sparse csr, ~1.4 GB resident
    lab = np.asarray(A.obs[PERT_COL].astype(str).values, dtype="U64")
    ci = np.flatnonzero(lab == CTRL)
    var = np.asarray([str(g) for g in A.var_names])
    gpos = {g: i for i, g in enumerate(var)}
    # 9 of the 132 targets sit BELOW the 5-CPM gate, so they are not columns of the scored
    # axis -- but we still need their own expression as the regressor.  Carry them as an
    # appended block; `n_gated` marks where the scored axis ends.
    extra = np.asarray([t for t in targets if t not in set(genes)])
    axis = np.concatenate([genes, extra])
    gcol = np.array([gpos[g] for g in axis])
    Xs = A.X[ci]                                       # CONTROLS ONLY -- perturbed rows dropped
    del A
    print(f"controls {ci.size}, gated genes {genes.size}, targets {targets.size}", flush=True)

    X = np.empty((ci.size, axis.size), dtype=np.float32)
    lib = np.empty(ci.size, dtype=np.float64)
    CH = 2000
    for a in range(0, ci.size, CH):
        m = np.asarray(Xs[a:a + CH].todense(), dtype=np.float64)
        lib[a:a + m.shape[0]] = m.sum(1)
        X[a:a + m.shape[0]] = norm_log(m)[:, gcol].astype(np.float32)
        del m
    del Xs
    return X, lib, genes, targets, axis


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, nargs="+", default=[0, 20, 50],
                    help="numbers of PCs to residualise out (library size is always removed)")
    a = ap.parse_args()

    cache = os.path.join("/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/aa534728-8641-48d7-a125-17e75d9c9399/scratchpad", "ctrl155.npz")
    if os.path.exists(cache):
        d = np.load(cache, allow_pickle=True)
        X, lib = d["X"], d["lib"]
        genes = np.asarray([str(g) for g in d["genes"]])
        targets = np.asarray([str(t) for t in d["targets"]])
        axis = np.asarray([str(g) for g in d["axis"]])
        print(f"cached control matrix {X.shape}", flush=True)
    else:
        X, lib, genes, targets, axis = build_control_matrix()
        np.savez(cache, X=X, lib=lib, genes=genes, targets=targets, axis=axis)

    apos = {g: i for i, g in enumerate(axis)}
    tcol = np.array([apos[t] for t in targets])                 # every target has a column
    n_gated = genes.size
    print(f"axis {axis.size} = {n_gated} scored genes + {axis.size - n_gated} sub-gate targets")

    # Expression level of each target in the controls -- the estimator's power budget.
    cpm = np.expm1(X[:, tcol].astype(np.float64)).mean(0)
    print(f"target control CPM: min {cpm.min():.1f}  q25 {np.percentile(cpm, 25):.1f}  "
          f"median {np.median(cpm):.1f}  max {cpm.max():.1f}")
    frac0 = (X[:, tcol] == 0).mean(0)
    print(f"target zero-fraction: median {np.median(frac0):.3f}  max {frac0.max():.3f}")

    # ---- PCA basis, computed once at the largest k and truncated ------------------------
    kmax = max(a.k)
    pcs = None
    if kmax > 0:
        from sklearn.decomposition import PCA
        p = PCA(n_components=kmax, svd_solver="randomized", random_state=0)
        pcs = p.fit_transform(X).astype(np.float64)
        print(f"PCA {kmax} comps, variance explained {p.explained_variance_ratio_.sum():.3f} "
              f"(PC1 {p.explained_variance_ratio_[0]:.3f})", flush=True)

    loglib = np.log(lib)
    for k in a.k:
        cols = [np.ones(X.shape[0]), (loglib - loglib.mean()) / loglib.std()]
        if k:
            cols += [pcs[:, j] / pcs[:, j].std() for j in range(k)]
        Z = np.column_stack(cols)
        # R = X - Z (Z^+ X), done in float64 through the small (k+2) basis
        Q = np.linalg.qr(Z)[0].astype(np.float32)                # (N, k+2) orthonormal
        R = X.copy()
        R -= Q @ (Q.T @ R)                                       # residualise, in place

        T64 = R[:, tcol].astype(np.float64)                       # (N, 132)
        tss = (T64 ** 2).sum(0)
        gss = np.empty(R.shape[1])
        cross = np.empty((targets.size, R.shape[1]))
        CH = 1000
        for c0 in range(0, R.shape[1], CH):
            B = R[:, c0:c0 + CH].astype(np.float64)
            gss[c0:c0 + B.shape[1]] = (B ** 2).sum(0)
            cross[:, c0:c0 + B.shape[1]] = T64.T @ B
            del B
        beta = cross / tss[:, None]
        corr = cross / np.sqrt(np.outer(tss, gss))
        del R, T64, cross

        self_beta = np.array([beta[i, tcol[i]] for i in range(targets.size)])
        assert np.allclose(self_beta, 1.0, atol=1e-3), "self-beta must be 1 by construction"
        out = os.path.join(HERE, "work", f"coef_k{k}.npz")
        np.savez_compressed(out, beta=beta[:, :n_gated].astype(np.float32),
                            corr=corr[:, :n_gated].astype(np.float32),
                            genes=genes, targets=targets,
                            target_cpm=cpm, target_resid_sd=np.sqrt(tss / X.shape[0]))
        print(f"k={k:3d}  self-beta {self_beta.min():.3f}-{self_beta.max():.3f}  "
              f"|corr| median {np.median(np.abs(corr)):.4f} p99 "
              f"{np.percentile(np.abs(corr), 99):.4f}  -> {out}", flush=True)
        del beta, corr


if __name__ == "__main__":
    main()
