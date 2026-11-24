"""224 — target features for the ladder's rung 1/2/5.

What is actually on disk (checked 2026-09-04):
  * DepMap co-essentiality: `experiments/114-depmap-coessentiality/work/source-features/
    depmap-features.npz`, 17,916 genes x 1,178 columns.  **This is 132's own feature
    artifact** — reusing it is what makes rung 1 a reproduction rather than a re-derivation.
  * STRING network embeddings: `experiments/100-sign-confidence/work/string-pca.npz`,
    21,688 genes x 64 unit-norm dims (from `data/external/slim/gene_string_embeddings.v0.3.h5`).
  * ESM-2 and GenePT: **NOT on disk.**  Neither is downloaded here; 087 already closed the
    STRING target kernel and 123 closed the SLIM reproduction, so the marginal value of a
    third sequence-only embedding on top of DepMap+STRING is small and the rung-1 verdict
    does not turn on it.  Stated, not silently skipped.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

DEPMAP = ROOT / "experiments/114-depmap-coessentiality/work/source-features/depmap-features.npz"
STRING = ROOT / "experiments/100-sign-confidence/work/string-pca.npz"


def target_features(targets: list[str]) -> tuple[np.ndarray, np.ndarray, dict]:
    """(n_targets, d) feature matrix + a bool mask of targets with ANY real feature."""
    blocks, meta = [], {}

    z = np.load(DEPMAP, allow_pickle=True)
    idx = {str(g): i for i, g in enumerate(z["genes"])}
    M = z["matrix"]
    have_d = np.asarray([t in idx for t in targets])
    D = np.zeros((len(targets), M.shape[1]), dtype=np.float64)
    rows = np.asarray([idx[t] for t in targets if t in idx])
    if rows.size:
        D[have_d] = np.nan_to_num(np.asarray(M[rows], dtype=np.float64))
    blocks.append(D)
    meta["depmap_dims"] = int(M.shape[1])
    meta["depmap_coverage"] = int(have_d.sum())

    z = np.load(STRING, allow_pickle=True)
    idx = {str(g): i for i, g in enumerate(z["genes"])}
    U = z["unit"]
    avail = z["available"]
    have_s = np.asarray([t in idx and bool(avail[idx[t]]) for t in targets])
    S = np.zeros((len(targets), U.shape[1]), dtype=np.float64)
    rows = np.asarray([idx[t] for t in targets if t in idx and bool(avail[idx[t]])])
    if rows.size:
        S[have_s] = np.nan_to_num(np.asarray(U[rows], dtype=np.float64))
    blocks.append(S)
    meta["string_dims"] = int(U.shape[1])
    meta["string_coverage"] = int(have_s.sum())

    X = np.hstack(blocks)
    have = have_d | have_s
    meta["any_coverage"] = int(have.sum())
    meta["n_targets"] = len(targets)
    meta["total_dims"] = int(X.shape[1])
    return X, have, meta


def standardise(X: np.ndarray, fit_rows: np.ndarray) -> np.ndarray:
    mu = X[fit_rows].mean(0, keepdims=True)
    sd = X[fit_rows].std(0, keepdims=True)
    return (X - mu) / np.maximum(sd, 1e-8)


def pca(X: np.ndarray, fit_rows: np.ndarray, n: int) -> np.ndarray:
    """PCA fitted on `fit_rows` only, applied to all rows.  132 used PCA(50)."""
    Z = standardise(X, fit_rows)
    A = Z[fit_rows]
    A = A - A.mean(0, keepdims=True)
    _, _, Vt = np.linalg.svd(A, full_matrices=False)
    V = Vt[:min(n, Vt.shape[0])].T
    return (Z - Z[fit_rows].mean(0, keepdims=True)) @ V


if __name__ == "__main__":
    import json
    import sys
    ts = sys.argv[1:] or ["TP53", "MYC", "RPL5", "NOTAREALGENE"]
    X, have, meta = target_features(ts)
    print(json.dumps(meta, indent=1))
    print("X", X.shape, "have", have.sum())
