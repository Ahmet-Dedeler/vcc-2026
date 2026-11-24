#!/usr/bin/env python3
"""383 — chunked TRUE log1p-mean deltas (never densify a full matrix).

truth_val   per-pert mean(X_log1p) - mean(control X_log1p) on the recipe val template
            (50 perts, 18,080 genes) -> work/truth_val.npz
truth_h1    same on h1-test.h5ad REAL cells (85 perts; raw counts -> log1p here, so the
            STATE-vs-truth cosine uses the same space as STATE predictions)
            -> work/truth_h1.npz
Saves keys: perts (str), genes (str), delta (P x G float32), n (cells per pert).
"""
from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
from build_templates383 import (  # noqa: E402
    H1_TEST,
    VAL_TEMPLATE,
    read_cats,
    read_sparse_block,
    support_gene_list,
)

CHUNK = 4000


def sparse_group_means(h5_path: Path, genes: list[str], log1p_raw: bool) -> dict:
    """Chunked per-target sums/counts over a sparse-count h5ad. Returns dict."""
    with h5py.File(h5_path, "r") as h:
        tg = read_cats(h, "obs/target_gene")
        perts = sorted(set(tg) - {"non-targeting"})
        code_of = {p: i for i, p in enumerate(perts + ["non-targeting"])}
        codes = np.array([code_of[t] for t in tg], dtype=np.int64)
        n_g = len(genes)
        sums = np.zeros((len(perts) + 1, n_g), dtype=np.float64)
        cnts = np.zeros(len(perts) + 1, dtype=np.int64)
        X = h["X"]
        if isinstance(X, h5py.Group):  # csr
            indptr, data, idx = X["indptr"][:], X["data"], X["indices"]
            n_obs = len(tg)
            for a in range(0, n_obs, CHUNK):
                b = min(a + CHUNK, n_obs)
                block = np.zeros((b - a, n_g), dtype=np.float64)
                for i, r in enumerate(range(a, b)):
                    s, e = int(indptr[r]), int(indptr[r + 1])
                    block[i, idx[s:e]] = data[s:e]
                if log1p_raw:
                    np.log1p(block, out=block)
                for k in range(len(perts) + 1):
                    m = codes[a:b] == k
                    if m.any():
                        sums[k] += block[m].sum(0)
                        cnts[k] += int(m.sum())
                if a % 20000 == 0:
                    print("  rows", a, "/", n_obs, flush=True)
        else:  # dense log1p already
            n_obs = X.shape[0]
            for a in range(0, n_obs, CHUNK):
                b = min(a + CHUNK, n_obs)
                block = np.asarray(X[a:b], dtype=np.float64)
                for k in range(len(perts) + 1):
                    m = codes[a:b] == k
                    if m.any():
                        sums[k] += block[m].sum(0)
                        cnts[k] += int(m.sum())
                if a % 20000 == 0:
                    print("  rows", a, "/", n_obs, flush=True)
    means = sums / np.maximum(cnts, 1)[:, None]
    delta = (means[:-1] - means[-1]).astype(np.float32)
    return {"perts": np.array(perts), "delta": delta,
            "n": cnts[:-1].astype(np.int64), "n_ctrl": int(cnts[-1])}


def main(which: str) -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    genes = support_gene_list()
    if which == "truth_val":
        out = sparse_group_means(VAL_TEMPLATE, genes, log1p_raw=False)
    else:
        out = sparse_group_means(H1_TEST, genes, log1p_raw=True)
    out["genes"] = np.array(genes)
    np.savez(WORK / f"{which}.npz", **out)
    print(which, "perts:", len(out["perts"]), "n_ctrl:", out["n_ctrl"],
          "mean|delta|:", round(float(np.abs(out["delta"]).mean()), 4))


if __name__ == "__main__":
    main(sys.argv[1])
