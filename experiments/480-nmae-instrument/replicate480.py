"""480 §4d: calibrate each surface's gate composition against A/B/C with the scorer's own nmae replicate anchor.

The vcc2026 anchor r for de_wilcoxon_lfc_nmae (cell-eval2 lfc_nmae_ref.py, "raw" estimator) is, per target,
    nmae_ref_raw(p) = mean_{g in gate_full} |lfc_A - lfc_B| / mean_{g in gate_full} |lfc_full|
with disjoint halves of floor(n/2) cells per group (controls split too; each half against its own controls), a gene
that a half does not report (its own > 5 CPM control filter) or reports non-finite read as 0, targets with >= 10 gate
genes, mean over targets, averaged over 5 seeded splits. On A/B/C it is 0.369-0.431 (metrics BRIEF §8). A surface
whose r sits in that band has A/B/C's nmae difficulty (effect size relative to noise on the gate).

Streams the surface's CSR h5 once; per split, CPM sums per (group, half) come from one sparse product per row block.

    uv run python replicate480.py <name> <h5ad-or-h5> <gate-name|323> [--ctrl non-targeting]
"""
from __future__ import annotations

import argparse
import json
import time

import h5py
import numpy as np
import scipy.sparse as sp

from common480 import HERE, gate_from, load_h1
from gates480 import read_labels, read_str

N_SPLITS = 5


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("h5")
    ap.add_argument("gate")
    ap.add_argument("--ctrl", default="non-targeting")
    a = ap.parse_args()
    t0 = time.time()
    if a.gate == "323":
        h1 = load_h1()
        targets, genes, T = h1["targets"], h1["genes"], h1["T"]
    else:
        z = np.load(HERE / f"work/gate-{a.gate}.npz", allow_pickle=True)
        targets, genes = z["targets"].astype(str), z["genes"].astype(str)
        T = gate_from(z["q"].astype(float), z["lfc"].astype(float), genes, targets)
    h = h5py.File(a.h5, "r")
    lab = read_labels(h)
    fg = read_str(h["var"]["_index"])
    col = {g: i for i, g in enumerate(fg)}
    gcols = np.array([col[g] for g in genes])
    X = h["X"]
    shape = tuple(int(v) for v in X.attrs["shape"])
    indptr = X["indptr"][:].astype(np.int64)
    groups = [a.ctrl] + list(targets)
    gi = {g: i for i, g in enumerate(groups)}
    ng = len(groups)
    # half assignment per split: 0 = A, 1 = B, -1 = dropped (odd leftover)
    half = np.full((N_SPLITS, shape[0]), -1, dtype=np.int8)
    grp = np.array([gi.get(x, -1) for x in lab])
    for s in range(N_SPLITS):
        rng = np.random.default_rng(s)
        for k in range(ng):
            idx = np.flatnonzero(grp == k)
            perm = rng.permutation(idx)
            m = perm.size // 2
            half[s, perm[:m]] = 0
            half[s, perm[m:2 * m]] = 1
    sums = np.zeros((N_SPLITS, 2, ng, genes.size))
    ncell = np.zeros((N_SPLITS, 2, ng))
    BLK = 2048
    for lo in range(0, shape[0], BLK):
        hi = min(lo + BLK, shape[0])
        a0, a1 = int(indptr[lo]), int(indptr[hi])
        B = sp.csr_matrix((X["data"][a0:a1].astype(np.float64), X["indices"][a0:a1], indptr[lo:hi + 1] - a0),
                          shape=(hi - lo, shape[1]))
        lib = np.asarray(B.sum(1)).ravel()
        lib[lib == 0] = 1.0
        C = sp.diags(1e6 / lib) @ B[:, gcols]
        g = grp[lo:hi]
        for s in range(N_SPLITS):
            for hh in (0, 1):
                sel = (half[s, lo:hi] == hh) & (g >= 0)
                if not sel.any():
                    continue
                Ind = sp.csr_matrix((np.ones(sel.sum()), (g[sel], np.flatnonzero(sel))), shape=(ng, hi - lo))
                sums[s, hh] += (Ind @ C).toarray()
                ncell[s, hh] += np.bincount(g[sel], minlength=ng)
    h.close()
    G = np.isfinite(T)
    scored = G.sum(1) >= 10
    den = np.where(G, np.abs(np.nan_to_num(T)), 0).sum(1) / np.maximum(G.sum(1), 1)
    per_split = []
    for s in range(N_SPLITS):
        L = []
        for hh in (0, 1):
            mean = sums[s, hh] / np.maximum(ncell[s, hh], 1)[:, None]
            ctrl = mean[0]
            with np.errstate(divide="ignore", invalid="ignore"):
                lfc = np.log2((mean[1:] + 1e-9) / (ctrl + 1e-9))
            lfc[:, ctrl <= 5.0] = 0.0                           # the half's own > 5 CPM control filter
            lfc[~np.isfinite(lfc)] = 0.0
            L.append(lfc)
        num = np.where(G, np.abs(L[0] - L[1]), 0).sum(1) / np.maximum(G.sum(1), 1)
        r = num[scored] / den[scored]
        per_split.append(float(np.mean(r)))
    out = dict(name=a.name, r_raw=float(np.mean(per_split)), per_split=per_split, n_scored=int(scored.sum()),
               n_targets=len(targets), secs=time.time() - t0)
    print(json.dumps(out))
    rec = HERE / "work/replicate480.json"
    allr = json.load(open(rec)) if rec.exists() else {}
    allr[a.name] = out
    json.dump(allr, open(rec, "w"), indent=1)


if __name__ == "__main__":
    main()
