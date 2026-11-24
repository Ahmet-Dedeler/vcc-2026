"""190 — cache the control-derived co-expression GRN for h1-test (the D2R2-shaped prior).

C[t, g] = Pearson r(t, g) across the recipient's 18,400 unperturbed cells on log1p-CPM,
with the option of first projecting out the top-k principal components of the cell x gene
matrix.  174 measured that 2 PCs carry the entire raw signal, so pc0 / pc2 / pc5 are cached
and every downstream measurement is reported against all three.

The NETWORK's expected response sign under knockdown of t is  -sign(C[t, g]).

Loads `data/harness/h1-test.h5ad` exactly the way 174 did (same code path, same gene axis
assertion) so the two experiments cannot silently disagree about the matrix.
Peak RSS ~ 3 GB.
"""
from __future__ import annotations

import os
import sys
import time

import h5py
import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.join(ROOT, "experiments", "190-regulon-ordering")
WORK = os.path.join(HERE, "work")
OUT = os.path.join(WORK, "coex190.npz")
sys.path.insert(0, os.path.join(ROOT, "experiments", "153-de-bench"))
sys.path.insert(0, os.path.join(ROOT, "experiments", "168-sharpening"))
PCS = (0, 2, 5)


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    from debench import Surface
    from h1_tiebreak import labels

    t0 = time.time()
    S = Surface(os.path.join(ROOT, "experiments", "164-callset-artefact", "de164.a15.npz"))
    E = np.load(os.path.join(ROOT, "experiments", "169-sharpening-six-member", "work",
                             "h1test-effect.npz"), allow_pickle=True)
    eg = {str(g): i for i, g in enumerate(E["genes"])}
    gi = np.array([eg[g] for g in S.genes])

    with h5py.File(os.path.join(ROOT, "data", "harness", "h1-test.h5ad"), "r") as h:
        genes = labels(h, "var/_index")
        lab = labels(h, "obs/target_gene")
        ci = np.flatnonzero(lab == "non-targeting")
        n, G = (int(x) for x in h["X"].attrs["shape"])
        ip = h["X/indptr"][:]
        rows = []
        for r0 in range(0, ci.size, 3000):
            idx = ci[r0:r0 + 3000]
            blk = np.zeros((idx.size, G), dtype=np.float32)
            for k, r in enumerate(idx):
                s, e = int(ip[r]), int(ip[r + 1])
                blk[k, h["X/indices"][s:e]] = h["X/data"][s:e]
            tot = blk.sum(1, keepdims=True)
            tot[tot == 0] = 1
            rows.append(np.log1p(blk * (1e6 / tot)))
        X = np.concatenate(rows)
        del rows
    assert np.array_equal(genes, np.asarray([str(g) for g in E["genes"]]))
    print(f"controls {X.shape} log-CPM  [{time.time() - t0:.0f}s]", flush=True)

    X -= X.mean(0, keepdims=True)
    tpos = {g: i for i, g in enumerate(genes)}
    tidx = np.array([tpos[str(p)] for p in S.perts])

    from sklearn.utils.extmath import randomized_svd
    store = {"genes": S.genes, "perts": S.perts}
    for k in PCS:
        Z = X
        if k > 0:
            _, _, Vt = randomized_svd(X, n_components=k, n_iter=4, random_state=0)
            Z = X - (X @ Vt.T) @ Vt
        sd = Z.std(0)
        sd = np.where(sd == 0, 1.0, sd)
        Zn = Z / sd
        C = (Zn[:, tidx].T @ Zn) / Zn.shape[0]
        C = C[:, gi].astype(np.float32)
        store[f"C_pc{k}"] = C
        print(f"  pc{k}: |r| median {np.median(np.abs(C)):.4f} "
              f"q99 {np.quantile(np.abs(C), 0.99):.4f}  [{time.time() - t0:.0f}s]", flush=True)
        if k > 0:
            del Zn, Z
    np.savez_compressed(OUT, **store)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
