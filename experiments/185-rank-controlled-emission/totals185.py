"""185 — cache the SHIPPED per-target, per-gene panel totals of 095's h1-test a15 artifact.

Every arm in this experiment (including the known-answer a15 / a8-flat re-emissions, which go
through production's `emit_counts` exactly as `realize171.py` does) delivers these totals
bit-identically.  Reading them once, streaming, keeps the 3 GB prediction file out of RAM for
the rest of the experiment.

Also records each target's per-cell totals as the artifact shipped them (min/max/sd), which
is the quantity the exact-marginal emitter pins to a constant.
"""
from __future__ import annotations

import os
import sys
import time

import h5py
import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
PRED = os.environ.get(
    "PRED",
    f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/"
    f"g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
)
OUT = os.path.join(HERE, "work", "a15_totals.npz")


def main() -> None:
    t0 = time.time()
    import anndata as ad
    A = ad.read_h5ad(PRED, backed="r")
    genes = np.asarray([str(g) for g in A.var_names])
    lab = A.obs["target_gene"].astype(str).values
    A.file.close()
    cats, codes = np.unique(lab, return_inverse=True)
    with h5py.File(PRED, "r") as f:
        indptr = f["X"]["indptr"][:]
        n_cells, G = codes.size, genes.size
        T = cats.size
        totals = np.zeros((T, G), dtype=np.int64)
        cell_tot = np.zeros(n_cells, dtype=np.int64)
        step = 4000
        for r0 in range(0, n_cells, step):
            r1 = min(n_cells, r0 + step)
            a, b = int(indptr[r0]), int(indptr[r1])
            data = f["X"]["data"][a:b]
            idx = f["X"]["indices"][a:b]
            assert np.all(data == np.round(data)) and np.all(data >= 0)
            rows = np.repeat(np.arange(r0, r1), np.diff(indptr[r0:r1 + 1]))
            di = data.astype(np.int64)
            np.add.at(totals, (codes[rows], idx), di)
            np.add.at(cell_tot, rows, di)
            print(f"  rows {r1}/{n_cells} [{time.time()-t0:.0f}s]", flush=True)
    stats = np.zeros((T, 4))
    for k in range(T):
        s = cell_tot[codes == k]
        stats[k] = (s.min(), s.max(), s.mean(), s.std())
    np.savez_compressed(OUT, targets=cats, genes=genes, totals=totals, cell_stats=stats)
    keep = cats != "non-targeting"
    print(f"wrote {OUT}: {T} groups x {G} genes; {keep.sum()} perturbation targets")
    print(f"panel total S: median {np.median(totals[keep].sum(1)):.0f}, "
          f"per-cell {np.median(totals[keep].sum(1))/400:.1f}")
    print(f"per-cell total sd across cells, median over targets: "
          f"{np.median(stats[keep, 3]):.1f}   (min {stats[keep,0].min():.0f}, "
          f"max {stats[keep,1].max():.0f})")


if __name__ == "__main__":
    main()
