"""169 — cache h1-test's MEAN-PER-CELL control CPM, streamed.

`realize169.py` (following `realize164.py`, following the scorer) computes both the reference
and the predicted `lfc` against `ref_lin`, the mean over control cells of that cell's CPM.
Any deterministic predicted `lfc` on this surface must use the same denominator or it is not
comparable to the table's `real_lfc`.  Streamed with h5py so nothing large is held.
"""
from __future__ import annotations

import os

import h5py
import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
REAL = f"{ROOT}/data/harness/h1-test.h5ad"
COL, CTRL = "target_gene", "non-targeting"


def _cat(h, key):
    g = h[f"obs/{key}"]
    if isinstance(g, h5py.Group):
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c)
                           for c in g["categories"][:]])
        return cats[g["codes"][:]]
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in g[:]])


def main() -> None:
    out = os.path.join(WORK, "h1test-cellcpm.npy")
    if os.path.exists(out):
        print(f"{out} exists")
        return
    with h5py.File(REAL, "r") as h:
        G = int(h["X"].attrs["shape"][1])
        lab = _cat(h, COL)
        isc = lab == CTRL
        indptr = h["X/indptr"][:]
        acc = np.zeros(G)
        n = 0
        step = 4000
        for i in range(0, lab.size, step):
            j = min(i + step, lab.size)
            ip = indptr[i:j + 1]
            d = h["X/data"][ip[0]:ip[-1]].astype(np.float64)
            ix = h["X/indices"][ip[0]:ip[-1]].astype(np.int64)
            nn = np.diff(ip)
            rowid = np.repeat(np.arange(i, j), nn)
            lib = np.bincount(rowid - i, weights=d, minlength=j - i)
            lib[lib == 0] = 1.0
            m = isc[rowid]
            np.add.at(acc, ix[m], d[m] * 1e6 / lib[rowid[m] - i])
            n += int(isc[i:j].sum())
        acc /= n
    np.save(out, acc)
    print(f"wrote {out}: {n} control cells, {(acc > 5.0).sum()} genes over 5 CPM")


if __name__ == "__main__":
    main()
