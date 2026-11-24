"""228 — the scorer's DE universe and CPM ranking for each challenge context, cached.

`rc_emit.CallSet` ranks over `control.keep` (mean per-cell CPM > 5 over that context's own
NTC control cells) and orders its second block by `control.ref_lin` (the per-gene mean of
per-cell CPM, which is also the scorer's lfc denominator).  Both are properties of the
CONTROL cells alone, so both are invariant to any arm -- which is exactly why the
`n_cpm = 2,200` half of the shipped call set cannot move and only the `n_head = 400` half can.

    uv run python ctrl228.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
SCORER_MIN_CPM = 5.0
CHUNK = 2000


def control_block(context: str) -> tuple[np.ndarray, np.ndarray]:
    """(keep, ref_lin) for one challenge context, streamed off the CSR."""
    cache = WORK / f"ctrl-{context}.npz"
    if cache.exists():
        z = np.load(cache)
        return np.asarray(z["keep"], bool), np.asarray(z["ref_lin"], np.float64)
    path = ROOT / f"data/challenge/context_{context}.h5ad"
    with h5py.File(path, "r") as h:
        indptr = np.asarray(h["X/indptr"][:])
        n_cells = len(indptr) - 1
        n_genes = int(h["X"].attrs["shape"][1])
        acc = np.zeros(n_genes, dtype=np.float64)
        for s in range(0, n_cells, CHUNK):
            e = min(s + CHUNK, n_cells)
            lo, hi = int(indptr[s]), int(indptr[e])
            idx = np.asarray(h["X/indices"][lo:hi])
            dat = np.asarray(h["X/data"][lo:hi], dtype=np.float64)
            rows = np.repeat(np.arange(s, e), np.diff(indptr[s:e + 1]))
            tot = np.bincount(rows - s, weights=dat, minlength=e - s)
            acc += np.bincount(idx, weights=1e6 * dat / tot[rows - s], minlength=n_genes)
    ref_lin = acc / n_cells
    keep = ref_lin > SCORER_MIN_CPM
    np.savez_compressed(cache, keep=keep, ref_lin=ref_lin)
    print(f"[{context}] {n_cells} control cells, {int(keep.sum())} genes in the "
          f"{SCORER_MIN_CPM}-CPM scorer DE universe", flush=True)
    return keep, ref_lin


if __name__ == "__main__":
    for c in ("A", "B", "C"):
        control_block(c)
