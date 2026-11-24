"""293 -- shared helpers.  All response arithmetic and statistics come from 285's `lib285.py` verbatim
(which itself imports 218's `_per_target` verbatim); nothing mathematical is new here.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/293-challenge-construct-ceiling"
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(ROOT / "experiments/285-reliability-benchmark/work"))
import lib285 as L  # noqa: E402

CTRL = "non-targeting"
CHUNK = 2048

SPLITS = {
    "train": ROOT / "data/external/vcc2025/train/adata_Training.h5ad",
    "validation": ROOT / "data/external/vcc2025/validation/adata_Validation.h5ad",
    "test": ROOT / "data/external/vcc2025/test/adata_Test.h5ad",
}
H1TEST = ROOT / "data/harness/h1-test.h5ad"

# frozen by PREREGISTRATION.md sec 3
RUNGS = (25, 50, 100, 200)
NDRAW = 8
POOL_CAP = 600       # cells kept per target for the matched-n centring / distractor pool
CONS_CAP = 1600      # cells kept per two-construct construct


def labels_of(h, path):
    o = h[path]
    if isinstance(o, h5py.Group) and "categories" in o:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in o["categories"][:]])
        return cats[o["codes"][:]]
    v = o["values"][:] if isinstance(o, h5py.Group) else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def construct_of(gid: str) -> str:
    """'GENE_P1_A|GENE_P1_B' -> 'GENE_P1'; 'GENE_P1P2_A|GENE_P1P2_B' -> 'GENE_P1P2' (ONE construct)."""
    parts = gid.split("|")
    stems = [(m.group(1) if (m := re.match(r"^(.*)_([AB])$", p)) else p) for p in parts]
    return stems[0] if len(set(stems)) == 1 else "|".join(parts)


def stream_group_sums(path: Path, row_groups: np.ndarray, n_groups: int, chunk: int = CHUNK):
    """Exact per-group integer count sums (row_groups[i] = group of row i, -1 = skip).  273's `stream_sums`."""
    with h5py.File(path, "r") as h:
        n, g = (int(x) for x in h["X"].attrs["shape"])
        indptr = h["X/indptr"][:]
        acc = np.zeros((n_groups, g), dtype=np.float64)
        for r0 in range(0, n, chunk):
            r1 = min(r0 + chunk, n)
            s, e = int(indptr[r0]), int(indptr[r1])
            Xb = sp.csr_matrix((h["X/data"][s:e].astype(np.float64), h["X/indices"][s:e],
                                indptr[r0:r1 + 1] - indptr[r0]), shape=(r1 - r0, g))
            grp = row_groups[r0:r1]
            keep = grp >= 0
            if keep.any():
                oh = sp.csr_matrix((np.ones(int(keep.sum())), (grp[keep], np.flatnonzero(keep))),
                                   shape=(n_groups, r1 - r0))
                acc += np.asarray((oh @ Xb).todense())
    cnt = np.bincount(row_groups[row_groups >= 0], minlength=n_groups)
    return acc, cnt


def stream_cells(path: Path, want: np.ndarray, chunk: int = CHUNK, log=print):
    """Extract the rows `want` (sorted row indices) as a dense uint16 matrix, streaming CSR row blocks."""
    want = np.sort(np.asarray(want, dtype=np.int64))
    with h5py.File(path, "r") as h:
        n, g = (int(x) for x in h["X"].attrs["shape"])
        indptr = h["X/indptr"][:]
        out = np.zeros((want.size, g), dtype=np.uint16)
        mark = np.zeros(n, dtype=bool)
        mark[want] = True
        pos = np.cumsum(mark) - 1                      # destination row of each source row
        mx = 0
        for r0 in range(0, n, chunk):
            r1 = min(r0 + chunk, n)
            m = mark[r0:r1]
            if not m.any():
                continue
            s, e = int(indptr[r0]), int(indptr[r1])
            Xb = sp.csr_matrix((h["X/data"][s:e].astype(np.float32), h["X/indices"][s:e],
                                indptr[r0:r1 + 1] - indptr[r0]), shape=(r1 - r0, g))
            D = np.asarray(Xb[np.flatnonzero(m)].todense())
            mx = max(mx, float(D.max()) if D.size else 0.0)
            assert mx < 65535, f"count {mx} overflows uint16"
            assert np.array_equal(D, np.rint(D)), "non-integer counts"
            out[pos[r0:r1][m]] = D.astype(np.uint16)
    return out, want, mx
