"""305 -- shared helpers.  All response arithmetic is 285's `lib285` (imported, not re-implemented);
the streaming group-sum pattern is 293's `lib293.stream_group_sums` (273's `stream_sums`), re-typed here
so that nothing in 293's directory is invoked or written.
"""
from __future__ import annotations

import hashlib
import json
import re
import resource
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent          # experiments/305-t7-expression-test/work
sys.path.insert(0, str(ROOT / "experiments/285-reliability-benchmark/work"))
import lib285 as L  # noqa: E402  -- prof(), TS, sha256_file, boot_*, Logger

CTRL = "non-targeting"
CHUNK = 2048                                     # 293's block size
SPLITS = {
    "train": ROOT / "data/external/vcc2025/train/adata_Training.h5ad",
    "validation": ROOT / "data/external/vcc2025/validation/adata_Validation.h5ad",
    "test": ROOT / "data/external/vcc2025/test/adata_Test.h5ad",
}


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def labels_of(h, path):                          # lib293.labels_of, verbatim
    o = h[path]
    if isinstance(o, h5py.Group) and "categories" in o:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in o["categories"][:]])
        return cats[o["codes"][:]]
    v = o["values"][:] if isinstance(o, h5py.Group) else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def construct_of(gid: str) -> str:               # lib293.construct_of, verbatim
    """'GENE_P1_A|GENE_P1_B' -> 'GENE_P1'; a non-targeting pair keeps its full id."""
    parts = gid.split("|")
    stems = [(m.group(1) if (m := re.match(r"^(.*)_([AB])$", p)) else p) for p in parts]
    return stems[0] if len(set(stems)) == 1 else "|".join(parts)


def read_obs(path: Path):
    """obs/guide_id, obs/target_gene, obs/batch, var/_index, X.attrs['shape'] -- no X element read."""
    with h5py.File(path, "r") as h:
        gid = labels_of(h, "obs/guide_id")
        tgt = labels_of(h, "obs/target_gene")
        bat = labels_of(h, "obs/batch")
        genes = labels_of(h, "var/_index")
        shape = tuple(int(v) for v in h["X"].attrs["shape"])
    return gid, tgt, bat, genes, shape


def halves(rows: np.ndarray, batches: np.ndarray) -> np.ndarray:
    """Deterministic, batch-balanced half assignment (0/1) for one construct's cells.

    Within each batch the construct's cells are taken in ascending row order and assigned
    alternately 0,1,0,1,...  No RNG: the assignment is a pure function of the file's row order.
    """
    out = np.empty(rows.size, dtype=np.int8)
    order = np.argsort(rows, kind="stable")
    for b in np.unique(batches):
        m = order[batches[order] == b]
        out[m] = (np.arange(m.size) % 2).astype(np.int8)
    return out


def stream_group_sums(path: Path, groups: np.ndarray, n_groups: int, chunk: int = CHUNK, log=print):
    """Exact per-group integer count sums over ALL genes (groups[i] = group of row i, -1 = skip).

    293's `lib293.stream_group_sums`; the accumulator is (n_groups x n_genes) float64, which holds
    integers exactly up to 2**53.  Returns (sums, cell counts, n_genes).
    """
    t0 = time.time()
    with h5py.File(path, "r") as h:
        n, g = (int(v) for v in h["X"].attrs["shape"])
        indptr = h["X/indptr"][:]
        acc = np.zeros((n_groups, g), dtype=np.float64)
        for r0 in range(0, n, chunk):
            r1 = min(r0 + chunk, n)
            grp = groups[r0:r1]
            keep = grp >= 0
            if not keep.any():
                continue
            s, e = int(indptr[r0]), int(indptr[r1])
            Xb = sp.csr_matrix((h["X/data"][s:e].astype(np.float64), h["X/indices"][s:e],
                                indptr[r0:r1 + 1] - indptr[r0]), shape=(r1 - r0, g))
            oh = sp.csr_matrix((np.ones(int(keep.sum())), (grp[keep], np.flatnonzero(keep))),
                               shape=(n_groups, r1 - r0))
            acc += np.asarray((oh @ Xb).todense())
            if (r0 // chunk) % 40 == 0:
                log(f"    rows {r0:>8,}/{n:,}  {time.time() - t0:5.0f}s  peak {peak_gib():.2f} GiB")
    cnt = np.bincount(groups[groups >= 0], minlength=n_groups)
    assert np.array_equal(acc, np.rint(acc)), "non-integer count sums"
    return acc, cnt, g


def effect(S: np.ndarray, S_ctrl: np.ndarray) -> np.ndarray:
    """THE FROZEN RESPONSE FORMULA (PREREGISTRATION sec 4), 285's `lib285.delta` on one row pair:

        e(g) = log1p(50000 * S[g] / sum_j S[j]) - log1p(50000 * Sc[g] / sum_j Sc[j])

    The normalising sum runs over ALL genes of the library (the full library total), whatever
    subset of columns is afterwards inspected.
    """
    return L.prof(np.asarray(S, dtype=np.float64)) - L.prof(np.asarray(S_ctrl, dtype=np.float64))


def dump(obj, path: Path) -> None:
    path.write_text(json.dumps(obj, indent=1, default=float))
