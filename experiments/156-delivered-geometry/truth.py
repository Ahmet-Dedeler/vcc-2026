"""156 — build exact-space pseudobulk truth for the two Orion held surfaces.

WHY
---
`pds_cosine` and `expr_mse` read the profile
    b_p = log1p(TS * P_p / sum_g P_p),      TS = 5e4
with P_p the summed raw counts of the group's cells, and the effect
    delta_p = b_p - b_ctrl.

Every prior effect-level screen in this repo used the DE table's `real_lfc`
(log2 of per-cell CPM arithmetic means) as a stand-in for that. It is a different
functional of the same counts -- 153's own docstring says the two disagree, and
155 could not reproduce a quoted `pds` to better than 0.005 because of it.
Since this experiment is specifically about a *profile-space* transform, the
truth has to be in the profile space the scorer uses. Nothing here is a proxy.

Output per surface: targets x genes exact-space effect, plus the control profile.
"""
from __future__ import annotations

import os
import sys

import h5py
import numpy as np
import scipy.sparse as sp

TS = 5.0e4
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")

SURFACES = {
    "hct116": "data/harness/orion-hct116-thinned.h5ad",
    "hek293t": "data/harness/orion-hek293t-thinned.h5ad",
}
CTRL_LABEL = "non-targeting"
CHUNK = 4000


def _read_str(h, path):
    o = h[path]
    if hasattr(o, "keys") and "categories" in o:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c)
                           for c in o["categories"][:]])
        return cats[o["codes"][:]]
    v = o["values"][:] if hasattr(o, "keys") else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def group_sums(path: str):
    """Sum raw counts per target group, streaming CSR row blocks."""
    with h5py.File(path, "r") as h:
        n, g = (int(x) for x in h["X"].attrs["shape"])
        genes = _read_str(h, "var/_index")
        labels = _read_str(h, "obs/target_gene")
        uniq = np.asarray(sorted(set(labels)))
        gid = {u: i for i, u in enumerate(uniq)}
        codes = np.asarray([gid[x] for x in labels], dtype=np.int64)

        indptr = h["X/indptr"][:]
        acc = np.zeros((len(uniq), g), dtype=np.float64)
        for r0 in range(0, n, CHUNK):
            r1 = min(r0 + CHUNK, n)
            s, e = int(indptr[r0]), int(indptr[r1])
            data = h["X/data"][s:e].astype(np.float64)
            idx = h["X/indices"][s:e]
            ip = indptr[r0:r1 + 1] - indptr[r0]
            X = sp.csr_matrix((data, idx, ip), shape=(r1 - r0, g))
            rows = codes[r0:r1]
            onehot = sp.csr_matrix(
                (np.ones(r1 - r0), (rows, np.arange(r1 - r0))),
                shape=(len(uniq), r1 - r0))
            acc += np.asarray((onehot @ X).todense())
            print(f"    rows {r1}/{n}", end="\r", flush=True)
        cells = np.bincount(codes, minlength=len(uniq))
    print()
    return uniq, genes, acc, cells


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    for name, rel in SURFACES.items():
        out = os.path.join(WORK, f"truth-{name}.npz")
        if os.path.exists(out):
            print(f"[{name}] exists, skipping")
            continue
        print(f"[{name}] {rel}")
        uniq, genes, acc, cells = group_sums(os.path.join(REPO, rel))
        assert CTRL_LABEL in set(uniq), f"no {CTRL_LABEL} group in {name}"

        prof = np.log1p(TS * acc / acc.sum(1, keepdims=True))
        ci = int(np.flatnonzero(uniq == CTRL_LABEL)[0])
        b_ctrl = prof[ci]
        keep = np.asarray([i for i in range(len(uniq)) if i != ci])
        targets = uniq[keep]
        effect = prof[keep] - b_ctrl[None, :]

        np.savez_compressed(
            out, targets=targets, genes=genes,
            effect=effect.astype(np.float32),
            b_ctrl=b_ctrl.astype(np.float32),
            ctrl_counts=acc[ci].astype(np.float64),
            cells=cells[keep].astype(np.int64),
            ctrl_cells=np.int64(cells[ci]), ts=np.float64(TS))
        print(f"[{name}] {len(targets)} targets x {len(genes)} genes"
              f" | median |effect| L2 {np.median(np.linalg.norm(effect, axis=1)):.4f}"
              f" | ctrl cells {cells[ci]} | median target cells {np.median(cells[keep]):.0f}")
        print(f"[{name}] -> {out}")


if __name__ == "__main__":
    sys.exit(main())
