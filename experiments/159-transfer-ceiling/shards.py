"""159 §2 -- per-target COUNT SHARDS, streamed, so any disjoint sub-sample of a target's
cells can be pseudobulked without re-reading the h5ad.

Each target's cells are randomly partitioned (deterministic, SEED below) into K near-equal
shards; the raw counts of each shard are summed.  Any two disjoint cell subsets of size
m = j * shard_size are then just sums of j shards each, so the whole split-half /
sub-sample-size curve is a few numpy adds.

Why shards and not one 50/50 split: the brief needs r as a FUNCTION of cells per target, in
order to extrapolate the 150-cell Orion harnesses to the challenge's 400.  Sub-sampling
inside a surface gives that curve from the surface's own data, with no cross-target
confound (every Orion target has exactly 150 cells, so binning by cell count is impossible).

Control cells are sharded the same way, so a half-effect can be built against a DISJOINT
half of the control.  Sharing one control between the two halves correlates their noise and
inflates the split-half cosine; that inflation is measured in `reliability.py` by running
both variants.

Memory: streams CSR row blocks with h5py; peak is the accumulator,
(n_groups * K) x G float64.  Largest here is h1-full at 255*20+20 rows x 18,080 -> 0.74 GB.
"""
from __future__ import annotations

import os
import sys

import h5py
import numpy as np
import scipy.sparse as sp

SEED = 20260902
CHUNK = 4000
CTRL = "non-targeting"
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")

# surface -> (path, target shards, control shards)
SURFACES = {
    "hct116": ("data/harness/orion-hct116-thinned.h5ad", 10, 20),
    "hek293t": ("data/harness/orion-hek293t-thinned.h5ad", 10, 20),
    "h1test": ("data/harness/h1-test.h5ad", 20, 20),
    "h1full": ("data/harness/h1-full.h5ad", 20, 20),
}


def _read_str(h, path):
    o = h[path]
    if hasattr(o, "keys") and "categories" in o:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c)
                           for c in o["categories"][:]])
        return cats[o["codes"][:]]
    v = o["values"][:] if hasattr(o, "keys") else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def build(name: str, rel: str, kt: int, kc: int) -> None:
    out = os.path.join(WORK, f"shards-{name}.npz")
    if os.path.exists(out):
        print(f"[{name}] exists, skipping")
        return
    path = os.path.join(REPO, rel)
    with h5py.File(path, "r") as h:
        n, g = (int(x) for x in h["X"].attrs["shape"])
        genes = _read_str(h, "var/_index")
        labels = _read_str(h, "obs/target_gene")
        assert CTRL in set(labels), f"{name}: no {CTRL} group"
        uniq = np.asarray(sorted(set(labels)))
        gid = {u: i for i, u in enumerate(uniq)}
        gcode = np.asarray([gid[x] for x in labels], dtype=np.int64)
        ci = gid[CTRL]

        # deterministic shard assignment, per group, seeded once
        rng = np.random.default_rng(SEED)
        shard = np.full(n, -1, dtype=np.int64)
        for u in uniq:
            j = int(gid[u])
            rows = np.flatnonzero(gcode == j)
            k = kc if j == ci else kt
            perm = rng.permutation(rows.size)
            shard[rows[perm]] = np.arange(rows.size) % k

        kmax = max(kt, kc)
        row = gcode * kmax + shard          # accumulator row index
        nrow = len(uniq) * kmax
        acc = np.zeros((nrow, g), dtype=np.float64)

        indptr = h["X/indptr"][:]
        for r0 in range(0, n, CHUNK):
            r1 = min(r0 + CHUNK, n)
            s, e = int(indptr[r0]), int(indptr[r1])
            X = sp.csr_matrix((h["X/data"][s:e].astype(np.float64), h["X/indices"][s:e],
                               indptr[r0:r1 + 1] - indptr[r0]), shape=(r1 - r0, g))
            oh = sp.csr_matrix((np.ones(r1 - r0), (row[r0:r1], np.arange(r1 - r0))),
                               shape=(nrow, r1 - r0))
            acc += np.asarray((oh @ X).todense())
            print(f"    rows {r1}/{n}", end="\r", flush=True)
        print()

    ncell = np.bincount(row, minlength=nrow)
    acc = acc.reshape(len(uniq), kmax, g)
    ncell = ncell.reshape(len(uniq), kmax)

    keep = np.asarray([i for i in range(len(uniq)) if i != ci])
    np.savez(out,
             targets=uniq[keep], genes=genes,
             tsum=acc[keep, :kt].astype(np.float32), tcell=ncell[keep, :kt],
             csum=acc[ci, :kc].astype(np.float32), ccell=ncell[ci, :kc],
             kt=np.int64(kt), kc=np.int64(kc), seed=np.int64(SEED))
    print(f"[{name}] {len(keep)} targets x {g} genes | {kt} target shards of "
          f"{ncell[keep, :kt].mean():.1f} cells | {kc} control shards of "
          f"{ncell[ci, :kc].mean():.1f} cells -> {out}")


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    want = sys.argv[1:] or list(SURFACES)
    for name in want:
        rel, kt, kc = SURFACES[name]
        print(f"[{name}] {rel}")
        build(name, rel, kt, kc)


if __name__ == "__main__":
    main()
