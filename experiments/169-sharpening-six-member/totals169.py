"""169 — cache the h1-test prediction's per-target gene totals, streamed.

`realize164.py` loads the whole 2.0 GB `pred.h5ad` with `anndata` and densifies each target
block just to take a column sum.  That is the 8.8 GB peak RSS 164 reported, and two sibling
agents are running.  Everything downstream needs only the per-gene panel totals (the quasi
emitter preserves them exactly, which is the whole reason 093/095/157/164 can re-emit from a
built prediction), so stream them once with h5py and cache.

Also caches the harness control column sums, which is what production's `b_ctrl` is built
from (`make_submission_pm_cd4.py:2056`: `b_ctrl = log1p(TS * totals / totals.sum())`).

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import os
import sys
import time

import h5py
import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")

REAL = f"{ROOT}/data/harness/h1-test.h5ad"
PRED = (f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/"
        f"g0.65_q1.0_a1.08_spread0.25/pred.h5ad")
COL, CTRL = "target_gene", "non-targeting"


def _cat(h, key):
    g = h[f"obs/{key}"]
    if isinstance(g, h5py.Group):
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c)
                           for c in g["categories"][:]])
        return cats[g["codes"][:]]
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in g[:]])


def _genes(h):
    gv = h["var/_index"]
    gv = gv["values"] if isinstance(gv, h5py.Group) else gv
    return np.asarray([g.decode() if isinstance(g, bytes) else str(g) for g in gv[:]])


def group_sums(path, want=None):
    """Per-label column sums, streamed group by group out of a CSR h5ad."""
    t0 = time.time()
    with h5py.File(path, "r") as h:
        G = int(h["X"].attrs["shape"][1])
        genes = _genes(h)
        lab = _cat(h, COL)
        labels, codes = np.unique(lab, return_inverse=True)
        if want is not None:
            keep = np.isin(labels, list(want))
            labels = labels[keep]
        indptr = h["X/indptr"][:]
        data_ds, ind_ds = h["X/data"], h["X/indices"]
        order = np.argsort(codes, kind="stable")
        allb = np.unique(lab)
        bounds = np.searchsorted(codes[order], np.arange(allb.size + 1))
        pos = {p: i for i, p in enumerate(allb)}
        sums = np.zeros((labels.size, G), dtype=np.int64)
        ncell = np.zeros(labels.size, dtype=np.int64)
        for k, name in enumerate(labels):
            rows = np.sort(order[bounds[pos[name]]:bounds[pos[name] + 1]])
            ncell[k] = rows.size
            assert rows[-1] - rows[0] + 1 == rows.size, f"{name}: rows not contiguous"
            i0, i1 = int(rows[0]), int(rows[-1]) + 1
            ip = indptr[i0:i1 + 1]
            d = data_ds[ip[0]:ip[-1]]
            ix = ind_ds[ip[0]:ip[-1]]
            np.add.at(sums[k], ix.astype(np.int64), d.astype(np.int64))
            if k % 20 == 0:
                print(f"  {k}/{labels.size}  {time.time() - t0:.0f}s", flush=True)
    return labels, genes, sums, ncell


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    out = os.path.join(WORK, "h1test-totals.npz")
    if os.path.exists(out):
        print(f"{out} exists; nothing to do")
        return

    print("streaming the harness control block ...", flush=True)
    with h5py.File(REAL, "r") as h:
        rlab = _cat(h, COL)
        rgenes = _genes(h)
        Gr = int(h["X"].attrs["shape"][1])
        ci = np.flatnonzero(rlab == CTRL)
        indptr = h["X/indptr"][:]
        ctrl_totals = np.zeros(Gr, dtype=np.float64)
        # Control rows are NOT contiguous in this harness, so stream every row block and
        # accumulate only the control rows rather than assuming a span.
        isctrl = rlab == CTRL
        step = 20000
        for i in range(0, rlab.size, step):
            j = min(i + step, rlab.size)
            ip = indptr[i:j + 1]
            d = h["X/data"][ip[0]:ip[-1]]
            ix = h["X/indices"][ip[0]:ip[-1]]
            nn = np.diff(ip)
            rowid = np.repeat(np.arange(i, j), nn)
            m = isctrl[rowid]
            np.add.at(ctrl_totals, ix[m].astype(np.int64), d[m].astype(np.float64))
    print(f"  {ci.size} control cells, {Gr} genes, {ctrl_totals.sum():.3e} counts", flush=True)

    print("streaming the prediction ...", flush=True)
    labels, pgenes, sums, ncell = group_sums(PRED)
    np.savez_compressed(
        out, pred_targets=labels, pred_genes=pgenes, pred_sums=sums, pred_ncell=ncell,
        harness_genes=rgenes, ctrl_totals=ctrl_totals, n_ctrl=np.int64(ci.size),
    )
    print(f"wrote {out}: {labels.size} prediction groups x {pgenes.size} genes; "
          f"cells/group {np.unique(ncell)}")


if __name__ == "__main__":
    sys.exit(main())
