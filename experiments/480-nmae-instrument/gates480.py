"""480: reference DE tables with the scorer's rules on a cell-level harness (Wilcoxon rank-sum vs the harness's own
control cells, BH per perturbation over genes > 5 CPM in control, lfc = log2 of arithmetic mean CPM ratio).

Memory-light: the CSR matrix is loaded once (1.4-1.7 GB for the Orion harnesses) and processed in gene chunks, so no
dense cells x genes matrix is ever built. Wilcoxon z from 092's scanpy-identical kernel against pre-sorted controls.

    uv run python gates480.py <name> <h5ad> [--ctrl non-targeting] [--max-cells 150] [--seed 0]

Writes work/gate-<name>.npz: targets, genes (tested axis), ctrl_cpm, lib_med, ncell, n_ctrl, p, q, lfc (float32).
Known answer (orion-hct116-thinned): q and lfc equal 098's refde_orion.npz (built by an independent dense script)."""
from __future__ import annotations

import argparse
import sys
import time

import h5py
import numpy as np
import scipy.sparse as sp
from scipy.stats import false_discovery_control, norm

from common480 import HERE, ROOT

sys.path.insert(0, str(ROOT / "experiments/092-call-set-truncation"))
from de_kernel import ctrl_tiesums, wilcoxon_z  # noqa: E402


def read_str(o):
    """anndata string column: categorical, nullable-string-array or plain."""
    if isinstance(o, h5py.Group):
        if "categories" in o:
            return o["categories"][:].astype(str)[o["codes"][:]]
        return o["values"][:].astype(str)
    return o[:].astype(str)


def read_labels(h, key="target_gene"):
    return read_str(h["obs"][key])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("h5ad")
    ap.add_argument("--ctrl", default="non-targeting")
    ap.add_argument("--max-cells", type=int, default=0, help="subsample each perturbation to this many cells (0 = all)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--chunk", type=int, default=1000)
    ap.add_argument("--blk", type=int, default=4096, help="rows per streamed HDF5 block (peak memory)")
    a = ap.parse_args()
    t0 = time.time()
    h = h5py.File(a.h5ad, "r")
    lab = read_labels(h)
    genes = read_str(h["var"]["_index"])
    X = h["X"]
    shape = tuple(int(v) for v in X.attrs["shape"])
    indptr = X["indptr"][:].astype(np.int64)
    dset, iset = X["data"], X["indices"]
    BLK = a.blk

    def row_blocks():
        """Yield (lo, hi, CSR block) streamed from HDF5; peak memory is one block."""
        for lo in range(0, shape[0], BLK):
            hi = min(lo + BLK, shape[0])
            a0, a1 = int(indptr[lo]), int(indptr[hi])
            yield lo, hi, sp.csr_matrix((dset[a0:a1], iset[a0:a1], indptr[lo:hi + 1] - a0), shape=(hi - lo, shape[1]))

    print(f"[{a.name}] {shape[0]} cells x {shape[1]} genes, nnz {int(indptr[-1]):,}", flush=True)
    lib = np.zeros(shape[0])
    for lo, hi, B in row_blocks():
        lib[lo:hi] = np.asarray(B.sum(1), dtype=np.float64).ravel()
    lib[lib == 0] = 1.0
    cidx = np.flatnonzero(lab == a.ctrl)
    rng = np.random.default_rng(a.seed)
    perts = sorted(set(lab) - {a.ctrl})
    rows = {}
    for p in perts:
        r = np.flatnonzero(lab == p)
        if a.max_cells and r.size > a.max_cells:
            r = np.sort(rng.choice(r, a.max_cells, replace=False))
        rows[p] = r
    w = np.zeros(shape[0])
    w[cidx] = 1e6 / lib[cidx] / cidx.size
    cpm_all = np.zeros(shape[1])
    for lo, hi, B in row_blocks():
        cpm_all += B.T @ w[lo:hi]
    tested = np.flatnonzero(cpm_all > 5.0)
    print(f"[{a.name}] {cidx.size} controls, {tested.size} genes > 5 CPM, {len(perts)} perturbations "
          f"({time.time() - t0:.0f}s)", flush=True)
    keep_rows = np.concatenate([cidx] + [rows[p] for p in perts])
    pos = np.full(shape[0], -1)
    pos[keep_rows] = np.arange(keep_rows.size)
    c_pos = pos[cidx]
    G = tested.size
    P = np.ones((len(perts), G))
    LFC = np.zeros((len(perts), G), dtype=np.float32)
    for lo in range(0, G, a.chunk):
        cols = tested[lo:lo + a.chunk]
        D = np.zeros((keep_rows.size, cols.size), dtype=np.float64)       # counts of the kept cells, this gene chunk
        for r0, r1, B in row_blocks():
            sel = pos[r0:r1] >= 0
            if sel.any():
                D[pos[r0:r1][sel]] = B[np.flatnonzero(sel)][:, cols].toarray()
        D *= (1e6 / lib[keep_rows])[:, None]
        Cc = D[c_pos]
        cmean = Cc.mean(0)
        cs = np.ascontiguousarray(np.sort(np.log1p(Cc).T, axis=1))
        tie = ctrl_tiesums(cs)
        del Cc
        for i, p in enumerate(perts):
            Bv = D[pos[rows[p]]]
            z = wilcoxon_z(cs, tie, np.ascontiguousarray(np.log1p(Bv).T))
            P[i, lo:lo + cols.size] = 2 * norm.sf(np.abs(z))
            LFC[i, lo:lo + cols.size] = np.log2((Bv.mean(0) + 1e-9) / (cmean + 1e-9))
        del D
        print(f"  genes {lo + cols.size}/{G} ({time.time() - t0:.0f}s)", flush=True)
    h.close()
    Q = np.vstack([false_discovery_control(P[i], method="bh") for i in range(len(perts))])
    gk = genes[tested]
    ctrl_cpm = cpm_all[tested]
    ncell = np.array([rows[p].size for p in perts])
    sig = [(Q[i] < 0.05) & np.isfinite(LFC[i]) & (LFC[i] != 0) & (gk != p) for i, p in enumerate(perts)]
    sz = np.array([s.sum() for s in sig])
    print(f"[{a.name}] gate sizes median {np.median(sz):.0f} mean {sz.mean():.0f} min {sz.min()} max {sz.max()}; "
          f">=10: {(sz >= 10).sum()}/{len(sz)}; median lib ctrl {np.median(lib[cidx]):.0f}", flush=True)
    np.savez_compressed(HERE / f"work/gate-{a.name}.npz", targets=np.array(perts), genes=gk, ctrl_cpm=ctrl_cpm,
                        lib_med=np.median(lib[cidx]), lib_pert_med=np.median(lib[np.concatenate(list(rows.values()))]),
                        ncell=ncell, n_ctrl=cidx.size, p=P.astype(np.float32), q=Q.astype(np.float32), lfc=LFC)
    print(f"[{a.name}] wrote work/gate-{a.name}.npz ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
