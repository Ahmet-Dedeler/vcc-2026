"""Per-(context, perturbation) group sums and delete-1 jackknife, streamed from a big h5ad.

A shipped prediction is 360,000 cells x 18,533 genes at ~2.5e9 nonzeros -- 30 GB dense-ish in
CSR, which anndata will happily try to hold. Everything the expression metric needs is
per-group, so read one group's CSR rows at a time straight out of h5py and never materialize
the panel.

    uv run python mech_chal_stats.py <pred.h5ad> <out.npz> [--ctx-col context]

Writes `perts` (context-tagged), `sums` [P, G] int64, `counts` [P], `jk` [P], `genes` [G].
"""
from __future__ import annotations

import sys
import time

import h5py
import numpy as np
from scipy.sparse import csr_matrix

sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/096-mse-floor")
from cell_eval2.moments import jackknife_correction  # noqa: E402

TS = 50_000.0


def _cat(h, key):
    g = h[f"obs/{key}"]
    if isinstance(g, h5py.Group):
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c)
                           for c in g["categories"][:]])
        return cats[g["codes"][:]]
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in g[:]])


def main(path, out, pert_col="target_gene", ctx_col="context", chunk=64):
    with h5py.File(path, "r") as h:
        shape = h["X"].attrs["shape"]
        G = int(shape[1])
        gv = h["var/_index"]
        gv = gv["values"] if isinstance(gv, h5py.Group) else gv
        genes = np.asarray([g.decode() if isinstance(g, bytes) else str(g) for g in gv[:]])
        lab = _cat(h, pert_col)
        ctx = _cat(h, ctx_col) if ctx_col and ctx_col in h["obs"] else np.full(lab.size, "")
        key = np.asarray([f"{c}|{t}" if c else t for c, t in zip(ctx, lab)])
        perts, codes = np.unique(key, return_inverse=True)
        P = perts.size
        sums = np.zeros((P, G), dtype=np.int64)
        counts = np.zeros(P, dtype=np.float64)
        jk = np.zeros(P, dtype=np.float64)
        indptr_all = h["X/indptr"][:]
        data_ds, ind_ds = h["X/data"], h["X/indices"]
        order = np.argsort(codes, kind="stable")
        bounds = np.searchsorted(codes[order], np.arange(P + 1))
        t0 = time.time()
        for p in range(P):
            rows = np.sort(order[bounds[p]:bounds[p + 1]])
            n = rows.size
            counts[p] = n
            # rows of one group are contiguous in every prediction this repo writes; fall back
            # to a gather if not, rather than silently reading the wrong span.
            if rows[-1] - rows[0] + 1 == n:
                i0, i1 = int(rows[0]), int(rows[-1]) + 1
                ip = indptr_all[i0:i1 + 1]
                d = data_ds[ip[0]:ip[-1]]
                ix = ind_ds[ip[0]:ip[-1]]
                ip = ip - ip[0]
            else:
                ds, ixs, ips = [], [], [0]
                for r in rows:
                    a, b = indptr_all[r], indptr_all[r + 1]
                    ds.append(data_ds[a:b]); ixs.append(ind_ds[a:b]); ips.append(ips[-1] + (b - a))
                d = np.concatenate(ds); ix = np.concatenate(ixs); ip = np.asarray(ips)
            X = csr_matrix((d.astype(np.float64), ix.astype(np.int64), ip.astype(np.int64)),
                           shape=(n, G))
            sums[p] = np.asarray(X.sum(axis=0)).ravel().astype(np.int64)
            jk[p] = jackknife_correction(X, np.zeros(n, np.intp), 1, TS, chunk=chunk)[0]
            if p % 50 == 0:
                print(f"  {p}/{P}  {time.time()-t0:.0f}s", flush=True)
    np.savez_compressed(out, perts=perts, sums=sums, counts=counts, jk=jk, genes=genes)
    print(f"wrote {out}  P={P} G={G}  {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], *(sys.argv[3:]))
