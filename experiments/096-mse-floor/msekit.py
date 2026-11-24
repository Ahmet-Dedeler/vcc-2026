"""Scorer-exact `expr_mse_unbiased_capped_norm`, from group sums, with the legs exposed.

Everything here routes through `cell_eval2`'s OWN primitives -- `prep.bulk_lognorm_means`,
`moments.jackknife_correction`, `metrics.delta.mse_unbiased_capped`,
`metrics.delta.distance_unbiased` -- so this is not a reimplementation of the metric, it is
the metric called on cached sufficient statistics. The only thing added is the ability to
perturb the PREDICTION's group sums (e.g. rescale its effect by alpha) and re-score without
re-running the whole `031-somatic` pipeline.

Two cached artifacts per h5ad, keyed by path + mtime:
    sums[P, G]  int64 group count sums     -> bulk_lognorm_means gives the profile
    jk[P]       fp64 delete-1 jackknife    -> C_p / C_hat_p, the expensive half

`sums` is all the profile needs; `jk` is what costs O(n*G) dense.
"""
from __future__ import annotations

import hashlib
import os
import sys

import numpy as np
from scipy.sparse import csr_matrix, issparse

ROOT = "~/Code/virtual-cell-challenge-2026"
CACHE = f"{ROOT}/experiments/096-mse-floor/cache"
os.makedirs(CACHE, exist_ok=True)

from cell_eval2.metrics.delta import distance_unbiased, mse_unbiased_capped  # noqa: E402
from cell_eval2.moments import GroupMoments, jackknife_correction  # noqa: E402
from cell_eval2.prep import bulk_lognorm_means  # noqa: E402

TS = 50_000.0
CTRL = "non-targeting"


def _key(path: str, tag: str) -> str:
    st = os.stat(path)
    h = hashlib.sha1(f"{path}|{st.st_size}|{st.st_mtime_ns}|{tag}".encode()).hexdigest()[:16]
    return f"{CACHE}/{os.path.basename(path)}.{tag}.{h}.npz"


def group_stats(path: str, pert_col: str = "target_gene", *, with_jk: bool = True,
                chunk: int = 64, verbose: bool = True):
    """`(perts, sums[P, G], counts[P], jk[P] or None, genes[G])`, cached on disk.

    Reads the h5ad ONCE. `sums` is exact int64 (counts are integral by the submission rules,
    and every harness here is counts too), so a profile can be rebuilt at any alpha without
    touching the cells again.
    """
    ck = _key(path, f"gs.{pert_col}.{int(with_jk)}")
    if os.path.exists(ck):
        z = np.load(ck, allow_pickle=True)
        jk = z["jk"] if "jk" in z.files and z["jk"].size else None
        return z["perts"].astype(str), z["sums"], z["counts"], jk, z["genes"].astype(str)

    import anndata as ad
    if verbose:
        print(f"[msekit] reading {path}", flush=True)
    A = ad.read_h5ad(path)
    labels = np.asarray(A.obs[pert_col].astype(str).values)
    genes = np.asarray([str(g) for g in A.var_names])
    perts, codes = np.unique(labels, return_inverse=True)
    X = A.X.tocsr() if issparse(A.X) else csr_matrix(A.X)
    P, G = perts.size, X.shape[1]
    sums = np.zeros((P, G), dtype=np.int64)
    counts = np.zeros(P, dtype=np.float64)
    order = np.argsort(codes, kind="stable")
    bounds = np.searchsorted(codes[order], np.arange(P + 1))
    for p in range(P):
        rows = order[bounds[p]:bounds[p + 1]]
        counts[p] = rows.size
        sums[p] = np.asarray(X[rows].sum(axis=0)).ravel().astype(np.int64)
    jk = None
    if with_jk:
        if verbose:
            print(f"[msekit] jackknife over {P} groups x {G} genes", flush=True)
        jk = jackknife_correction(X, codes, P, TS, chunk=chunk)
    np.savez_compressed(ck, perts=perts, sums=sums, counts=counts,
                        jk=jk if jk is not None else np.zeros(0), genes=genes)
    return perts, sums, counts, jk, genes


def score(pred_perts, pred_sums, pred_counts, pred_jk,
          real_perts, real_sums, real_counts, real_jk, genes, *, pert_col="target_gene"):
    """`(mse, num_dict, den_dict)` -- `expr_mse_unbiased_capped_norm` and both legs.

    The aggregation is `ratio_of_sums`, exactly as `catalog.py` declares it.
    """
    pm = bulk_lognorm_means(pred_sums.astype(np.float64), TS)
    rm = bulk_lognorm_means(real_sums.astype(np.float64), TS)
    pmom = GroupMoments(perts=pred_perts, counts=pred_counts,
                        sumsq=np.zeros_like(pred_counts), jk=pred_jk)
    rmom = GroupMoments(perts=real_perts, counts=real_counts,
                        sumsq=np.zeros_like(real_counts), jk=real_jk)
    kw = dict(comparator="bulk_lognorm", pert_col=pert_col, control=CTRL,
              genes=genes, target_gene_map=None)
    num = mse_unbiased_capped(pred_bulk=(pred_perts, pm), real_bulk=(real_perts, rm),
                              pred_moments=pmom, real_moments=rmom, **kw)
    den = distance_unbiased(real_bulk=(real_perts, rm), real_moments=rmom, **kw)
    # denominator is real-side-only and spans EVERY non-control real row; the numerator spans
    # the predicted rows. `ratio_of_sums` in `partition.py` sums the rows each leg emits.
    keys = [k for k in num if k in den]
    n = float(np.sum([num[k] for k in keys]))
    d = float(np.sum([den[k] for k in keys]))
    return n / d, num, den


if __name__ == "__main__":
    real, pred = sys.argv[1], sys.argv[2]
    rp, rs, rc, rj, rg = group_stats(real)
    pp, ps, pc, pj, pg = group_stats(pred)
    assert np.array_equal(rg, pg), (rg.size, pg.size)
    m, num, den = score(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    print(f"expr_mse_unbiased_capped_norm = {m:.10f}   "
          f"num_mean={np.mean(list(num.values())):.10g} den_mean={np.mean(list(den.values())):.10g}")
