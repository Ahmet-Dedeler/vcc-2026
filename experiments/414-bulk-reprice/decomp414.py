"""414 (c): read A_c and B_c straight off the shipped artifacts (395, 410 welded, 410b ballast), in the scorer's own space.

Same construction as 393's decomp393.py (which it reproduces on 319p2 as the known-answer row): stream each prediction
per (context, target) block with h5py bincount (never densified), mu = log1p(5e4 * tot / sum), d = mu - b_ctrl(context),
w = cell_eval2's `_row_weights`, A_c = sum_p w_p (|d_p|^2 - d_p[target]^2), B_c = cell_eval2's `_across_pert_budget`.
These two need OUR file and the real control only; claim_c / Kc_c / sumD / X come from 187's retrodicted model.

    uv run python decomp414.py            (writes work/decomp414.json)
"""
from __future__ import annotations

import json
import os
import resource
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
os.environ.setdefault("NUMBA_CACHE_DIR", str(HERE / "work" / "numba-cache"))
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import h5py
import numpy as np

ROOT = HERE.parents[1]
WORK = HERE / "work"
TS = 50000.0
ARTS = {
    "395": ROOT / "experiments/395-protected-signs-arm/work/build/the-thermostat-has-a-union-rep.h5ad",
    "410": ROOT / "experiments/410-pertscale-revisit/work/build/the-houseplant-filed-its-taxes.h5ad",
    "410b": ROOT / "experiments/410-pertscale-revisit/work/build_b/a-spreadsheet-achieved-enlightenment.h5ad",
}
# 393's 319p2 artifact was swept, so its exact KA row cannot be re-run. Known-answer rows instead:
#   KA-a  block streaming == an independent anndata-backed dense read of 3 blocks, exactly
#   KA-b  sum(w) == 300/18532 exactly (all 300 targets resolve; 393's control)
#   KA-c  (band, not exact) B/A within 393's 319p2 band 0.9913-0.9918 +- 0.003 for 395 (same family, earlier profile)
KA393 = json.loads((ROOT / "experiments/393-mse-shared-component/work/decomp393.json").read_text())["per_context"]


def str_col(h, path):
    obj = h[path]
    if isinstance(obj, h5py.Group):
        if "categories" in obj:
            return np.asarray(obj["categories"][:]).astype(str)[np.asarray(obj["codes"][:])]
        return np.asarray(obj["values"][:]).astype(str)
    return np.asarray(obj[:]).astype(str)


def controls(genes):
    out = {}
    for ctx in "ABC":
        with h5py.File(ROOT / f"data/challenge/context_{ctx}.h5ad", "r") as h:
            assert np.array_equal(str_col(h, "var/_index"), genes), f"gene axis mismatch on {ctx}"
            pooled = np.bincount(np.asarray(h["X"]["indices"][:], dtype=np.int64),
                                 weights=np.asarray(h["X"]["data"][:], dtype=np.float64), minlength=genes.size)
        out[ctx] = np.log1p(TS * pooled / pooled.sum())
    return out


def decompose(path, bctrl_cache):
    from cell_eval2.metrics.delta import _across_pert_budget, _row_weights
    t0 = time.time()
    with h5py.File(path, "r") as h:
        genes = str_col(h, "var/_index")
        ctx_lab = str_col(h, "obs/context")
        tgt_lab = str_col(h, "obs/target_gene")
        indptr = np.asarray(h["X"]["indptr"][:])
        lab = list(zip(ctx_lab.tolist(), tgt_lab.tolist()))
        blocks, start = [], 0
        for i in range(1, len(lab) + 1):
            if i == len(lab) or lab[i] != lab[start]:
                blocks.append((lab[start], start, i))
                start = i
        totals = np.zeros((len(blocks), genes.size), dtype=np.float64)
        for k, (_, lo_row, hi_row) in enumerate(blocks):
            lo, hi = int(indptr[lo_row]), int(indptr[hi_row])
            totals[k] = np.bincount(h["X"]["indices"][lo:hi].astype(np.int64),
                                    weights=h["X"]["data"][lo:hi].astype(np.float64), minlength=genes.size)
    if "b" not in bctrl_cache:
        bctrl_cache["b"] = controls(genes)
    b_ctrl = bctrl_cache["b"]
    keys = [k for k, _, _ in blocks]
    gpos = {g: i for i, g in enumerate(genes)}
    out = {}
    for ctx in "ABC":
        rows = [i for i, k in enumerate(keys) if k[0] == ctx]
        tot = totals[rows]
        mu = np.log1p(TS * tot / tot.sum(1, keepdims=True))
        dh = mu - b_ctrl[ctx][None, :]
        cols = {j: gpos[keys[i][1]] for j, i in enumerate(rows) if keys[i][1] in gpos}
        w = _row_weights(len(rows), genes.size, cols)
        a = np.einsum("ij,ij->i", dh, dh)
        for j, col in cols.items():
            a[j] -= dh[j, col] ** 2
        A = float((a * w).sum())
        B = float(_across_pert_budget(mu, cols, w))
        assert abs(float(w.sum()) - 300 / 18532) < 1e-12 and len(cols) == 300, "KA-b FAILED"
        out[ctx] = dict(A=A, B=B, BoverA=B / A, n_rows=len(rows), n_resolved=len(cols), wsum=float(w.sum()),
                        row_energy_q=np.quantile(a, [0.1, 0.5, 0.9]).tolist(),
                        block_total_median=float(np.median(tot.sum(1))))
        print(f"  {ctx}: A {A:.6f} B {B:.6f} B/A {B / A:.5f} rows {len(rows)} "
              f"[{time.time() - t0:.0f}s]", flush=True)
    # KA-a: three blocks re-read independently through anndata (backed), dense, summed
    import anndata as ad
    adata = ad.read_h5ad(path, backed="r")
    rng = np.random.default_rng(414)
    worst = 0.0
    for k in rng.choice(len(blocks), 3, replace=False):
        _, lo_row, hi_row = blocks[int(k)]
        dense = np.asarray(adata.X[lo_row:hi_row].toarray(), dtype=np.float64).sum(0)
        worst = max(worst, float(np.abs(dense - totals[int(k)]).max()))
    adata.file.close()
    print(f"  KA-a block stream vs anndata dense read, 3 blocks: worst |d| {worst:.1e}", flush=True)
    assert worst == 0.0, "KA-a FAILED"
    out["ka_block_worst"] = worst
    return out


def main() -> None:
    res, cache = {}, {}
    for name, path in ARTS.items():
        print(f"[{name}] {path.name}", flush=True)
        res[name] = decompose(path, cache)
        if name == "395":
            band = {c: (res[name][c]["BoverA"], KA393[c]["B"] / KA393[c]["A"]) for c in "ABC"}
            print(f"  KA-c band (395 vs 393's 319p2 B/A): {band}", flush=True)
            res[name]["ka_band_vs_319p2"] = band
    res["peak_rss_gb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    (WORK / "decomp414.json").write_text(json.dumps(res, indent=1))
    print(f"peak RSS {res['peak_rss_gb']:.2f} GB", flush=True)


if __name__ == "__main__":
    sys.exit(main())
