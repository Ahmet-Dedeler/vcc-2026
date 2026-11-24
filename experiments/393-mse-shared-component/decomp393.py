#!/usr/bin/env python3
"""393 Task 1: shared-vs-centred decomposition of the shipped 319p2 artifact.

Streams experiments/319-cap-variant/work/build_p2/the-photocopier-demands-a-raise.h5ad
per (context, target) block (h5py bincount; never densifies), builds the scorer-space
pseudobulk effects d = log1p(TS*tot/sum) - b_ctrl (TS=5e4, b_ctrl from the pooled real
controls of that context), and decomposes with the scorer's own _row_weights /
_across_pert_budget. Writes only into this experiment's work/.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", str(Path(__file__).resolve().parent / "work" / "numba-cache"))
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "src"))

TS = 50000.0
ART = ROOT / "experiments/319-cap-variant/work/build_p2/the-photocopier-demands-a-raise.h5ad"


def str_col(h, path):
    obj = h[path]
    if isinstance(obj, h5py.Group):
        if "categories" in obj:
            cats = np.asarray(obj["categories"][:]).astype(str)
            return cats[np.asarray(obj["codes"][:])]
        if "values" in obj:
            return np.asarray(obj["values"][:]).astype(str)
        raise ValueError(f"string group without categories/values: {path}")
    return np.asarray(obj[:]).astype(str)


def main() -> None:
    t0 = time.time()
    lines = []
    P = lambda *a: (print(*a, flush=True), lines.append(" ".join(str(x) for x in a)))
    import resource

    from cell_eval2.metrics.delta import _across_pert_budget, _row_weights

    with h5py.File(ART, "r") as h:
        art_genes = str_col(h, "var/_index")
        ctx_lab = str_col(h, "obs/context")
        tgt_lab = str_col(h, "obs/target_gene")
        indptr = np.asarray(h["X"]["indptr"][:])
        shape = tuple(h["X"].attrs["shape"])
        Xdata = h["X"]["data"]
        Xind = h["X"]["indices"]
        P(f"artifact shape {shape}, genes {art_genes.size}, "
          f"contexts {sorted(set(ctx_lab.tolist()))}")
        # block boundaries on (context, target) runs
        lab = list(zip(ctx_lab.tolist(), tgt_lab.tolist()))
        blocks, start = [], 0
        for i in range(1, len(lab) + 1):
            if i == len(lab) or lab[i] != lab[start]:
                blocks.append((lab[start], start, i))
                start = i
        sizes = {hi - lo for _, lo, hi in blocks}
        P(f"blocks {len(blocks)}, block sizes {sizes}")
        n_blocks, G = len(blocks), art_genes.size
        totals = np.zeros((n_blocks, G), dtype=np.float64)
        for k, (_key, lo_row, hi_row) in enumerate(blocks):
            lo, hi = int(indptr[lo_row]), int(indptr[hi_row])
            idx = Xind[lo:hi].astype(np.int64)
            val = Xdata[lo:hi].astype(np.float64)
            if idx.size:
                totals[k] = np.bincount(idx, weights=val, minlength=G)
            if (k + 1) % 300 == 0:
                P(f"  streamed {k + 1}/{n_blocks} [{time.time() - t0:.0f}s]")
        keys = [k for k, _, _ in blocks]
    P(f"stream done [{time.time() - t0:.0f}s]; totals {totals.shape}, "
      f"sum {totals.sum():.4g}")

    # b_ctrl per context from the pooled REAL controls; check gene-axis alignment
    # (h5py bincount straight off the CSR store; never densifies)
    b_ctrl, ctrl_n = {}, {}
    for ctx in ("A", "B", "C"):
        with h5py.File(ROOT / f"data/challenge/context_{ctx}.h5ad", "r") as h:
            cgenes = str_col(h, "var/_index")
            same_order = bool(np.array_equal(cgenes, art_genes))
            P(f"ctx {ctx}: control genes {cgenes.size}, identical order to artifact: {same_order}")
            if not same_order:
                raise SystemExit(f"gene axis mismatch on context {ctx}; refusing to proceed unaligned")
            idx = np.asarray(h["X"]["indices"][:], dtype=np.int64)
            val = np.asarray(h["X"]["data"][:], dtype=np.float64)
            pooled = np.bincount(idx, weights=val, minlength=G)
            ctrl_n[ctx] = int(h["X"].attrs["shape"][0])
            b_ctrl[ctx] = np.log1p(TS * pooled / pooled.sum())
            del idx, val, pooled
    P(f"controls n_cells {ctrl_n} [{time.time() - t0:.0f}s]")

    gpos = {g: i for i, g in enumerate(art_genes)}
    out = {}
    dbars = {}
    for ctx in ("A", "B", "C"):
        rows = [i for i, k in enumerate(keys) if k[0] == ctx]
        P(f"\nctx {ctx}: {len(rows)} rows")
        tot = totals[rows]
        row_sums = tot.sum(axis=1)
        P(f"  per-block totals: min {row_sums.min():.0f} med {np.median(row_sums):.0f} "
          f"max {row_sums.max():.0f}")
        mu = np.log1p(TS * tot / row_sums[:, None])
        dh = mu - b_ctrl[ctx][None, :]
        tgts = [keys[i][1] for i in rows]
        cols = {j: gpos[t] for j, t in enumerate(tgts) if t in gpos}
        P(f"  resolved {len(cols)}/{len(rows)} targets on the gene axis")
        w = _row_weights(len(rows), G, cols)
        a = np.einsum("ij,ij->i", dh, dh)
        for j, col in cols.items():
            a[j] -= dh[j, col] ** 2
        A = float((a * w).sum())
        B = float(_across_pert_budget(mu, cols, w))
        # shared component: weight-mean effect, and the A-B identity check
        wsum = float(w.sum())
        dbar = (dh * w[:, None]).sum(axis=0) / wsum
        P_shared = float(wsum * float(dbar @ dbar))
        share_identity = (A - B) / A
        share_direct = P_shared / A
        out[ctx] = dict(A=A, B=B, wsum=wsum, n_rows=len(rows), n_resolved=len(cols),
                        A_minus_B=A - B, P_dbar2=P_shared,
                        shared_share_identity=share_identity,
                        shared_share_direct=share_direct,
                        dbar_norm=float(np.linalg.norm(dbar)))
        dbars[ctx] = dbar
        P(f"  A {A:.6f}  B {B:.6f}  B/A {B / A:.5f}  wsum {wsum:.8f}")
        P(f"  A-B {A - B:.6f}  P||dbar||^2 {P_shared:.6f}  "
          f"shared share: identity {(A - B) / A:.4%} vs direct {P_shared / A:.4%}")
        del tot, mu, dh

    P("\ncos(dbar) across lines:")
    cos = {}
    for x, y in (("A", "B"), ("A", "C"), ("B", "C")):
        c = float(dbars[x] @ dbars[y] / (np.linalg.norm(dbars[x]) * np.linalg.norm(dbars[y])))
        cos[f"{x}{y}"] = c
        P(f"  cos(dbar_{x}, dbar_{y}) = {c:.4f}")

    # KA-adjacent sanity vs 211 rc_callset
    ref = {"A": (0.236429781363753, 0.23433559941175167),
           "B": (0.22642509487478957, 0.2243262456510523),
           "C": (0.23662010621111196, 0.23447576430467285)}
    P("\nsanity vs 211 rc_callset A_c/B_c (319p2 = later profile, expect few %):")
    for ctx in ("A", "B", "C"):
        dA = out[ctx]["A"] / ref[ctx][0] - 1
        dB = out[ctx]["B"] / ref[ctx][1] - 1
        P(f"  {ctx}: dA {dA:+.2%}  dB {dB:+.2%}")

    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    P(f"\npeak RSS {rss:.2f} GB, wall {(time.time() - t0) / 60:.1f} min")
    res = dict(per_context={c: dict(out[c]) for c in out},
                 cos=cos, peak_rss_gb=rss)
    (WORK / "decomp393.json").write_text(json.dumps(res, indent=1, default=float))
    # dbars are 3x18533 floats (1.3 MB) — keep for the pricing step
    np.savez(WORK / "dbars393.npz", **{f"dbar_{c}": dbars[c] for c in dbars},
             genes=np.asarray(art_genes))
    (WORK / "decomp393.log").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
