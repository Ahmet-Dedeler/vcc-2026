"""285 -- CALIBRATION: how much of a true CELL split's disagreement is Poisson counting noise?

On the two sources that have cells (h1-test and K562-GWPS), the same targets are split two ways at the SAME
total-UMI budget per side:
  * CELL split   -- two disjoint sets of n cells
  * COUNT split  -- the pooled 2n-cell sum, halved by a Binomial(count, 1/2) per gene
The ratio cell/count is what makes the CD4 bank's count-split row (the only cell-sampling reference a
pseudobulk bank can carry) readable as a cell split.

    PYTHONDONTWRITEBYTECODE=1 uv run python csplit_285.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib285 as L  # noqa: E402

P = L.Logger(L.WORK / "csplit_285.log")


def labels_of(h, path):
    o = h[path]
    if isinstance(o, h5py.Group) and "categories" in o:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in o["categories"][:]])
        return cats[o["codes"][:]]
    v = o["values"][:] if isinstance(o, h5py.Group) else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def run(name, dense, targets, scored, ns, out):
    P(f"\n{name}: cell split vs count split at matched total UMIs, 8 draws")
    out[name] = {}
    for n in ns:
        ts = [t for t in targets if dense[t].shape[0] >= 2 * n]
        if len(ts) < 20:
            continue
        cs = np.full((L.NDRAW, len(ts)), np.nan); ks = np.full((L.NDRAW, len(ts)), np.nan)
        for d in range(L.NDRAW):
            r = np.random.default_rng(L.SEED + 17 * d + n)
            A, B, C1, C2 = [], [], [], []
            for t in ts:
                p = r.permutation(dense[t].shape[0])
                a = dense[t][p[:n]].sum(0); b = dense[t][p[n:2 * n]].sum(0)
                A.append(a); B.append(b)
                h = a + b                                     # the SAME 2n cells, split by molecule instead
                x = r.binomial(np.rint(h).astype(np.int64), 0.5).astype(np.float64)
                C1.append(x); C2.append(h - x)
            PA, PB = L.prof(np.stack(A))[:, scored], L.prof(np.stack(B))[:, scored]
            QA, QB = L.prof(np.stack(C1))[:, scored], L.prof(np.stack(C2))[:, scored]
            cs[d] = L.row_cos(PA - PA.mean(0, keepdims=True), PB - PB.mean(0, keepdims=True))
            ks[d] = L.row_cos(QA - QA.mean(0, keepdims=True), QB - QB.mean(0, keepdims=True))
        xc, xk = np.nanmean(cs, 0), np.nanmean(ks, 0)
        rat = xc / np.where(np.abs(xk) > 1e-9, xk, np.nan)
        out[name][n] = dict(n_targets=len(ts), cell=dict(**L.stats(xc), boot=L.boot_median(xc)),
                            count=dict(**L.stats(xk), boot=L.boot_median(xk)),
                            ratio_of_medians=float(np.nanmedian(xc) / np.nanmedian(xk)),
                            median_ratio=float(np.nanmedian(rat)),
                            diff=L.boot_median_diff(xk, xc))
        P(f"  n={n:4d}/side ({len(ts)} targets): CELL {L.fmt(out[name][n]['cell']['boot'])}   COUNT "
          f"{L.fmt(out[name][n]['count']['boot'])}   ratio of medians {out[name][n]['ratio_of_medians']:.3f}   "
          f"COUNT − CELL {L.fmt(out[name][n]['diff'])}")
        out[name][n]["arrays"] = dict(cell=xc.tolist(), count=xk.tolist(), targets=list(ts))


def main() -> None:
    t0 = time.time()
    P(f"285 csplit_285.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    out: dict = {}

    # h1-test
    H1 = L.ROOT / "data/harness/h1-test.h5ad"
    with h5py.File(H1, "r") as h:
        lab = labels_of(h, "obs/target_gene"); genes = labels_of(h, "var/_index")
        n_cells, n_genes = (int(x) for x in h["X"].attrs["shape"])
        X = sp.csr_matrix((h["X/data"][:].astype(np.float32), h["X/indices"][:], h["X/indptr"][:]),
                          shape=(n_cells, n_genes))
    targets = sorted(set(lab.tolist()) - {"non-targeting"})
    gpos = {g: i for i, g in enumerate(genes)}
    scored = np.asarray(sorted(set(range(n_genes)) - {gpos[t] for t in targets}), dtype=np.int64)
    dense = {t: np.asarray(X[np.flatnonzero(lab == t)].todense(), dtype=np.float64) for t in targets}
    del X
    run("h1_test", dense, targets, scored, (25, 50, 100, 200), out)
    del dense

    # K562-GWPS panel cells (already extracted)
    m = np.load(L.WORK / "k562_gwps_meta.npz", allow_pickle=True)
    genes = m["genes"].astype(str); pan_t = m["panel_targets"].astype(str); pert = m["pan_pert"].astype(str)
    gpos = {g: i for i, g in enumerate(genes)}
    scored = np.asarray(sorted(set(range(len(genes))) - {gpos[t] for t in pan_t if t in gpos}), dtype=np.int64)
    Xd = np.asarray(np.load(L.WORK / "k562_gwps_panelX.npy", mmap_mode="r"), dtype=np.float32)
    dense = {t: np.asarray(Xd[np.flatnonzero(pert == t)], dtype=np.float64) for t in pan_t}
    run("k562_gwps", dense, list(pan_t), scored, (25, 50, 100), out)

    out["seconds"] = round(time.time() - t0, 1)
    out["peak_rss_gib"] = L.peak_gib()
    L.Logger.dump(P, out, L.WORK / "csplit_285.json")
    P(f"\ndone {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB")
    P.save()


if __name__ == "__main__":
    main()
