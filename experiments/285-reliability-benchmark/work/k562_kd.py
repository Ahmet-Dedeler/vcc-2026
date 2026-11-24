"""285 -- is the GUIDE-CONSTRUCT term knockdown EFFICIENCY, or is it off-target biology?

If two constructs against one gene disagree only because one of them knocks the gene down weakly, the guide
term is a dose artefact and better guide selection removes it.  If they disagree at matched, strong
knockdown, it is not removable.  Measured on K562-GWPS's 775 two-construct perturbations: each construct's
ON-TARGET knockdown is its own log-fold-change at the perturbed gene's own column (a column excluded from
every scored statistic), and the guide agreement is stratified by the WEAKER of the two constructs' KD.

    PYTHONDONTWRITEBYTECODE=1 uv run python k562_kd.py
"""
from __future__ import annotations
import sys, time
from collections import defaultdict
from pathlib import Path
import h5py, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib285 as L

PATH = L.ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
P = L.Logger(L.WORK / "k562_kd.log")
BLOCK = 8192
N = 25


def labels_of(h, path):
    o = h[path]
    if isinstance(o, h5py.Group) and "categories" in o:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in o["categories"][:]])
        return cats[o["codes"][:]]
    v = o["values"][:] if isinstance(o, h5py.Group) else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def main():
    t0 = time.time()
    P(f"285 k562_kd.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    with h5py.File(PATH, "r") as h:
        pert = labels_of(h, "obs/perturbation"); guide = labels_of(h, "obs/guide_id")
        genes = labels_of(h, "var/gene_name")
        n_cells, n_genes = h["X"].shape
    pg = defaultdict(set)
    for t, g in zip(pert, guide):
        pg[t].add(g)
    gpos = {g: i for i, g in enumerate(genes)}
    multi = sorted([t for t in pg if t != "control" and len(pg[t]) == 2 and t in gpos])
    keep = np.isin(pert, multi) | (pert == "control")
    idx = np.flatnonzero(keep)
    P(f"  {len(multi)} two-construct perturbations whose own gene is on the axis; extracting {idx.size:,} cells")
    with h5py.File(PATH, "r") as h:
        X = h["X"]; blocks = []
        for r0 in range(0, n_cells, BLOCK):
            r1 = min(r0 + BLOCK, n_cells)
            m = keep[r0:r1]
            if m.any():
                blocks.append(X[r0:r1][m].astype(np.uint16))
        Xd = np.concatenate(blocks); del blocks
    pert, guide = pert[idx], guide[idx]
    is_ctrl = pert == "control"
    ctrl = Xd[is_ctrl].sum(0, dtype=np.float64)
    pctrl = L.prof(ctrl)
    scored = np.asarray(sorted(set(range(n_genes)) - {gpos[t] for t in multi}), dtype=np.int64)
    info = {}
    for t in multi:
        ix = np.flatnonzero(pert == t)
        cs = sorted([ix[guide[ix] == g] for g in sorted(set(guide[ix].tolist()))], key=lambda a: -a.size)
        info[t] = (cs[0], cs[1])
    P(f"  extracted {Xd.shape}, {time.time()-t0:.0f}s, peak {L.peak_gib():.2f} GiB")

    ts = [t for t in multi if info[t][0].size >= 2 * N and info[t][1].size >= N]
    kd = np.full((len(ts), 2), np.nan)
    for i, t in enumerate(ts):
        for k in (0, 1):
            S = Xd[np.sort(info[t][k])].sum(0, dtype=np.float64)
            kd[i, k] = L.prof(S)[gpos[t]] - pctrl[gpos[t]]      # the construct's own on-target log delta
    weak = np.nanmax(kd, axis=1)     # the LESS effective construct (closer to 0 = weaker knockdown)
    P(f"  on-target delta per construct: median {np.nanmedian(kd):+.3f}; the weaker construct's median "
      f"{np.nanmedian(weak):+.3f}; both constructs negative on {float(np.mean((kd < 0).all(1))):.3f} of "
      f"{len(ts)} perturbations")

    gu = np.full((L.NDRAW, len(ts)), np.nan); ce = np.full((L.NDRAW, len(ts)), np.nan)
    for d in range(L.NDRAW):
        r = np.random.default_rng(L.SEED + 11 * d + N)
        G1, G2, C1, C2 = [], [], [], []
        for t in ts:
            big, small = info[t]
            pb = r.permutation(big)
            C1.append(Xd[np.sort(pb[:N])].sum(0, dtype=np.float64)); C2.append(Xd[np.sort(pb[N:2*N])].sum(0, dtype=np.float64))
            G1.append(C1[-1]); G2.append(Xd[np.sort(r.choice(small, N, replace=False))].sum(0, dtype=np.float64))
        gu[d] = L.views(np.stack(G1), np.stack(G2), ctrl, ctrl, scored)["cos_resid"]
        ce[d] = L.views(np.stack(C1), np.stack(C2), ctrl, ctrl, scored)["cos_resid"]
    xg, xc = np.nanmean(gu, 0), np.nanmean(ce, 0)

    out = {"n": len(ts), "kd_median": float(np.nanmedian(kd)), "weak_kd_median": float(np.nanmedian(weak)),
           "strata": {}}
    q = np.nanpercentile(weak, [100/3, 200/3])
    ter = {"T1 weakest KD (least negative)": weak > q[1], "T2": (weak > q[0]) & (weak <= q[1]),
           "T3 strongest KD (most negative)": weak <= q[0]}
    P(f"\n  guide agreement stratified by the WEAKER construct's on-target knockdown (n={N} cells/side):")
    for k, m in ter.items():
        rng = np.random.default_rng(L.SEED)
        mm = m & np.isfinite(xg) & np.isfinite(xc)
        ii = rng.integers(0, int(mm.sum()), size=(L.NBOOT, int(mm.sum())))
        rb = np.median(xg[mm][ii], 1) / np.median(xc[mm][ii], 1)
        out["strata"][k] = dict(n=int(mm.sum()), weak_kd_median=float(np.nanmedian(weak[mm])),
                                guide=dict(**L.stats(xg[mm]), boot=L.boot_median(xg[mm])),
                                cells=dict(**L.stats(xc[mm]), boot=L.boot_median(xc[mm])),
                                diff=L.boot_median_diff(xc[mm], xg[mm]),
                                rho=dict(point=float(np.nanmedian(xg[mm]) / np.nanmedian(xc[mm])),
                                         lo=float(np.percentile(rb, 2.5)), hi=float(np.percentile(rb, 97.5))))
        v = out["strata"][k]
        P(f"    {k:34s} n {v['n']:3d}  weaker-KD median {v['weak_kd_median']:+.3f}  GUIDE "
          f"{L.fmt(v['guide']['boot'])}  CELLS {L.fmt(v['cells']['boot'])}  rho {v['rho']['point']:+.3f} "
          f"[{v['rho']['lo']:+.3f}, {v['rho']['hi']:+.3f}]")
    from scipy.stats import spearmanr
    ok = np.isfinite(xg) & np.isfinite(weak)
    out["spearman_guide_vs_weakKD"] = float(spearmanr(xg[ok], -weak[ok]).statistic)
    out["spearman_cells_vs_weakKD"] = float(spearmanr(xc[ok], -weak[ok]).statistic)
    P(f"  Spearman(guide agreement, knockdown strength) {out['spearman_guide_vs_weakKD']:+.3f}; "
      f"(cell-split agreement, knockdown strength) {out['spearman_cells_vs_weakKD']:+.3f}")
    out["seconds"] = round(time.time() - t0, 1); out["peak_rss_gib"] = L.peak_gib()
    np.savez_compressed(L.WORK / "k562_kd.npz", targets=np.asarray(ts), kd=kd, guide=xg, cells=xc)
    L.Logger.dump(P, out, L.WORK / "k562_kd.json")
    P(f"\ndone {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB"); P.save()


if __name__ == "__main__":
    main()
