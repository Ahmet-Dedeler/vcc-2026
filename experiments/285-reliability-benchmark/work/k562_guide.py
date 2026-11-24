"""285 -- the K562-GWPS GUIDE-CONSTRUCT row on a LARGE cohort.

Only 8 of the 300 panel targets carry two guide constructs in K562-GWPS, which is too few for an interval.
The file carries 775 two-construct perturbations in total (non-panel, mostly essential genes), and the guide
unit is a property of the CRISPRi library, not of the panel -- so the contrast is measured there:
  (b) construct 1's cells vs construct 2's cells, n per side
  (a) two disjoint halves of construct 1's cells, n per side  [the matched cell-sampling reference]
Both readouts; the centring is the matched-n, matched-depth mean over the same cohort, per side.

    PYTHONDONTWRITEBYTECODE=1 uv run python k562_guide.py
"""
from __future__ import annotations
import sys, time
from collections import defaultdict
from pathlib import Path
import h5py, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib285 as L

PATH = L.ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
P = L.Logger(L.WORK / "k562_guide.log")
BLOCK = 8192


def labels_of(h, path):
    o = h[path]
    if isinstance(o, h5py.Group) and "categories" in o:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in o["categories"][:]])
        return cats[o["codes"][:]]
    v = o["values"][:] if isinstance(o, h5py.Group) else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def main():
    t0 = time.time()
    P(f"285 k562_guide.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    with h5py.File(PATH, "r") as h:
        pert = labels_of(h, "obs/perturbation"); guide = labels_of(h, "obs/guide_id")
        genes = labels_of(h, "var/gene_name")
        n_cells, n_genes = h["X"].shape
    pg = defaultdict(set)
    for t, g in zip(pert, guide):
        pg[t].add(g)
    multi = sorted([t for t in pg if t != "control" and len(pg[t]) == 2])
    keep = np.isin(pert, multi) | (pert == "control")
    idx = np.flatnonzero(keep)
    P(f"  {len(multi)} two-construct perturbations; extracting {idx.size:,} cells "
      f"({int((pert == 'control').sum()):,} controls)")
    with h5py.File(PATH, "r") as h:
        X = h["X"]; blocks = []
        for r0 in range(0, n_cells, BLOCK):
            r1 = min(r0 + BLOCK, n_cells)
            m = keep[r0:r1]
            if m.any():
                blocks.append(X[r0:r1][m].astype(np.uint16))
        Xd = np.concatenate(blocks); del blocks
    P(f"  extracted {Xd.shape} uint16, {time.time()-t0:.0f}s, peak {L.peak_gib():.2f} GiB")
    pert, guide = pert[idx], guide[idx]
    is_ctrl = pert == "control"
    ctrl = Xd[is_ctrl].sum(0, dtype=np.float64)
    gpos = {g: i for i, g in enumerate(genes)}
    scored = np.asarray(sorted(set(range(n_genes)) - {gpos[t] for t in multi if t in gpos}), dtype=np.int64)
    info = {}
    for t in multi:
        ix = np.flatnonzero(pert == t)
        cs = sorted([ix[guide[ix] == g] for g in sorted(set(guide[ix].tolist()))], key=lambda a: -a.size)
        info[t] = (cs[0], cs[1])
    out = {"n_multi": len(multi), "rows": {}}
    arrays = {}

    def sums(sel):
        return Xd[np.sort(sel)].sum(0, dtype=np.float64)

    for n in (25, 50):
        ts = [t for t in multi if info[t][0].size >= 2 * n and info[t][1].size >= n]
        if len(ts) < 15:
            continue
        gu = np.full((L.NDRAW, len(ts)), np.nan); ce = np.full((L.NDRAW, len(ts)), np.nan)
        gp_ = np.full((L.NDRAW, len(ts)), np.nan); cp = np.full((L.NDRAW, len(ts)), np.nan)
        agg = {"guide": [], "cells": []}
        for d in range(L.NDRAW):
            r = np.random.default_rng(L.SEED + 11 * d + n)
            G1, G2, C1, C2 = [], [], [], []
            for t in ts:
                big, small = info[t]
                pb = r.permutation(big)
                C1.append(sums(pb[:n])); C2.append(sums(pb[n:2 * n]))
                G1.append(sums(pb[:n])); G2.append(sums(r.choice(small, n, replace=False)))
            vg = L.views(np.stack(G1), np.stack(G2), ctrl, ctrl, scored)
            vc = L.views(np.stack(C1), np.stack(C2), ctrl, ctrl, scored)
            gu[d], ce[d] = vg["cos_resid"], vc["cos_resid"]
            gp_[d], cp[d] = vg["pds"], vc["pds"]
            agg["guide"].append(vg["cos_agg"]); agg["cells"].append(vc["cos_agg"])
        xg, xc = np.nanmean(gu, 0), np.nanmean(ce, 0)
        arrays[f"guide|{n}"] = xg; arrays[f"cells|{n}"] = xc
        rng = np.random.default_rng(L.SEED)
        m = np.isfinite(xg) & np.isfinite(xc)
        ii = rng.integers(0, int(m.sum()), size=(L.NBOOT, int(m.sum())))
        rb = np.median(xg[m][ii], 1) / np.median(xc[m][ii], 1)
        out["rows"][n] = dict(n_targets=len(ts), guide=dict(**L.stats(xg), boot=L.boot_median(xg)),
                              cells=dict(**L.stats(xc), boot=L.boot_median(xc)),
                              cells_minus_guide=L.boot_median_diff(xc, xg),
                              rho=dict(point=float(np.nanmedian(xg) / np.nanmedian(xc)),
                                       lo=float(np.percentile(rb, 2.5)), hi=float(np.percentile(rb, 97.5))),
                              pds_guide=float(np.nanmean(gp_)), pds_cells=float(np.nanmean(cp)),
                              cos_agg_guide=float(np.mean(agg["guide"])), cos_agg_cells=float(np.mean(agg["cells"])))
        v = out["rows"][n]
        P(f"  n={n}/side, {len(ts)} two-construct perturbations: GUIDE {L.fmt(v['guide']['boot'])} "
          f"(<=0 {v['guide']['share_le0']:.3f})  CELLS in one construct {L.fmt(v['cells']['boot'])} "
          f"(<=0 {v['cells']['share_le0']:.3f})  CELLS − GUIDE {L.fmt(v['cells_minus_guide'])}  "
          f"rho {v['rho']['point']:+.3f} [{v['rho']['lo']:+.3f}, {v['rho']['hi']:+.3f}]")
        P(f"      scorer's statistic: mean pds GUIDE {v['pds_guide']:.4f} vs CELLS {v['pds_cells']:.4f}; "
          f"cos_agg {v['cos_agg_guide']:.4f} vs {v['cos_agg_cells']:.4f}")
    out["seconds"] = round(time.time() - t0, 1); out["peak_rss_gib"] = L.peak_gib()
    np.savez_compressed(L.WORK / "k562_guide.npz", targets=np.asarray(multi), **arrays)
    L.Logger.dump(P, out, L.WORK / "k562_guide.json")
    P(f"\ndone {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB"); P.save()


if __name__ == "__main__":
    main()
