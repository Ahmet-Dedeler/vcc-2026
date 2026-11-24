"""285 -- the Nadig 2025 rows (Jurkat = the challenge's line A tissue, and HepG2): (a) the cell ladder,
(b) the guide construct, (d) the gemgroup, on the 2,393 DepMap-essential targets (0 of the 300 panel).

X is a DENSE UNCOMPRESSED float32 matrix, so the cells of the chosen targets are read directly.

    PYTHONDONTWRITEBYTECODE=1 uv run python nadig_285.py jurkat|hepg2
"""
from __future__ import annotations
import sys, time
from collections import Counter, defaultdict
from pathlib import Path
import h5py, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib285 as L

WHICH = sys.argv[1] if len(sys.argv) > 1 else "jurkat"
PATH = L.ROOT / f"data/external/nadig2025/{WHICH}_raw_singlecell.h5ad"
P = L.Logger(L.WORK / f"nadig_{WHICH}.log")
RUN = (25, 50)


def col(h, frame, key):
    g = h[frame]
    if "__categories" in g and key in g["__categories"]:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in g["__categories"][key][:]])
        return cats[g[key][:]]
    v = g[key][:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def main():
    t0 = time.time()
    P(f"285 nadig_285.py [{WHICH}] -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    with h5py.File(PATH, "r") as h:
        gene = col(h, "obs", "gene"); sg = col(h, "obs", "sgID_AB"); gg = col(h, "obs", "gem_group")
        gnames = col(h, "var", "gene_name")
        n_cells, n_genes = h["X"].shape
    cg = Counter(gene.tolist())
    ctrl_lab = [t for t in cg if t.lower().startswith("non-target")]
    real = sorted(set(cg) - set(ctrl_lab))
    pgu = defaultdict(set)
    for t, g in zip(gene, sg):
        pgu[t].add(g)
    multi = sorted([t for t in real if len(pgu[t]) > 1 and cg[t] >= 60])
    rng = np.random.default_rng(L.SEED)
    others = sorted(rng.choice([t for t in real if t not in set(multi) and cg[t] >= 60],
                               size=min(400, sum(1 for t in real if t not in set(multi) and cg[t] >= 60)),
                               replace=False).tolist())
    use = sorted(set(multi) | set(others))
    keep = np.isin(gene, use) | np.isin(gene, ctrl_lab)
    idx = np.flatnonzero(keep)
    P(f"  {n_cells:,} cells x {n_genes} genes; {len(real)} targets, {len(multi)} multi-guide (>=60 cells), "
      f"{len(others)} single-guide sampled; extracting {idx.size:,} cells")
    with h5py.File(PATH, "r") as h:
        X = h["X"]
        blocks = []
        for r0 in range(0, n_cells, 16384):
            r1 = min(r0 + 16384, n_cells)
            m = keep[r0:r1]
            if m.any():
                blocks.append(np.asarray(X[r0:r1][m], dtype=np.float32))
        Xd = np.concatenate(blocks); del blocks
    P(f"  extracted {Xd.shape}, {time.time()-t0:.0f}s, peak {L.peak_gib():.2f} GiB")
    gene, sg, gg = gene[idx], sg[idx], gg[idx]
    is_ctrl = np.isin(gene, ctrl_lab)
    ctrl_full = Xd[is_ctrl].sum(0, dtype=np.float64)
    gpos = {g: i for i, g in enumerate(gnames)}
    scored = np.asarray(sorted(set(range(n_genes)) - {gpos[t] for t in use if t in gpos}), dtype=np.int64)
    idx_t = {t: np.flatnonzero(gene == t) for t in use}
    gems = sorted(set(gg.tolist())); gpar = np.asarray([gems.index(x) % 2 for x in gg])
    out = {"n_targets": len(use), "n_multi": len(multi), "n_gem": len(gems), "budgets": {}}
    arrays = {}
    P(f"  {len(use)} targets used, {scored.size} scored columns, {len(gems)} gem groups")

    def sums(sel):
        return Xd[np.sort(sel)].sum(0, dtype=np.float64)

    for n in RUN:
        lad = np.full((L.NDRAW, len(use)), np.nan); gemr = np.full((L.NDRAW, len(use)), np.nan)
        gud = np.full((L.NDRAW, len(multi)), np.nan); gcl = np.full((L.NDRAW, len(multi)), np.nan)
        for d in range(L.NDRAW):
            r = np.random.default_rng(L.SEED + 11 * d + n)
            keepj, SA, SB = [], [], []
            for j, t in enumerate(use):
                ix = idx_t[t]
                if ix.size < 2 * n:
                    continue
                p = r.permutation(ix)
                keepj.append(j); SA.append(sums(p[:n])); SB.append(sums(p[n:2 * n]))
            SA, SB = np.stack(SA), np.stack(SB)
            PA, PB = L.prof(SA)[:, scored], L.prof(SB)[:, scored]
            pmA, pmB = PA.mean(0), PB.mean(0)
            lad[d, np.asarray(keepj)] = L.row_cos(PA - pmA, PB - pmB)
            k2, GA, GB = [], [], []
            for j, t in enumerate(use):
                ix = idx_t[t]
                a, b = ix[gpar[ix] == 0], ix[gpar[ix] == 1]
                if a.size < n or b.size < n:
                    continue
                k2.append(j); GA.append(sums(r.choice(a, n, replace=False))); GB.append(sums(r.choice(b, n, replace=False)))
            if len(k2) >= 15:
                QA, QB = L.prof(np.stack(GA))[:, scored], L.prof(np.stack(GB))[:, scored]
                gemr[d, np.asarray(k2)] = L.row_cos(QA - pmA, QB - pmB)
            k3, H1, H2, C1, C2 = [], [], [], [], []
            for j, t in enumerate(multi):
                ix = idx_t[t]
                cs = sorted([ix[sg[ix] == g] for g in sorted(set(sg[ix].tolist()))], key=lambda a: -a.size)
                if cs[0].size < 2 * n or cs[1].size < n:
                    continue
                pb = r.permutation(cs[0])
                k3.append(j); H1.append(sums(pb[:n])); H2.append(sums(r.choice(cs[1], n, replace=False)))
                C1.append(sums(pb[:n])); C2.append(sums(pb[n:2 * n]))
            if len(k3) >= 15:
                RA, RB = L.prof(np.stack(H1))[:, scored], L.prof(np.stack(H2))[:, scored]
                gud[d, np.asarray(k3)] = L.row_cos(RA - pmA, RB - pmB)
                UA, UB = L.prof(np.stack(C1))[:, scored], L.prof(np.stack(C2))[:, scored]
                gcl[d, np.asarray(k3)] = L.row_cos(UA - pmA, UB - pmB)
        xa, xg, xq, xc = (np.nanmean(z, 0) for z in (lad, gemr, gud, gcl))
        arrays[f"a|{n}"] = xa; arrays[f"gem|{n}"] = xg; arrays[f"guide|{n}"] = xq; arrays[f"guide_cellctrl|{n}"] = xc
        out["budgets"][n] = dict(cells=dict(**L.stats(xa), boot=L.boot_median(xa)),
                                 gem=dict(**L.stats(xg), boot=L.boot_median(xg)),
                                 guide=dict(**L.stats(xq), boot=L.boot_median(xq)),
                                 guide_cellctrl=dict(**L.stats(xc), boot=L.boot_median(xc)),
                                 cells_minus_gem=L.boot_median_diff(xa, xg),
                                 cells_minus_guide=L.boot_median_diff(xc, xq))
        P(f"  n={n:3d}/side: CELLS {L.fmt(out['budgets'][n]['cells']['boot'])} (n_t {out['budgets'][n]['cells']['n']}) | "
          f"GEMGROUP {L.fmt(out['budgets'][n]['gem']['boot'])} (CELLS − GEM {L.fmt(out['budgets'][n]['cells_minus_gem'])}) | "
          f"GUIDE {L.fmt(out['budgets'][n]['guide']['boot'])} vs same-construct CELLS "
          f"{L.fmt(out['budgets'][n]['guide_cellctrl']['boot'])} (CELLS − GUIDE {L.fmt(out['budgets'][n]['cells_minus_guide'])}, "
          f"n_t {out['budgets'][n]['guide']['n']})")
    out["seconds"] = round(time.time() - t0, 1); out["peak_rss_gib"] = L.peak_gib()
    np.savez_compressed(L.WORK / f"nadig_{WHICH}.npz", targets=np.asarray(use), multi=np.asarray(multi), **arrays)
    L.Logger.dump(P, out, L.WORK / f"nadig_{WHICH}.json")
    P(f"\ndone {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB"); P.save()


if __name__ == "__main__":
    main()
