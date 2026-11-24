"""296 -- K562-GWPS (CRISPRi, Replogle 2022): the knockdown-stratified construct term.

One streaming pass gives (i) the capped per-construct cells for the matched-n ladder and (ii) EXACT
full-depth on-target statistics (the construct's own gene column and its total UMIs) for EVERY perturbation
in the file -- so the stratifier and the library-wide failure rate come from all cells, while the ladder
uses a capped, memory-bounded subset.

Per two-construct perturbation, at n cells/side and 8 draws (285's design, 293's contrast set):
  (a)  the LARGER construct's two disjoint cell halves      (the matched cell-sampling reference)
  (a2) the SMALLER construct's own two halves               (its own reliability)
  (b)  larger construct's side A vs the smaller construct   (the construct term)
  (x)  side A vs the NEXT perturbation's smaller construct  (the cross-target null)
plus a per-CONSTRUCT response matrix (one n-cell draw each) for the failing-construct phenomenon.

    PYTHONDONTWRITEBYTECODE=1 uv run python k562_296.py
"""
from __future__ import annotations

import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib296 as X  # noqa: E402
from lib296 import L  # noqa: E402

LIB = "K562-GWPS"
PATH = X.ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
P = L.Logger(X.WORK / "k562_296.log")
BLOCK = 8192
CONS_CAP = 300
CTRL_CAP = 25000
RUNGS = (25, 50)


def labels_of(h, path):
    o = h[path]
    if isinstance(o, h5py.Group) and "categories" in o:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in o["categories"][:]])
        return cats[o["codes"][:]]
    v = o["values"][:] if isinstance(o, h5py.Group) else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def main() -> None:
    t0 = time.time()
    P(f"296 k562_296.py [{LIB}] -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    with h5py.File(PATH, "r") as h:
        pert = labels_of(h, "obs/perturbation")
        guide = labels_of(h, "obs/guide_id")
        genes = labels_of(h, "var/gene_name")
        n_cells, n_genes = h["X"].shape
    gpos = {g: i for i, g in enumerate(genes)}
    pg = defaultdict(set)
    for t, g in zip(pert, guide):
        pg[t].add(g)
    is_ctrl = pert == "control"
    multi = sorted([t for t in pg if t != "control" and len(pg[t]) == 2])
    single = sorted([t for t in pg if t != "control" and len(pg[t]) == 1])
    P(f"  {n_cells:,} cells x {n_genes} genes; {len(pg) - 1} perturbations "
      f"({len(multi)} two-construct, {len(single)} single-construct), {int(is_ctrl.sum()):,} control cells")

    # ------------------------------------------------------------------ what to extract (capped)
    rng0 = np.random.default_rng(L.SEED + 5)
    cons_rows: dict[str, np.ndarray] = {}
    want = np.zeros(n_cells, bool)
    for t in multi:
        for g in sorted(pg[t]):
            ix = np.flatnonzero((pert == t) & (guide == g))
            if ix.size > CONS_CAP:
                ix = np.sort(rng0.choice(ix, CONS_CAP, replace=False))
            cons_rows[g] = ix
            want[ix] = True
    cidx = np.flatnonzero(is_ctrl)
    if cidx.size > CTRL_CAP:
        cidx = np.sort(rng0.choice(cidx, CTRL_CAP, replace=False))
    want[cidx] = True
    P(f"  extracting {int(want.sum()):,} cells (construct cap {CONS_CAP}, control cap {CTRL_CAP}) "
      f"-> {int(want.sum()) * n_genes * 2 / 1024**3:.2f} GiB uint16")

    # ------------------------------------------------------------------ full-depth on-target bookkeeping
    # group = 0..M-1 the two-construct CONSTRUCTS, M.. the single-construct TARGETS, last = controls
    allc = [g for t in multi for g in sorted(pg[t])]
    cg = {g: i for i, g in enumerate(allc)}
    M = len(allc)
    grp = np.full(n_cells, -1, np.int64)
    for t in multi:
        for g in sorted(pg[t]):
            grp[(pert == t) & (guide == g)] = cg[g]
    spos = {t: M + i for i, t in enumerate(single)}
    for t in single:
        grp[pert == t] = spos[t]
    G_CTRL = M + len(single)
    grp[is_ctrl] = G_CTRL
    # the gene column whose counts each group's on-target statistic needs (-1 = gene not on the axis)
    grp_col = np.full(G_CTRL + 1, -1, np.int64)
    for t in multi:
        for g in sorted(pg[t]):
            grp_col[cg[g]] = gpos.get(t, -1)
    for t in single:
        grp_col[spos[t]] = gpos.get(t, -1)
    on = np.zeros(G_CTRL + 1)          # summed counts at the group's own gene column
    tot = np.zeros(G_CTRL + 1)         # summed total counts
    ncl = np.zeros(G_CTRL + 1, np.int64)
    ctrl_sum = np.zeros(n_genes)

    out = np.zeros((int(want.sum()), n_genes), np.uint16)
    dst = np.cumsum(want) - 1
    with h5py.File(PATH, "r") as h:
        Xh = h["X"]
        for r0 in range(0, n_cells, BLOCK):
            r1 = min(r0 + BLOCK, n_cells)
            gb = grp[r0:r1]
            wb = want[r0:r1]
            if not ((gb >= 0).any() or wb.any()):
                continue
            B = np.asarray(Xh[r0:r1], dtype=np.float32)
            k = gb >= 0
            if k.any():
                gk = gb[k]
                np.add.at(tot, gk, B[k].sum(1, dtype=np.float64))
                np.add.at(ncl, gk, 1)
                col = grp_col[gk]
                ok = col >= 0
                if ok.any():
                    np.add.at(on, gk[ok], B[k][ok, col[ok]].astype(np.float64))
            cb = is_ctrl[r0:r1]
            if cb.any():
                ctrl_sum += B[cb].sum(0, dtype=np.float64)
            if wb.any():
                D = B[wb]
                assert D.max(initial=0.0) < 65535.0 and np.array_equal(D, np.rint(D)), "count overflow"
                out[dst[r0:r1][wb]] = D.astype(np.uint16)
            del B
    Xd = out
    P(f"  extracted {Xd.shape} uint16, {time.time() - t0:.0f}s, peak {L.peak_gib():.2f} GiB")

    # ------------------------------------------------------------------ full-depth knockdown per construct
    pc = L.prof(ctrl_sum)
    ctrl_p = ctrl_sum / ctrl_sum.sum()

    def kd_of(gi: int) -> tuple[float, float]:
        """(on-target log delta, knockdown fraction) of one group, at FULL depth."""
        j = grp_col[gi]
        if j < 0 or tot[gi] <= 0:
            return np.nan, np.nan
        p = on[gi] / tot[gi]
        return float(np.log1p(L.TS * p) - pc[j]), float(1.0 - p / ctrl_p[j]) if ctrl_p[j] > 0 else np.nan

    kd_cons = {g: kd_of(cg[g]) for g in allc}
    kd_single = {t: kd_of(spos[t]) for t in single}
    P(f"  full depth: {int(ncl[G_CTRL]):,} control cells; on-target lfc measurable on "
      f"{sum(np.isfinite(v[0]) for v in kd_cons.values())}/{len(allc)} constructs and "
      f"{sum(np.isfinite(v[0]) for v in kd_single.values())}/{len(single)} single-construct targets")

    # ------------------------------------------------------------------ the matched ladder
    ctrl_cells = np.flatnonzero(is_ctrl[want])          # positions of the kept control cells inside Xd
    ctrl_bank = Xd[ctrl_cells].sum(0, dtype=np.float64)
    scored = np.asarray(sorted(set(range(n_genes)) - {gpos[t] for t in multi if t in gpos}), np.int64)
    pos = {g: dst[cons_rows[g]] for g in allc}          # construct -> row positions inside Xd

    def sums(sel):
        return Xd[np.sort(sel)].sum(0, dtype=np.float64)

    res: dict = {"library": LIB, "n_multi": len(multi), "n_single": len(single),
                 "n_control_cells": int(is_ctrl.sum()), "cons_cap": CONS_CAP, "ctrl_cap": CTRL_CAP,
                 "scored_columns": int(scored.size), "rungs": {},
                 "kd_cons": {g: dict(lfc=kd_cons[g][0], kfrac=kd_cons[g][1], cells=int(ncl[cg[g]]))
                             for g in allc},
                 "kd_single": {t: dict(lfc=kd_single[t][0], kfrac=kd_single[t][1], cells=int(ncl[spos[t]]))
                               for t in single}}
    store: dict[str, np.ndarray] = {}

    for n in RUNGS:
        cohort = []
        for t in multi:
            g1, g2 = sorted(pg[t])
            big, small = (g1, g2) if pos[g1].size >= pos[g2].size else (g2, g1)
            if pos[big].size >= 2 * n and pos[small].size >= n:
                cohort.append((t, big, small))
        if len(cohort) < 15:
            continue
        K = len(cohort)
        keys = ("a", "a2", "b", "x")
        cosr = {k: np.full((X.NDRAW, K), np.nan) for k in keys}
        pdsr = {k: np.full((X.NDRAW, K), np.nan) for k in keys}
        # per-construct rows (one n-cell draw each) for the failing-construct phenomenon
        clist = [g for _, big, small in cohort for g in (big, small)]
        cmean = np.full((X.NDRAW, len(clist)), np.nan)     # cos to the library's mean response (uncentred)
        cross_r: list[np.ndarray] = []                     # residual-space pairwise cosines, per draw
        cross_u: list[np.ndarray] = []                     # uncentred-space pairwise cosines, per draw
        for d in range(X.NDRAW):
            r = np.random.default_rng(L.SEED + 11 * d + n)
            SA, SB = {}, {}
            for t, big, small in cohort:
                pb = r.permutation(pos[big])
                A = sums(pb[:n])
                SA["a"] = SA.get("a", []) + [A]
                SB["a"] = SB.get("a", []) + [sums(pb[n:2 * n])]
                SA["b"] = SA.get("b", []) + [A]
                SB["b"] = SB.get("b", []) + [sums(r.choice(pos[small], n, replace=False))]
                if pos[small].size >= 2 * n:
                    ps = r.permutation(pos[small])
                    SA["a2"] = SA.get("a2", []) + [sums(ps[:n])]
                    SB["a2"] = SB.get("a2", []) + [sums(ps[n:2 * n])]
                else:
                    SA["a2"] = SA.get("a2", []) + [None]
                    SB["a2"] = SB.get("a2", []) + [None]
            # (x): side A of t against the NEXT perturbation's smaller construct
            SA["x"] = list(SA["a"])
            SB["x"] = [SB["b"][(i + 1) % K] for i in range(K)]
            for k in keys:
                ok = np.asarray([v is not None for v in SA[k]])
                if ok.sum() < 15:
                    continue
                MA = np.stack([SA[k][i] if ok[i] else SA["a"][i] for i in range(K)])
                MB = np.stack([SB[k][i] if ok[i] else SB["a"][i] for i in range(K)])
                v = L.views(MA, MB, ctrl_bank, ctrl_bank, scored)
                cosr[k][d, ok] = v["cos_resid"][ok]
                pdsr[k][d, ok] = v["pds"][ok]
            # per-construct matrix
            CM = np.stack([sums(r.choice(pos[g], n, replace=False)) for g in clist])
            DC = L.delta(CM, ctrl_bank)[:, scored]
            mu = DC.mean(0)
            nrm = np.linalg.norm(DC, axis=1)
            with np.errstate(invalid="ignore", divide="ignore"):
                cmean[d] = (DC @ mu) / (nrm * np.linalg.norm(mu))
            RC = DC - mu[None, :]
            for Mx, acc in ((RC, cross_r), (DC, cross_u)):
                nn = np.linalg.norm(Mx, axis=1)
                with np.errstate(invalid="ignore", divide="ignore"):
                    acc.append((Mx @ Mx.T) / np.outer(nn, nn))
            del CM, DC, RC
        m = {k: np.nanmean(cosr[k], 0) for k in keys}
        mp = {f"pds_{k}": np.nanmean(pdsr[k], 0) for k in keys}
        row = dict(n=n, k=K, targets=[t for t, _, _ in cohort],
                   c_big=[b for _, b, _ in cohort], c_small=[s for _, _, s in cohort],
                   cells_big=[int(ncl[cg[b]]) for _, b, _ in cohort],
                   cells_small=[int(ncl[cg[s]]) for _, _, s in cohort])
        for k in keys:
            row[k] = dict(**L.stats(m[k]), boot=L.boot_median(m[k]))
        row["rho_guide"] = X.rho_boot(m["b"], m["a"])
        row["a_minus_b"] = L.boot_median_diff(m["a"], m["b"])
        row["b_minus_x"] = L.boot_median_diff(m["b"], m["x"])
        row["pds"] = {k: float(np.nanmean(mp[f"pds_{k}"])) for k in keys}
        res["rungs"][n] = row
        for k in keys:
            store[f"{n}|{k}"] = m[k]
            store[f"{n}|pds_{k}"] = mp[f"pds_{k}"]
        store[f"{n}|targets"] = np.asarray([t for t, _, _ in cohort])
        store[f"{n}|c_big"] = np.asarray([b for _, b, _ in cohort])
        store[f"{n}|c_small"] = np.asarray([s for _, _, s in cohort])
        store[f"{n}|constructs"] = np.asarray(clist)
        store[f"{n}|cos_to_mean"] = np.nanmean(cmean, 0)
        store[f"{n}|cross_resid"] = np.nanmean(np.stack(cross_r), 0).astype(np.float32)
        store[f"{n}|cross_uncent"] = np.nanmean(np.stack(cross_u), 0).astype(np.float32)
        P(f"  n={n}/side, {K} pairs:  (a) {L.fmt(row['a']['boot'])}  (a2) {L.fmt(row['a2']['boot'])}  "
          f"(b) {L.fmt(row['b']['boot'])}  (x) {L.fmt(row['x']['boot'])}  rho_guide {X.fmt_rho(row['rho_guide'])}")
        P(f"           pds (a) {row['pds']['a']:.4f} / (b) {row['pds']['b']:.4f} / (x) {row['pds']['x']:.4f}")
        del cross_r, cross_u

    res["seconds"] = round(time.time() - t0, 1)
    res["peak_rss_gib"] = L.peak_gib()
    np.savez_compressed(X.WORK / "k562_296.npz", **store)
    (X.WORK / "k562_296.json").write_text(json.dumps(res, indent=1, default=float))
    P(f"\ndone {res['seconds']}s, peak {res['peak_rss_gib']:.2f} GiB")
    P.save()


if __name__ == "__main__":
    main()
