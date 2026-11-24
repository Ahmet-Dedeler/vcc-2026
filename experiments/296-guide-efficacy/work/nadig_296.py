"""296 -- Nadig 2025 (CRISPRi, Jurkat / HepG2): the knockdown-stratified construct term.

Same contrast set as 293 and as `k562_296.py`; the centring / distractor POOL follows 285's `nadig_285.py`
(the two-construct targets plus up to 400 sampled single-construct targets), so the ALL row is directly
comparable with 285's published rho_guide for this library.

    PYTHONDONTWRITEBYTECODE=1 uv run python nadig_296.py jurkat|hepg2
"""
from __future__ import annotations

import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib296 as X  # noqa: E402
from lib296 import L  # noqa: E402

WHICH = sys.argv[1] if len(sys.argv) > 1 else "jurkat"
PATH = X.ROOT / f"data/external/nadig2025/{WHICH}_raw_singlecell.h5ad"
P = L.Logger(X.WORK / f"nadig_{WHICH}_296.log")
BLOCK = 16384
CONS_CAP = 300
POOL_CAP = 300
CTRL_CAP = 20000
RUNGS = (25, 50)
MIN_COHORT = 15


def col(h, frame, key):
    g = h[frame]
    if "__categories" in g and key in g["__categories"]:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in g["__categories"][key][:]])
        return cats[g[key][:]]
    v = g[key][:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def main() -> None:
    t0 = time.time()
    LIB = f"Nadig-{WHICH}"
    P(f"296 nadig_296.py [{LIB}] -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    with h5py.File(PATH, "r") as h:
        gene = col(h, "obs", "gene")
        sg = col(h, "obs", "sgID_AB")
        gnames = col(h, "var", "gene_name")
        n_cells, n_genes = h["X"].shape
    gpos = {g: i for i, g in enumerate(gnames)}
    cnt = Counter(gene.tolist())
    ctrl_lab = [t for t in cnt if t.lower().startswith("non-target")]
    is_ctrl = np.isin(gene, ctrl_lab)
    real = sorted(set(cnt) - set(ctrl_lab))
    pgu = defaultdict(set)
    for t, g in zip(gene, sg):
        pgu[t].add(g)
    multi = sorted([t for t in real if len(pgu[t]) > 1 and cnt[t] >= 60])
    singles = sorted([t for t in real if len(pgu[t]) == 1])
    rng0 = np.random.default_rng(L.SEED)
    others = sorted(rng0.choice([t for t in real if t not in set(multi) and cnt[t] >= 60],
                                size=min(400, sum(1 for t in real if t not in set(multi) and cnt[t] >= 60)),
                                replace=False).tolist())
    pool = sorted(set(multi) | set(others))
    P(f"  {n_cells:,} cells x {n_genes} genes; {len(real)} targets ({len(multi)} multi-construct >=60 cells, "
      f"{len(singles)} single-construct), {int(is_ctrl.sum()):,} control cells; pool {len(pool)}")

    # ------------------------------------------------------------------ what to extract (capped)
    idx_t = {t: np.flatnonzero(gene == t) for t in pool}
    cons_rows: dict[str, np.ndarray] = {}
    want = np.zeros(n_cells, bool)
    for t in multi:
        for g in sorted(pgu[t]):
            ix = idx_t[t][sg[idx_t[t]] == g]
            if ix.size > CONS_CAP:
                ix = np.sort(rng0.choice(ix, CONS_CAP, replace=False))
            cons_rows[f"{t}||{g}"] = ix
            want[ix] = True
    pool_rows: dict[str, np.ndarray] = {}
    for t in pool:
        ix = idx_t[t]
        if ix.size > POOL_CAP:
            ix = np.sort(rng0.choice(ix, POOL_CAP, replace=False))
        pool_rows[t] = ix
        want[ix] = True
    cidx = np.flatnonzero(is_ctrl)
    if cidx.size > CTRL_CAP:
        cidx = np.sort(rng0.choice(cidx, CTRL_CAP, replace=False))
    want[cidx] = True

    # ------------------------------------------------------------------ full-depth on-target bookkeeping
    allc = [f"{t}||{g}" for t in multi for g in sorted(pgu[t])]
    cg = {c: i for i, c in enumerate(allc)}
    M = len(allc)
    grp = np.full(n_cells, -1, np.int64)
    for c in allc:
        t, g = c.split("||")
        grp[(gene == t) & (sg == g)] = cg[c]
    # a SECOND group vector for the pooled per-target statistic (all cells of a target, constructs included)
    tpos = {t: M + i for i, t in enumerate(real)}
    tg = np.full(n_cells, -1, np.int64)
    for t in real:
        tg[gene == t] = tpos[t]
    G = M + len(real)
    grp_col = np.full(G, -1, np.int64)
    for c in allc:
        grp_col[cg[c]] = gpos.get(c.split("||")[0], -1)
    for t in real:
        grp_col[tpos[t]] = gpos.get(t, -1)
    on = np.zeros(G); tot = np.zeros(G); ncl = np.zeros(G, np.int64)
    ctrl_sum = np.zeros(n_genes)

    out = np.zeros((int(want.sum()), n_genes), np.uint16)
    dst = np.cumsum(want) - 1
    P(f"  extracting {int(want.sum()):,} cells -> {int(want.sum()) * n_genes * 2 / 1024**3:.2f} GiB uint16")
    with h5py.File(PATH, "r") as h:
        Xh = h["X"]
        for r0 in range(0, n_cells, BLOCK):
            r1 = min(r0 + BLOCK, n_cells)
            B = np.asarray(Xh[r0:r1], dtype=np.float32)
            rs = B.sum(1, dtype=np.float64)
            for gb in (grp[r0:r1], tg[r0:r1]):
                k = gb >= 0
                if not k.any():
                    continue
                gk = gb[k]
                np.add.at(tot, gk, rs[k])
                np.add.at(ncl, gk, 1)
                cc = grp_col[gk]
                ok = cc >= 0
                if ok.any():
                    np.add.at(on, gk[ok], B[k][ok, cc[ok]].astype(np.float64))
            cb = is_ctrl[r0:r1]
            if cb.any():
                ctrl_sum += B[cb].sum(0, dtype=np.float64)
            wb = want[r0:r1]
            if wb.any():
                D = B[wb]
                assert D.max(initial=0.0) < 65535.0 and np.array_equal(D, np.rint(D)), "count overflow"
                out[dst[r0:r1][wb]] = D.astype(np.uint16)
            del B
    Xd = out
    P(f"  extracted {Xd.shape} uint16, {time.time() - t0:.0f}s, peak {L.peak_gib():.2f} GiB")

    pc = L.prof(ctrl_sum)
    ctrl_p = ctrl_sum / ctrl_sum.sum()

    def kd_of(gi: int):
        j = grp_col[gi]
        if j < 0 or tot[gi] <= 0:
            return np.nan, np.nan
        p = on[gi] / tot[gi]
        return float(np.log1p(L.TS * p) - pc[j]), (float(1.0 - p / ctrl_p[j]) if ctrl_p[j] > 0 else np.nan)

    kd_cons = {c: kd_of(cg[c]) for c in allc}
    kd_target = {t: kd_of(tpos[t]) for t in real}

    # ------------------------------------------------------------------ the ladder
    ctrl_bank = Xd[np.flatnonzero(is_ctrl[want])].sum(0, dtype=np.float64)
    scored = np.asarray(sorted(set(range(n_genes)) - {gpos[t] for t in pool if t in gpos}), np.int64)
    pos = {c: dst[cons_rows[c]] for c in allc}
    ppos = {t: dst[pool_rows[t]] for t in pool}
    P(f"  {len(pool)} pool targets, {scored.size} scored columns, {len(allc)} constructs; "
      f"on-target measurable on {sum(np.isfinite(v[0]) for v in kd_cons.values())}/{len(allc)} constructs, "
      f"{sum(np.isfinite(v[0]) for v in kd_target.values())}/{len(real)} targets")

    def sums(sel):
        return Xd[np.sort(sel)].sum(0, dtype=np.float64)

    res: dict = {"library": LIB, "n_targets": len(real), "n_multi": len(multi), "n_single": len(singles),
                 "n_control_cells": int(is_ctrl.sum()), "pool": len(pool), "scored_columns": int(scored.size),
                 "kd_cons": {c: dict(lfc=kd_cons[c][0], kfrac=kd_cons[c][1], cells=int(ncl[cg[c]]))
                             for c in allc},
                 "kd_target": {t: dict(lfc=kd_target[t][0], kfrac=kd_target[t][1], cells=int(ncl[tpos[t]]),
                                       n_constructs=len(pgu[t])) for t in real},
                 "rungs": {}}
    store: dict[str, np.ndarray] = {}

    for n in RUNGS:
        pool_ok = [t for t in pool if ppos[t].size >= 2 * n]
        cohort = []
        for t in multi:
            cs = sorted([f"{t}||{g}" for g in sorted(pgu[t])], key=lambda c: -pos[c].size)
            big, small = cs[0], cs[1]
            if pos[big].size >= 2 * n and pos[small].size >= n:
                cohort.append((t, big, small))
        if len(cohort) < MIN_COHORT:
            P(f"  n={n}: cohort {len(cohort)} < {MIN_COHORT}, skipped")
            continue
        K = len(cohort)
        keys = ("a", "a2", "b", "x")
        cosr = {k: np.full((X.NDRAW, K), np.nan) for k in keys}
        pdsr = {k: np.full((X.NDRAW, K), np.nan) for k in keys}
        clist = [c for _, big, small in cohort for c in (big, small)]
        cmean = np.full((X.NDRAW, len(clist)), np.nan)
        cross_r, cross_u = [], []
        bi = {t: i for i, t in enumerate(pool_ok)}
        for d in range(X.NDRAW):
            r = np.random.default_rng(L.SEED + 11 * d + n)
            base, baseB = [], []
            for t in pool_ok:
                p = r.permutation(ppos[t])
                base.append(sums(p[:n])); baseB.append(sums(p[n:2 * n]))
            base, baseB = np.stack(base), np.stack(baseB)
            DA = (L.prof(base) - L.prof(ctrl_bank)[None, :])[:, scored]
            DB = (L.prof(baseB) - L.prof(ctrl_bank)[None, :])[:, scored]
            mA, mB = DA.mean(0), DB.mean(0)
            sides = {k: {} for k in keys}
            for t, big, small in cohort:
                pb = r.permutation(pos[big])
                A = sums(pb[:n])
                sides["a"][t] = (A, sums(pb[n:2 * n]))
                sides["b"][t] = (A, sums(r.choice(pos[small], n, replace=False)))
                if pos[small].size >= 2 * n:
                    ps = r.permutation(pos[small])
                    sides["a2"][t] = (sums(ps[:n]), sums(ps[n:2 * n]))
            for i_t, (t, big, small) in enumerate(cohort):
                u, _, su = cohort[(i_t + 1) % K]
                sides["x"][t] = (sides["a"][t][0], sums(r.choice(pos[su], n, replace=False)))
            for k in keys:
                have = [t for t, _, _ in cohort if t in sides[k]]
                if len(have) < MIN_COHORT:
                    continue
                MA, MB, rows = base.copy(), baseB.copy(), []
                for t in have:
                    if t in bi:
                        MA[bi[t]] = sides[k][t][0]; MB[bi[t]] = sides[k][t][1]; rows.append(bi[t])
                    else:
                        MA = np.vstack([MA, sides[k][t][0][None, :]])
                        MB = np.vstack([MB, sides[k][t][1][None, :]])
                        rows.append(MA.shape[0] - 1)
                v = L.views(MA, MB, ctrl_bank, ctrl_bank, scored, mA, mB)
                for jj, t in enumerate(have):
                    kk = [x[0] for x in cohort].index(t)
                    cosr[k][d, kk] = v["cos_resid"][rows[jj]]
                    pdsr[k][d, kk] = v["pds"][rows[jj]]
            CM = np.stack([sums(r.choice(pos[c], n, replace=False)) for c in clist])
            DC = L.delta(CM, ctrl_bank)[:, scored]
            mu = DC.mean(0)
            nn = np.linalg.norm(DC, axis=1)
            with np.errstate(invalid="ignore", divide="ignore"):
                cmean[d] = (DC @ mu) / (nn * np.linalg.norm(mu))
            RC = DC - mu[None, :]
            for Mx, acc in ((RC, cross_r), (DC, cross_u)):
                q = np.linalg.norm(Mx, axis=1)
                with np.errstate(invalid="ignore", divide="ignore"):
                    acc.append((Mx @ Mx.T) / np.outer(q, q))
            del CM, DC, RC, base, baseB
        m = {k: np.nanmean(cosr[k], 0) for k in keys}
        mp = {k: np.nanmean(pdsr[k], 0) for k in keys}
        row = dict(n=n, k=K, targets=[t for t, _, _ in cohort], n_pool=len(pool_ok),
                   c_big=[b for _, b, _ in cohort], c_small=[s for _, _, s in cohort],
                   cells_big=[int(ncl[cg[b]]) for _, b, _ in cohort],
                   cells_small=[int(ncl[cg[s]]) for _, _, s in cohort])
        for k in keys:
            row[k] = dict(**L.stats(m[k]), boot=L.boot_median(m[k]))
        row["rho_guide"] = X.rho_boot(m["b"], m["a"])
        row["a_minus_b"] = L.boot_median_diff(m["a"], m["b"])
        row["b_minus_x"] = L.boot_median_diff(m["b"], m["x"])
        row["pds"] = {k: float(np.nanmean(mp[k])) for k in keys}
        res["rungs"][n] = row
        for k in keys:
            store[f"{n}|{k}"] = m[k]; store[f"{n}|pds_{k}"] = mp[k]
        store[f"{n}|targets"] = np.asarray([t for t, _, _ in cohort])
        store[f"{n}|c_big"] = np.asarray([b for _, b, _ in cohort])
        store[f"{n}|c_small"] = np.asarray([s for _, _, s in cohort])
        store[f"{n}|constructs"] = np.asarray(clist)
        store[f"{n}|cos_to_mean"] = np.nanmean(cmean, 0)
        store[f"{n}|cross_resid"] = np.nanmean(np.stack(cross_r), 0).astype(np.float32)
        store[f"{n}|cross_uncent"] = np.nanmean(np.stack(cross_u), 0).astype(np.float32)
        P(f"  n={n}/side, {K} pairs (pool {len(pool_ok)}):  (a) {L.fmt(row['a']['boot'])}  "
          f"(a2) {L.fmt(row['a2']['boot'])}  (b) {L.fmt(row['b']['boot'])}  (x) {L.fmt(row['x']['boot'])}  "
          f"rho_guide {X.fmt_rho(row['rho_guide'])}")
        P(f"           pds (a) {row['pds']['a']:.4f} / (b) {row['pds']['b']:.4f} / (x) {row['pds']['x']:.4f}")
        del cross_r, cross_u

    res["seconds"] = round(time.time() - t0, 1)
    res["peak_rss_gib"] = L.peak_gib()
    np.savez_compressed(X.WORK / f"nadig_{WHICH}_296.npz", **store)
    (X.WORK / f"nadig_{WHICH}_296.json").write_text(json.dumps(res, indent=1, default=float))
    P(f"\ndone {res['seconds']}s, peak {res['peak_rss_gib']:.2f} GiB")
    P.save()


if __name__ == "__main__":
    main()
