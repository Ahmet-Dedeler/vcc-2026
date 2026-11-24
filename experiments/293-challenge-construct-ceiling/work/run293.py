"""293 -- the matched construct design on the 2025 VCC splits (PREREGISTRATION.md sec 3).

One split per invocation:  PYTHONDONTWRITEBYTECODE=1 uv run python run293.py <train|validation|test>

Reads (read-only) the split's h5ad in CSR row blocks: ONE pass gives both the extracted cells (uint16,
capped per PREREGISTRATION sec 5) and exact full-depth group sums (each target, each two-construct
construct, the two NTC-parity control halves).  Writes work/run293_<split>.{json,npz,log}.
"""
from __future__ import annotations

import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib293 as X  # noqa: E402
from lib293 import L  # noqa: E402


def stream_pass(path: Path, want_mask: np.ndarray, groups: np.ndarray, n_groups: int, chunk: int = X.CHUNK):
    """One pass: dense uint16 rows where want_mask, plus exact float64 group sums where groups >= 0."""
    with h5py.File(path, "r") as h:
        n, g = (int(v) for v in h["X"].attrs["shape"])
        indptr = h["X/indptr"][:]
        out = np.zeros((int(want_mask.sum()), g), dtype=np.uint16)
        acc = np.zeros((n_groups, g), dtype=np.float64)
        pos = np.cumsum(want_mask) - 1
        mx = 0.0
        for r0 in range(0, n, chunk):
            r1 = min(r0 + chunk, n)
            wm, gm = want_mask[r0:r1], groups[r0:r1]
            if not (wm.any() or (gm >= 0).any()):
                continue
            s, e = int(indptr[r0]), int(indptr[r1])
            Xb = sp.csr_matrix((h["X/data"][s:e].astype(np.float64), h["X/indices"][s:e],
                                indptr[r0:r1 + 1] - indptr[r0]), shape=(r1 - r0, g))
            keep = gm >= 0
            if keep.any():
                oh = sp.csr_matrix((np.ones(int(keep.sum())), (gm[keep], np.flatnonzero(keep))),
                                   shape=(n_groups, r1 - r0))
                acc += np.asarray((oh @ Xb).todense())
            if wm.any():
                D = np.asarray(Xb[np.flatnonzero(wm)].todense())
                mx = max(mx, float(D.max()) if D.size else 0.0)
                assert mx < 65535.0 and np.array_equal(D, np.rint(D)), "count overflow / non-integer"
                out[pos[r0:r1][wm]] = D.astype(np.uint16)
    cnt = np.bincount(groups[groups >= 0], minlength=n_groups)
    return out, acc, cnt, mx


def main() -> None:
    split = sys.argv[1]
    # optional 2nd arg: binomial thinning fraction applied to EVERY summed side (and the controls), so the
    # native 2025 depth (54,256 median UMI/cell) can be read at the 2026-faithful depth of h1-test (19,965).
    thin = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    soff = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    tag = ("" if thin >= 1.0 else f"_thin{thin:.3f}") + (f"_s{soff}" if soff else "")
    path = X.SPLITS[split]
    t0 = time.time()
    P = L.Logger(X.WORK / f"run293_{split}{tag}.log")
    P(f"293 run293.py {split} thin={thin:.4f} -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")

    with h5py.File(path, "r") as h:
        gid = X.labels_of(h, "obs/guide_id")
        tgt = X.labels_of(h, "obs/target_gene")
        bat = X.labels_of(h, "obs/batch")
        genes = X.labels_of(h, "var/_index")
        n_cells, n_genes = (int(v) for v in h["X"].attrs["shape"])
    cons = np.asarray([X.construct_of(g) for g in gid])
    is_ctrl = tgt == X.CTRL
    targets = np.asarray(sorted(set(tgt.tolist()) - {X.CTRL}))
    tpos = {t: i for i, t in enumerate(targets)}
    gpos = {g: i for i, g in enumerate(genes)}
    tcols = np.asarray([gpos[t] for t in targets], dtype=np.int64)
    scored = np.asarray(sorted(set(range(n_genes)) - set(tcols.tolist())), dtype=np.int64)
    spos = {c: i for i, c in enumerate(scored)}           # scored-column position (target cols absent)

    by_t = defaultdict(set)
    for t, c in zip(tgt, cons):
        by_t[t].add(c)
    multi = {t: sorted(by_t[t]) for t in targets if len(by_t[t]) > 1}
    mts = sorted(multi)
    lanes = sorted(set(bat.tolist()))
    lpar = np.asarray([lanes.index(b) % 2 for b in bat])
    ntc_cons = sorted(set(cons[is_ctrl].tolist()))
    ntc_par = {c: i % 2 for i, c in enumerate(ntc_cons)}
    ctrl_side = np.asarray([ntc_par[c] if is_ctrl[i] else -1 for i, c in enumerate(cons)])
    P(f"  {n_cells:,} cells x {n_genes} genes; {targets.size} targets, {len(set(cons.tolist()))} constructs, "
      f"{len(lanes)} lanes, {int(is_ctrl.sum()):,} NTC cells in {len(ntc_cons)} constructs "
      f"({int((ctrl_side == 0).sum()):,} / {int((ctrl_side == 1).sum()):,} by parity); "
      f"{len(mts)} two-construct targets {mts}; scored columns {scored.size}")

    # ---------------------------------------------------------------- what to extract, and the group map
    rng0 = np.random.default_rng(L.SEED + 7 * len(split))
    want = np.zeros(n_cells, dtype=bool)
    pool_src: dict[str, np.ndarray] = {}
    for t in targets:
        ix = np.flatnonzero(tgt == t)
        if ix.size > X.POOL_CAP:
            ix = np.sort(rng0.choice(ix, X.POOL_CAP, replace=False))
        pool_src[t] = ix
        want[ix] = True
    cons_src: dict[str, np.ndarray] = {}
    for t in mts:
        for c in multi[t]:
            ix = np.flatnonzero(cons == c)
            if ix.size > X.CONS_CAP:
                ix = np.sort(rng0.choice(ix, X.CONS_CAP, replace=False))
            cons_src[c] = ix
            want[ix] = True
    # group map: 0..T-1 each target FULL, T..T+C-1 each two-construct construct FULL, then the 2 control halves
    allc = [c for t in mts for c in multi[t]]
    T, C = targets.size, len(allc)
    G_CA, G_CB = T + C, T + C + 1
    groups = np.full(n_cells, -1, dtype=np.int64)
    for i, t in enumerate(targets):
        groups[tgt == t] = i
    cpos = {c: T + i for i, c in enumerate(allc)}
    for c in allc:                                          # constructs override their target's group
        groups[cons == c] = cpos[c]
    groups[ctrl_side == 0] = G_CA
    groups[ctrl_side == 1] = G_CB
    P(f"  extracting {int(want.sum()):,} cells (pool cap {X.POOL_CAP}/target, construct cap {X.CONS_CAP}); "
      f"{T + C + 2} full-depth groups")

    Xd, acc, gcnt, mx = stream_pass(path, want, groups, T + C + 2)
    dst = np.cumsum(want) - 1
    pool = {t: dst[pool_src[t]] for t in targets}
    cell = {c: dst[cons_src[c]] for c in allc}
    par = {c: lpar[cons_src[c]] for c in allc}
    P(f"  extracted {Xd.shape} uint16 (max count {mx:.0f}), {time.time() - t0:.0f}s, peak {L.peak_gib():.2f} GiB")

    # ---------------------------------------------------------------- full-depth reference quantities
    cA_full, cB_full = acc[G_CA], acc[G_CB]
    ctrl_pooled = cA_full + cB_full
    pc = L.prof(ctrl_pooled)
    # on-target knockdown per construct (its own gene column; excluded from every scored statistic)
    kd_full = {}
    for t in mts:
        for c in multi[t]:
            j = gpos[t]
            kd_full[c] = dict(cells=int(gcnt[cpos[c]]),
                              lfc=float(L.prof(acc[cpos[c]])[j] - pc[j]),
                              ctrl_prof=float(pc[j]))
    # on-target knockdown of EVERY target in the split (its pooled cells; the mechanism's base rate)
    kd_all = {t: float(L.prof(acc[tpos[t]] + sum(acc[cpos[c]] for c in multi.get(t, [])))[gpos[t]] - pc[gpos[t]])
              for t in targets}
    # full-depth response magnitude per target (residual row norm on scored columns, panel-centred)
    Dfull = L.prof(np.stack([acc[tpos[t]] + sum(acc[cpos[c]] for c in multi.get(t, [])) if t in multi
                             else acc[tpos[t]] for t in targets]))
    Dfull = Dfull[:, scored] - pc[scored][None, :]
    Rfull = Dfull - Dfull.mean(0, keepdims=True)
    mag_full = {t: float(np.linalg.norm(Rfull[tpos[t]])) for t in targets}
    P(f"  full-depth: controls {int(gcnt[G_CA]):,}/{int(gcnt[G_CB]):,}; knockdown lfc per construct "
      + ", ".join(f"{c} {v['lfc']:+.3f}" for c, v in kd_full.items()))

    # ---------------------------------------------------------------- the matched ladder
    if thin < 1.0:                      # thin the control sums once, with their own stream
        tr = np.random.default_rng(L.SEED + 991)
        cA_full, cB_full = L.thin(cA_full, thin, tr), L.thin(cB_full, thin, tr)

    def rowsum(idx, r=None):
        s = Xd[np.sort(idx)].sum(0, dtype=np.int64).astype(np.float64)
        return s if thin >= 1.0 else L.thin(s, thin, r)

    res: dict = {"split": split, "n_cells": n_cells, "n_targets": int(T), "n_multi": len(mts),
                 "thin": thin, "multi": multi, "scored_columns": int(scored.size), "lanes": len(lanes),
                 "controls": dict(A=int(gcnt[G_CA]), B=int(gcnt[G_CB])), "kd_full": kd_full,
                 "mag_full": mag_full, "kd_all": kd_all, "pool_cap": X.POOL_CAP, "cons_cap": X.CONS_CAP,
                 "cells_per_construct": {c: int(gcnt[cpos[c]]) for c in allc},
                 "cells_extracted_per_construct": {c: int(cell[c].size) for c in allc},
                 "rungs": {}}
    store: dict[str, np.ndarray] = {}

    for n in X.RUNGS:
        pool_ok = [t for t in targets if pool[t].size >= 2 * n]
        cohort = [t for t in mts
                  if max(cell[multi[t][0]].size, cell[multi[t][1]].size) >= 2 * n
                  and min(cell[multi[t][0]].size, cell[multi[t][1]].size) >= n]
        lane_ok = []
        for t in mts:
            c = max(multi[t], key=lambda c: cell[c].size)
            if (par[c] == 0).sum() >= n and (par[c] == 1).sum() >= n:
                lane_ok.append(t)
        if not cohort:
            P(f"  n={n}: no target qualifies"); continue
        K = len(cohort)
        arr = {k: np.full((X.NDRAW, K), np.nan) for k in
               ("a_cos", "b_cos", "d_cos", "bfix_cos", "a2_cos", "x_cos", "a_pds", "b_pds", "d_pds",
                "bfix_pds", "a2_pds", "x_pds", "a_own", "b_own", "d_own")}
        pan = np.full((X.NDRAW, len(pool_ok)), np.nan)      # the split's own panel-wide cell split
        panp = np.full((X.NDRAW, len(pool_ok)), np.nan)
        kd_n = {c: [] for t in cohort for c in multi[t]}
        for d in range(X.NDRAW):
            r = np.random.default_rng(L.SEED + 131 * d + 7 * n + len(split) + 100003 * soff)
            # (i) the matched-n centring / distractor pool, per side, rebuilt for this (n, draw)
            SA = {}; SB = {}
            for t in pool_ok:
                p = r.permutation(pool[t])
                SA[t] = rowsum(p[:n], r); SB[t] = rowsum(p[n:2 * n], r)
            base = np.stack([SA[t] for t in pool_ok]); baseB = np.stack([SB[t] for t in pool_ok])
            DA = (L.prof(base) - L.prof(cA_full)[None, :])[:, scored]
            DB = (L.prof(baseB) - L.prof(cB_full)[None, :])[:, scored]
            mA = DA.mean(0); mB = DB.mean(0)
            bi = {t: i for i, t in enumerate(pool_ok)}
            vp = L.views(base, baseB, cA_full, cB_full, scored, mA, mB)   # panel-wide baseline at this n
            pan[d] = vp["cos_resid"]; panp[d] = vp["pds"]
            # (ii) the three contrasts, sharing side A (one permutation of C_big per draw)
            sides = {k: {} for k in ("a", "b", "d", "bfix", "a2", "x")}
            for t in cohort:
                c1, c2 = multi[t]
                big, small = (c1, c2) if cell[c1].size >= cell[c2].size else (c2, c1)
                pb = r.permutation(cell[big])
                A = rowsum(pb[:n], r)
                sides["a"][t] = (A, rowsum(pb[n:2 * n], r))
                sides["b"][t] = (A, rowsum(r.choice(cell[small], n, replace=False), r))
                if cell[small].size >= 2 * n:          # the SMALLER construct's own cell halves
                    ps = r.permutation(cell[small])
                    sides["a2"][t] = (rowsum(ps[:n], r), rowsum(ps[n:2 * n], r))
                if t in lane_ok:
                    e = cell[big][par[big] == 0]; o = cell[big][par[big] == 1]
                    sides["d"][t] = (rowsum(r.choice(e, n, replace=False), r), rowsum(r.choice(o, n, replace=False), r))
                # sensitivity: construct 1 -> side A, construct 2 -> side B (fixed identity)
                if cell[c1].size >= n and cell[c2].size >= n:
                    sides["bfix"][t] = (rowsum(r.choice(cell[c1], n, replace=False), r),
                                        rowsum(r.choice(cell[c2], n, replace=False), r))
                for c in (c1, c2):
                    if cell[c].size >= n:
                        s = rowsum(r.choice(cell[c], n, replace=False), r)
                        kd_n[c].append(float(L.prof(s)[gpos[t]] - pc[gpos[t]]))
            # (iii) each contrast scored with the split's other targets as distractors
            # cross-target NULL: side A = target t's C_big, side B = the NEXT target's C_small
            for i_t, t in enumerate(cohort):
                if t not in sides["a"] or len(cohort) < 2:
                    continue
                u = cohort[(i_t + 1) % len(cohort)]
                c1u, c2u = multi[u]
                smallu = c2u if cell[c1u].size >= cell[c2u].size else c1u
                if cell[smallu].size >= n:
                    sides["x"][t] = (sides["a"][t][0], rowsum(r.choice(cell[smallu], n, replace=False), r))
            for key in ("a", "b", "d", "bfix", "a2", "x"):
                have = [t for t in cohort if t in sides[key]]
                if not have:
                    continue
                MA, MB = base.copy(), baseB.copy()
                rows = []
                for t in have:
                    if t not in bi:            # a two-construct target too small for the pool: append it
                        MA = np.vstack([MA, sides[key][t][0][None, :]])
                        MB = np.vstack([MB, sides[key][t][1][None, :]])
                        rows.append(MA.shape[0] - 1)
                    else:
                        MA[bi[t]] = sides[key][t][0]; MB[bi[t]] = sides[key][t][1]
                        rows.append(bi[t])
                v = L.views(MA, MB, cA_full, cB_full, scored, mA, mB)
                for j, t in enumerate(have):
                    k = cohort.index(t)
                    arr[f"{key}_cos"][d, k] = v["cos_resid"][rows[j]]
                    arr[f"{key}_pds"][d, k] = v["pds"][rows[j]]
                    if key not in ("bfix", "a2", "x"):
                        arr[f"{key}_own"][d, k] = v["cos_own"][rows[j]]
        # (iv) average the 8 draws, then the statistics
        m = {k: np.nanmean(arr[k], 0) for k in arr}
        row: dict = dict(n=n, cohort=cohort, n_cohort=K, n_pool=len(pool_ok), lane_cohort=lane_ok,
                         per_target={t: {k: (None if not np.isfinite(m[k][i]) else float(m[k][i]))
                                         for k in m} for i, t in enumerate(cohort)},
                         kd_n={c: (float(np.mean(v)) if v else None) for c, v in kd_n.items()})
        for k in m:
            row[k] = dict(**L.stats(m[k]), boot=L.boot_median(m[k]))
        row["a_minus_b"] = L.boot_median_diff(m["a_cos"], m["b_cos"])
        row["a_minus_d"] = L.boot_median_diff(m["a_cos"], m["d_cos"])
        # rho_guide = median(b)/median(a): ratio of bootstrapped medians (285 k562_guide.py) + per-target
        good = np.isfinite(m["a_cos"]) & np.isfinite(m["b_cos"])
        rb_rng = np.random.default_rng(L.SEED)
        ii = rb_rng.integers(0, int(good.sum()), size=(L.NBOOT, int(good.sum())))
        rr = np.median(m["b_cos"][good][ii], 1) / np.median(m["a_cos"][good][ii], 1)
        per_ratio = np.where(np.abs(m["a_cos"]) > 1e-9, m["b_cos"] / m["a_cos"], np.nan)
        row["rho_guide"] = dict(point=float(np.median(m["b_cos"][good]) / np.median(m["a_cos"][good])),
                                lo=float(np.percentile(rr, 2.5)), hi=float(np.percentile(rr, 97.5)),
                                n=int(good.sum()), median_of_ratios=float(np.nanmedian(per_ratio)),
                                per_target={t: (None if not np.isfinite(per_ratio[i]) else float(per_ratio[i]))
                                            for i, t in enumerate(cohort)})
        xp, xpp = np.nanmean(pan, 0), np.nanmean(panp, 0)
        row["panel"] = dict(n_targets=len(pool_ok), cos=dict(**L.stats(xp), boot=L.boot_median(xp)),
                            pds=dict(**L.stats(xpp), boot=L.boot_median(xpp)),
                            multi_rank={t: (float((xp < m["a_cos"][i]).mean()) if np.isfinite(m["a_cos"][i]) else None)
                                        for i, t in enumerate(cohort)},
                            per_target={t: float(xp[j]) for j, t in enumerate(pool_ok)})
        res["rungs"][n] = row
        for k in m:
            store[f"{n}|{k}"] = m[k]
        store[f"{n}|targets"] = np.asarray(cohort)
        store[f"{n}|panel_cos"] = xp; store[f"{n}|panel_pds"] = xpp
        store[f"{n}|panel_targets"] = np.asarray(pool_ok)
        P(f"  n={n:3d}/side  cohort {K} (pool {len(pool_ok)}/{T}, lanes {len(lane_ok)}):  "
          f"(a) {L.fmt(row['a_cos']['boot'])}  (b) {L.fmt(row['b_cos']['boot'])}  (d) {L.fmt(row['d_cos']['boot'])}  "
          f"rho_guide {row['rho_guide']['point']:+.3f} [{row['rho_guide']['lo']:+.3f}, {row['rho_guide']['hi']:+.3f}]")
        P(f"           pds (a) {row['a_pds']['median']:+.4f} / (b) {row['b_pds']['median']:+.4f} / "
          f"(d) {row['d_pds']['median']:+.4f}   (b) fixed-identity check {L.fmt(row['bfix_cos']['boot'])}")
        P(f"           panel baseline (all {len(pool_ok)} pool targets, cell split): cos "
          f"{L.fmt(row['panel']['cos']['boot'])} (<=0 {row['panel']['cos']['share_le0']:.3f}), pds median "
          f"{row['panel']['pds']['median']:+.4f}; un-centred (a) {np.nanmedian(m['a_own']):+.4f} / "
          f"(b) {np.nanmedian(m['b_own']):+.4f}")

    res["seconds"] = round(time.time() - t0, 1)
    res["peak_rss_gib"] = L.peak_gib()
    np.savez_compressed(X.WORK / f"run293_{split}{tag}.npz", **store)
    (X.WORK / f"run293_{split}{tag}.json").write_text(json.dumps(res, indent=1, default=float))
    P(f"\ndone {res['seconds']}s, peak {res['peak_rss_gib']:.2f} GiB")
    P.save()


if __name__ == "__main__":
    main()
