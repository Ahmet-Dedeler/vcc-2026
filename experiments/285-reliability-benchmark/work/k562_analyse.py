"""285 -- the K562 rows of the benchmark: (a) the cell ladder, (a') UMI thinning, (b) the guide construct,
(d) the gemgroup, (e) the SAME LINE IN A DIFFERENT STUDY (GWPS vs the essential screen), all at matched cells
per side, plus 267's cached reliability medians as KA-1.

CENTRING (fixed after a first run exposed the artefact).  Every row's residual is
    r = log1p(TS*S/sum S)  -  the mean of that statistic over a MATCHED-n, MATCHED-DEPTH random sample of
                              every target of the source, taken on the same side
so the centring vector is rebuilt for every (n, draw) and for every thinning fraction.  Centring a thinned
side on a FULL-DEPTH mean leaves a large shared depth artefact in both sides and inflates the cosine (in the
first run it made 1/8-depth read ABOVE full depth).  The convention is control-free (the control cancels in
the residual) and cohort-stable, so the 8-target guide rows sit in the same residual space as the 258-target
ladder.

    PYTHONDONTWRITEBYTECODE=1 uv run python k562_analyse.py
"""
from __future__ import annotations

import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib285 as L  # noqa: E402

P = L.Logger(L.WORK / "k562_analyse.log")
FAIL: list[str] = []
E267 = L.ROOT / "experiments/267-replicate-validated-oracle/work/half"
RUN = (25, 50, 100)
FRACS = (1.0, 0.5, 0.25, 0.125)


def main() -> None:
    t0 = time.time()
    P(f"285 k562_analyse.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    out: dict = {}
    ka: list[dict] = []

    # ---------------------------------------------------------------- KA-1: 267's cached reliability medians
    want = {"k562_gwps": 0.0471, "orion_hct116": 0.0276, "orion_hek293t": 0.0282}
    got = {}
    for line in want:
        z = np.load(E267 / f"truth-{line}.npz", allow_pickle=True)
        got[line] = float(np.round(np.nanmedian(np.asarray(z["rel|batch|cos_resid"], dtype=np.float64)), 4))
    ok = all(abs(got[k] - want[k]) < 1e-12 for k in want)
    ka.append(dict(row="KA-1 267's batch-split residual-cosine medians from its cached truth-<fold>.npz "
                       "(`rel|batch|cos_resid`)", required="0.0e+00", want=want, got=got, ok=bool(ok)))
    P(f"KA-1: 267 batch medians {got} vs {want} -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        FAIL.append("KA-1")

    # ---------------------------------------------------------------- the panel cells
    m = np.load(L.WORK / "k562_gwps_meta.npz", allow_pickle=True)
    genes = m["genes"].astype(str)
    pan_t = m["panel_targets"].astype(str)
    pert = m["pan_pert"].astype(str); batch = m["pan_batch"].astype(str); guide = m["pan_guide"].astype(str)
    ctrl = np.asarray(m["ctrl"], dtype=np.float64)
    gpos = {g: i for i, g in enumerate(genes)}
    tcols = np.asarray([gpos[t] for t in pan_t if t in gpos], dtype=np.int64)
    scored = np.asarray(sorted(set(range(len(genes))) - set(tcols.tolist())), dtype=np.int64)
    idx_t = {t: np.flatnonzero(pert == t) for t in pan_t}
    ncell = np.asarray([idx_t[t].size for t in pan_t])
    batches = sorted(set(batch.tolist())); bpar = np.asarray([batches.index(b) % 2 for b in batch])
    Xd = np.asarray(np.load(L.WORK / "k562_gwps_panelX.npy", mmap_mode="r"), dtype=np.float32)
    P(f"  panel: {len(pan_t)} targets, {Xd.shape[0]:,} cells, {len(genes)} genes ({scored.size} scored), "
      f"cells/target median {np.median(ncell):.0f}; controls {m['n_ctrl'].tolist()}; "
      f"median UMI/cell {float(np.median(np.asarray(m['pan_umi']))):.0f}")

    def sums(sel):
        return Xd[np.sort(sel)].sum(0, dtype=np.float64)

    full = np.stack([sums(idx_t[t]) for t in pan_t])
    mag = np.linalg.norm(L.resid(L.delta(full, ctrl[0])[:, scored]), axis=1)

    def cvec(pm, c):
        """centring vector in delta space; subtracting it leaves prof(S) - pm, which is control-free."""
        return pm - L.prof(c)[scored]

    arrays: dict[str, np.ndarray] = {}
    fixed = ncell >= 2 * max(RUN)
    P(f"\n  FIXED COHORT = the {int(fixed.sum())} targets with >= {2*max(RUN)} cells; all-available cohorts: "
      + ", ".join(f"n={n}: {int((ncell >= 2*n).sum())}" for n in RUN))

    pg = defaultdict(list)
    for t in pan_t:
        for g in sorted(set(guide[idx_t[t]].tolist())):
            pg[t].append(g)
    ginfo = {}
    for t, gs in pg.items():
        if len(gs) > 1:
            cs = sorted([idx_t[t][guide[idx_t[t]] == g] for g in gs], key=lambda a: -a.size)
            ginfo[t] = (cs[0], cs[1], int(min(cs[0].size // 2, cs[1].size)))
    P("  two-construct panel targets: " + "; ".join(f"{t} {ginfo[t][0].size}/{ginfo[t][1].size} (n_t {ginfo[t][2]})"
                                                    for t in sorted(ginfo)))

    lad = {n: np.full((L.NDRAW, len(pan_t)), np.nan) for n in RUN}
    ladp = {n: np.full((L.NDRAW, len(pan_t)), np.nan) for n in RUN}
    thn = {(n, f): np.full((L.NDRAW, len(pan_t)), np.nan) for n in RUN for f in FRACS}
    gem = {n: np.full((L.NDRAW, len(pan_t)), np.nan) for n in RUN}
    GN = (35, 75)
    gtar = {gn: sorted([t for t in ginfo if ginfo[t][2] >= gn]) for gn in GN}
    gud = {gn: np.full((L.NDRAW, len(gtar[gn])), np.nan) for gn in GN}
    gcl = {gn: np.full((L.NDRAW, len(gtar[gn])), np.nan) for gn in GN}

    for n in RUN:
        for d in range(L.NDRAW):
            r = np.random.default_rng(L.SEED + 11 * d + n)
            keep, SA, SB = [], [], []
            for j, t in enumerate(pan_t):
                ix = idx_t[t]
                if ix.size < 2 * n:
                    continue
                p = r.permutation(ix)
                keep.append(j); SA.append(sums(p[:n])); SB.append(sums(p[n:2 * n]))
            keep = np.asarray(keep); SA, SB = np.stack(SA), np.stack(SB)
            pmA, pmB = L.prof(SA)[:, scored].mean(0), L.prof(SB)[:, scored].mean(0)
            v = L.views(SA, SB, ctrl[0], ctrl[0], scored, mA=cvec(pmA, ctrl[0]), mB=cvec(pmB, ctrl[0]))
            lad[n][d, keep] = v["cos_resid"]; ladp[n][d, keep] = v["pds"]
            thn[(n, 1.0)][d, keep] = v["cos_resid"]
            for f in FRACS:
                if f == 1.0:
                    continue
                rt = np.random.default_rng(L.SEED + 313 * d + int(1000 * f) + n)
                TA, TB = L.thin(SA, f, rt), L.thin(SB, f, rt)
                qA, qB = L.prof(TA)[:, scored].mean(0), L.prof(TB)[:, scored].mean(0)
                cA_t, cB_t = L.thin(ctrl[0], f, rt), L.thin(ctrl[0], f, rt)
                thn[(n, f)][d, keep] = L.views(TA, TB, cA_t, cB_t, scored,
                                               mA=cvec(qA, cA_t), mB=cvec(qB, cB_t))["cos_resid"]
            keep2, GA, GB = [], [], []
            for j, t in enumerate(pan_t):
                ix = idx_t[t]
                a, b = ix[bpar[ix] == 0], ix[bpar[ix] == 1]
                if a.size < n or b.size < n:
                    continue
                keep2.append(j)
                GA.append(sums(r.choice(a, n, replace=False))); GB.append(sums(r.choice(b, n, replace=False)))
            if len(keep2) >= 5:
                gem[n][d, np.asarray(keep2)] = L.views(np.stack(GA), np.stack(GB), ctrl[1], ctrl[2], scored,
                                                       mA=cvec(pmA, ctrl[1]), mB=cvec(pmB, ctrl[2]))["cos_resid"]
        P(f"    n={n} done, {time.time()-t0:.0f}s, peak {L.peak_gib():.2f} GiB")

    # ---------------------------------------------------------------- the guide rows (their own n)
    for gn in GN:
        ts = gtar[gn]
        if not ts:
            continue
        for d in range(L.NDRAW):
            r = np.random.default_rng(L.SEED + 97 * d + gn)
            allA, allB = [], []
            for t in pan_t:                                  # matched-depth centring sample
                ix = idx_t[t]
                if ix.size < 2 * gn:
                    continue
                p = r.permutation(ix)
                allA.append(sums(p[:gn])); allB.append(sums(p[gn:2 * gn]))
            pmA = L.prof(np.stack(allA))[:, scored].mean(0); pmB = L.prof(np.stack(allB))[:, scored].mean(0)
            SG1, SG2, SA1, SB1 = [], [], [], []
            for t in ts:
                big, small, _ = ginfo[t]
                pb = r.permutation(big)
                SA1.append(sums(pb[:gn])); SB1.append(sums(pb[gn:2 * gn]))
                SG1.append(sums(pb[:gn])); SG2.append(sums(r.choice(small, gn, replace=False)))
            gud[gn][d] = L.views(np.stack(SG1), np.stack(SG2), ctrl[0], ctrl[0], scored,
                                 mA=cvec(pmA, ctrl[0]), mB=cvec(pmB, ctrl[0]))["cos_resid"]
            gcl[gn][d] = L.views(np.stack(SA1), np.stack(SB1), ctrl[0], ctrl[0], scored,
                                 mA=cvec(pmA, ctrl[0]), mB=cvec(pmB, ctrl[0]))["cos_resid"]

    del Xd

    # ---------------------------------------------------------------- report
    out["ladder_a"] = {}
    P("\n(a) K562-GWPS panel targets, same-guide random cell halves, 8 draws:")
    for n in RUN:
        x, xp = np.nanmean(lad[n], 0), np.nanmean(ladp[n], 0)
        arrays[f"a|{n}"] = x; arrays[f"a_pds|{n}"] = xp
        out["ladder_a"][n] = dict(all=dict(**L.stats(x), boot=L.boot_median(x), pds_mean=float(np.nanmean(xp))),
                                  fixed=dict(**L.stats(x[fixed]), boot=L.boot_median(x[fixed]),
                                             pds_mean=float(np.nanmean(xp[fixed]))))
        P(f"  n={n:4d}/side  ALL n_t {out['ladder_a'][n]['all']['n']:3d} median "
          f"{L.fmt(out['ladder_a'][n]['all']['boot'])} <=0 {out['ladder_a'][n]['all']['share_le0']:.3f} | "
          f"FIXED n_t {out['ladder_a'][n]['fixed']['n']:3d} median {L.fmt(out['ladder_a'][n]['fixed']['boot'])} "
          f"<=0 {out['ladder_a'][n]['fixed']['share_le0']:.3f} | pds(all) {out['ladder_a'][n]['all']['pds_mean']:.4f}")

    out["thin"] = {}
    P("\n(a') thinned counts at fixed cells per side (fixed cohort, centring rebuilt at the thinned depth):")
    for (n, f), a in sorted(thn.items()):
        x = np.nanmean(a, 0)
        arrays[f"thin|{n}|{f}"] = x
        out["thin"][f"n{n}_f{f}"] = dict(n_cells=n, frac=f, **L.stats(x[fixed]), boot=L.boot_median(x[fixed]),
                                         all=L.boot_median(x))
        P(f"  n={n:3d} f={f:<6g}  fixed-cohort median {L.fmt(out['thin'][f'n{n}_f{f}']['boot'])}  "
          f"all-cohort {L.fmt(out['thin'][f'n{n}_f{f}']['all'])}")

    out["gemgroup"] = {}
    P("\n(d) gemgroup-parity halves at matched cells per side:")
    for n in RUN:
        x = np.nanmean(gem[n], 0)
        arrays[f"gem|{n}"] = x
        out["gemgroup"][n] = dict(**L.stats(x), boot=L.boot_median(x),
                                  cells_minus_gemgroup=L.boot_median_diff(arrays[f"a|{n}"], x))
        if out["gemgroup"][n]["n"]:
            P(f"  n={n:4d}/side  median {L.fmt(out['gemgroup'][n]['boot'])}  n_t {out['gemgroup'][n]['n']}  "
              f"<=0 {out['gemgroup'][n]['share_le0']:.3f}  |  CELLS − GEMGROUP "
              f"{L.fmt(out['gemgroup'][n]['cells_minus_gemgroup'])}")

    out["guide"] = {}
    P("\n(b) guide construct vs cells of ONE construct, at the SAME n:")
    for gn in GN:
        ts = gtar[gn]
        if not ts:
            continue
        xg, xa = np.nanmean(gud[gn], 0), np.nanmean(gcl[gn], 0)
        arrays[f"guide|{gn}"] = xg; arrays[f"guide_cellctrl|{gn}"] = xa
        rat = xg / np.where(np.abs(xa) > 1e-9, xa, np.nan)
        out["guide"][gn] = dict(targets=ts, guide=L.stats(xg), cells=L.stats(xa),
                                diff=L.boot_median_diff(xa, xg), ratio_median=float(np.nanmedian(rat)),
                                per_target={t: dict(guide=float(xg[i]), cells=float(xa[i])) for i, t in enumerate(ts)})
        P(f"  n={gn}/side, {len(ts)} targets: GUIDE {L.fmt(L.boot_median(xg))}  CELLS {L.fmt(L.boot_median(xa))}  "
          f"CELLS − GUIDE {L.fmt(out['guide'][gn]['diff'])}  (guide/cells ratio median "
          f"{out['guide'][gn]['ratio_median']:+.3f})")
        for i, t in enumerate(ts):
            P(f"      {t:8s} guide {xg[i]:+.4f}  cells {xa[i]:+.4f}")

    # ---------------------------------------------------------------- (e) the same line in another study
    P("\n(e) DIFFERENT STUDY, SAME LINE: K562-GWPS vs the K562 essential screen (Replogle 2022, two "
      "independent experiments), on their shared perturbations:")
    ma = np.load(L.WORK / "k562_gwps_meta.npz", allow_pickle=True)
    mb = np.load(L.WORK / "k562_essential_meta.npz", allow_pickle=True)
    ga_, gb_ = ma["genes"].astype(str), mb["genes"].astype(str)
    shared = ma["shared"].astype(str)
    assert list(shared) == list(mb["shared"].astype(str))
    common_g = sorted(set(ga_.tolist()) & set(gb_.tolist()))
    pa = {g: i for i, g in enumerate(ga_)}; pb = {g: i for i, g in enumerate(gb_)}
    ia = np.asarray([pa[g] for g in common_g], dtype=np.int64)
    ib = np.asarray([pb[g] for g in common_g], dtype=np.int64)
    tset = set(shared.tolist())
    keep_col = np.flatnonzero([g not in tset for g in common_g])
    P(f"  {len(shared)} shared perturbations; common gene axis {len(common_g)} ({keep_col.size} scored after "
      f"removing the shared targets' own genes); UMI/cell median GWPS {float(ma['umi_median']):.0f} / "
      f"essential {float(mb['umi_median']):.0f}")
    BA = np.load(L.WORK / "k562_gwps_study_blocks.npy", mmap_mode="r")
    BB = np.load(L.WORK / "k562_essential_study_blocks.npy", mmap_mode="r")
    nA, nB = np.asarray(ma["ncell_shared"]), np.asarray(mb["ncell_shared"])
    cA = np.asarray(ma["ctrl"], dtype=np.float64)[0][ia]
    cB = np.asarray(mb["ctrl"], dtype=np.float64)[0][ib]
    EDG = (0, 25, 50, 100)
    nd = BA.shape[0]
    out["study"] = {}
    for n in (25, 50):
        cohort = np.flatnonzero((nA >= 2 * n) & (nB >= 2 * n))
        acc = {k: np.full((nd, cohort.size), np.nan) for k in ("cross1", "cross2", "withinA", "withinB")}
        for d in range(nd):
            def side(Bk, ii, lo, hi):
                S = np.zeros((cohort.size, len(common_g)))
                for b in range(len(EDG) - 1):
                    if EDG[b] >= lo and EDG[b + 1] <= hi:
                        S += np.asarray(Bk[d, cohort * (len(EDG) - 1) + b][:, ii], dtype=np.float64)
                return S
            A1, A2 = side(BA, ia, 0, n), side(BA, ia, n, 2 * n)
            B1, B2 = side(BB, ib, 0, n), side(BB, ib, n, 2 * n)
            acc["withinA"][d] = L.views(A1, A2, cA, cA, keep_col)["cos_resid"]
            acc["withinB"][d] = L.views(B1, B2, cB, cB, keep_col)["cos_resid"]
            acc["cross1"][d] = L.views(A1, B1, cA, cB, keep_col)["cos_resid"]
            acc["cross2"][d] = L.views(A2, B2, cA, cB, keep_col)["cos_resid"]
        r_ = {k: np.nanmean(v, 0) for k, v in acc.items()}
        cross = np.nanmean(np.vstack([r_["cross1"], r_["cross2"]]), 0)
        within = np.nanmean(np.vstack([r_["withinA"], r_["withinB"]]), 0)
        arrays[f"study_cross|{n}"] = cross; arrays[f"study_within|{n}"] = within
        arrays[f"study_cohort|{n}"] = cohort
        out["study"][n] = dict(n_targets=int(cohort.size), cross=dict(**L.stats(cross), boot=L.boot_median(cross)),
                               within_gwps=dict(**L.stats(r_["withinA"]), boot=L.boot_median(r_["withinA"])),
                               within_ess=dict(**L.stats(r_["withinB"]), boot=L.boot_median(r_["withinB"])),
                               within_mean=dict(**L.stats(within), boot=L.boot_median(within)),
                               within_minus_cross=L.boot_median_diff(within, cross),
                               ratio=float(np.nanmedian(cross) / np.nanmedian(within)))
        P(f"  n={n:3d}/side, {cohort.size} shared targets: WITHIN-study cell halves "
          f"{L.fmt(out['study'][n]['within_mean']['boot'])} (GWPS {np.nanmedian(r_['withinA']):+.4f} / essential "
          f"{np.nanmedian(r_['withinB']):+.4f})  vs  ACROSS-STUDY {L.fmt(out['study'][n]['cross']['boot'])}   "
          f"WITHIN − ACROSS {L.fmt(out['study'][n]['within_minus_cross'])}  (across/within "
          f"{out['study'][n]['ratio']:.3f})")

    out["ka"] = ka
    out["failures"] = FAIL
    out["panel_targets"] = pan_t.tolist()
    out["ncell"] = ncell.tolist()
    out["guide_targets_n"] = {t: ginfo[t][2] for t in ginfo}
    out["seconds"] = round(time.time() - t0, 1)
    out["peak_rss_gib"] = L.peak_gib()
    np.savez_compressed(L.WORK / "k562_285.npz", panel_targets=pan_t, ncell=ncell, magnitude=mag,
                        fixed=fixed, shared=shared, **arrays)
    L.Logger.dump(P, out, L.WORK / "k562_285.json")
    P(f"\ndone: {len(FAIL)} KA failures, {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB")
    P.save()
    if FAIL:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
