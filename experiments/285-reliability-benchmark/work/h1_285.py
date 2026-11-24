"""285 -- h1-test (H1 hESC, 10x Flex, 85 targets x exactly 400 cells, 48 lanes, 5 two-construct targets).

The cheapest source that carries (a) the cell ladder to 200/side, (a') UMI thinning, (b) the guide construct
and (d) the lane, all at MATCHED cells per side, on ONE target set.

    PYTHONDONTWRITEBYTECODE=1 uv run python h1_285.py

Reads (read-only): data/harness/h1-test.h5ad, data/external/vcc2025/test/adata_Test.h5ad (obs only, for
guide_id -- 273's join), experiments/273-.../work/reliability273.json (KA-2 / KA-2b reference values).
Writes work/h1_285.{npz,json,log}.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib285 as L  # noqa: E402

H1 = L.ROOT / "data/harness/h1-test.h5ad"
SRC = L.ROOT / "data/external/vcc2025/test/adata_Test.h5ad"
E273 = L.ROOT / "experiments/273-flex-reference-instrument/work"
P = L.Logger(L.WORK / "h1_285.log")
FAIL: list[str] = []


def labels_of(h, path):
    o = h[path]
    if isinstance(o, h5py.Group) and "categories" in o:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in o["categories"][:]])
        return cats[o["codes"][:]]
    v = o["values"][:] if isinstance(o, h5py.Group) else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def main() -> None:
    t0 = time.time()
    P(f"285 h1_285.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")

    with h5py.File(H1, "r") as h:
        lab = labels_of(h, "obs/target_gene")
        bc = labels_of(h, "obs/_index")
        genes = labels_of(h, "var/_index")
        n_cells, n_genes = (int(x) for x in h["X"].attrs["shape"])
        X = sp.csr_matrix((h["X/data"][:].astype(np.float32), h["X/indices"][:], h["X/indptr"][:]),
                          shape=(n_cells, n_genes))
    P(f"  h1-test {X.shape}, nnz {X.nnz:,} ({X.nnz/n_cells:.0f}/cell), peak {L.peak_gib():.2f} GiB, {time.time()-t0:.0f}s")
    with h5py.File(SRC, "r") as h:
        src_idx = labels_of(h, "obs/_index")
        src_guide = labels_of(h, "obs/guide_id")
    pos = {b: i for i, b in enumerate(src_idx)}
    hit = np.asarray([pos[b] for b in bc], dtype=np.int64)
    guide = src_guide[hit]
    lane = np.asarray([b.split("-", 1)[1] for b in bc])
    lanes = sorted(set(lane.tolist()))
    lane_par = np.asarray([lanes.index(x) % 2 for x in lane])

    targets = np.asarray(sorted(set(lab.tolist()) - {"non-targeting"}))
    tpos = {t: i for i, t in enumerate(targets)}
    is_ctrl = lab == "non-targeting"
    gpos = {g: i for i, g in enumerate(genes)}
    tcols = np.asarray([gpos[t] for t in targets], dtype=np.int64)
    scored = np.asarray(sorted(set(range(n_genes)) - set(tcols.tolist())), dtype=np.int64)
    ctrl_full = np.asarray(X[np.flatnonzero(is_ctrl)].sum(0), dtype=np.float64).ravel()
    P(f"  {targets.size} targets, {len(lanes)} lanes, {int(is_ctrl.sum())} controls, scored columns {scored.size}")

    # per-target cell index sets
    idx_t = {t: np.flatnonzero(lab == t) for t in targets}
    assert all(len(v) == 400 for v in idx_t.values())
    dense = {t: np.asarray(X[idx_t[t]].todense(), dtype=np.float64) for t in targets}   # 85 x (400 x 18,080)
    P(f"  per-target dense blocks built, peak {L.peak_gib():.2f} GiB, {time.time()-t0:.0f}s")

    full = np.stack([dense[t].sum(0) for t in targets])
    out: dict = {}
    ka: list[dict] = []

    def score(SA, SB, cA=ctrl_full, cB=ctrl_full):
        return L.views(SA, SB, cA, cB, scored)

    # ------------------------------------------------------------------ KA-2 / KA-2b: 273's lane parity 200/200
    A = np.stack([dense[t][lane_par[idx_t[t]] == 0].sum(0) for t in targets])
    B = np.stack([dense[t][lane_par[idx_t[t]] == 1].sum(0) for t in targets])
    nA = np.asarray([int((lane_par[idx_t[t]] == 0).sum()) for t in targets])
    cA = np.asarray(X[np.flatnonzero(is_ctrl & (lane_par == 0))].sum(0), dtype=np.float64).ravel()
    cB = np.asarray(X[np.flatnonzero(is_ctrl & (lane_par == 1))].sum(0), dtype=np.float64).ravel()
    v = L.views(A, B, cA, cB, scored)
    ref = json.load(open(E273 / "reliability273.json"))["reliability"]
    want = ref["batch|resid|disjoint"]["median"]
    got = float(np.nanmedian(v["cos_resid"]))
    ok = abs(got - want) <= 1e-4
    ka.append(dict(row="KA-2 273's h1-test lane-parity 200/200 residual median recomputed here",
                   required="|d| <= 1e-4 vs reliability273.json", want=want, got=got, diff=abs(got - want), ok=bool(ok)))
    P(f"KA-2: lane parity residual median {got:+.6f} vs 273's {want:+.6f} (|d| {abs(got-want):.2e}) "
      f"{'PASS' if ok else 'FAIL'}; cells/half {nA.min()}..{nA.max()} (median {np.median(nA):.0f})")
    if not ok:
        FAIL.append("KA-2")
    want2 = ref["rand24_m200|resid"]["mean"]
    ka.append(dict(row="KA-2b 159's resid_disjoint mean at 200 cells/side (cache read)", required="0.0e+00",
                   want=0.330039, got=want2, diff=abs(want2 - 0.330039), ok=bool(abs(want2 - 0.330039) < 1e-6)))
    P(f"KA-2b: 159's stored resid_disjoint mean at 200/side {want2:.6f} (159: 0.330039) "
      f"{'PASS' if abs(want2-0.330039) < 1e-6 else 'FAIL'}")
    if abs(want2 - 0.330039) >= 1e-6:
        FAIL.append("KA-2b")
    out["lane_parity_full"] = dict(cells_per_side=float(np.median(nA)), **L.stats(v["cos_resid"]),
                                   pds_mean=float(np.nanmean(v["pds"])), cos_agg=v["cos_agg"])

    # KA-4: the union of two halves reproduces the full row (integer identity)
    d4 = float(np.abs(A + B - full).max())
    ka.append(dict(row="KA-4 union of the two lane halves == the full per-target sums", required="0", got=d4, ok=bool(d4 == 0.0)))
    P(f"KA-4: max|A + B - full| = {d4:.1f} {'PASS' if d4 == 0 else 'FAIL'}")
    if d4 != 0.0:
        FAIL.append("KA-4")

    # KA-5: pds of a view against itself
    p_self, _ = L.per_target(L.delta(full, ctrl_full)[:, scored], L.delta(full, ctrl_full)[:, scored])
    ka.append(dict(row="KA-5 pds of a view against itself", required="1.0 on every target",
                   got=float(p_self.min()), ok=bool(np.allclose(p_self, 1.0))))
    P(f"KA-5: self-pds min {p_self.min():.6f} {'PASS' if np.allclose(p_self, 1.0) else 'FAIL'}")
    if not np.allclose(p_self, 1.0):
        FAIL.append("KA-5")

    # ------------------------------------------------------------------ response-magnitude strata (full depth)
    Rfull = L.resid(L.delta(full, ctrl_full)[:, scored])
    mag = np.linalg.norm(Rfull, axis=1)
    dec = np.searchsorted(np.percentile(mag, np.arange(10, 100, 10)), mag)
    out["magnitude"] = dict(median=float(np.median(mag)), deciles=[float(x) for x in np.percentile(mag, np.arange(10, 100, 10))])

    # ------------------------------------------------------------------ (a) the cell ladder, 8 draws
    rng = np.random.default_rng(L.SEED)
    RUN = [25, 50, 100, 200]
    lad = {n: [] for n in RUN}
    lad_pds = {n: [] for n in RUN}
    thin_rows: dict[tuple[int, float], list] = {}
    for d in range(L.NDRAW):
        perm = {t: rng.permutation(400) for t in targets}
        blocks = {}
        for n in RUN:
            SA = np.stack([dense[t][perm[t][:n]].sum(0) for t in targets])
            SB = np.stack([dense[t][perm[t][n:2 * n]].sum(0) for t in targets])
            blocks[n] = (SA, SB)
            v = score(SA, SB)
            lad[n].append(v["cos_resid"]); lad_pds[n].append(v["pds"])
        # (a') thinning at fixed n, independent binomial on each side
        for n in (100, 200):
            SA, SB = blocks[n]
            for f in (1.0, 0.5, 0.25, 0.125):
                r2 = np.random.default_rng(L.SEED + 1000 * d + int(1000 * f) + n)
                v = score(L.thin(SA, f, r2), L.thin(SB, f, r2))
                thin_rows.setdefault((n, f), []).append(v["cos_resid"])
        # the fixed-UMI-budget contrast: (2m, 1/2) vs (m, 1)
        for m in (25, 50, 100):
            SA, SB = blocks[2 * m]
            r2 = np.random.default_rng(L.SEED + 7000 + 13 * d + m)
            v = score(L.thin(SA, 0.5, r2), L.thin(SB, 0.5, r2))
            thin_rows.setdefault((2 * m, 0.5), []).append(v["cos_resid"])
    ladder = {n: np.nanmean(np.vstack(lad[n]), 0) for n in RUN}
    ladder_pds = {n: np.nanmean(np.vstack(lad_pds[n]), 0) for n in RUN}
    P("\n(a) same-guide random cell halves, 8 draws averaged, 85 targets, controls shared (the scorer's rule):")
    out["ladder_a"] = {}
    for n in RUN:
        b = L.boot_median(ladder[n])
        out["ladder_a"][n] = dict(**L.stats(ladder[n]), boot=b, pds_mean=float(np.nanmean(ladder_pds[n])),
                                  pds=L.stats(ladder_pds[n]))
        P(f"  n={n:4d}/side  cos_resid median {L.fmt(b)}  mean {out['ladder_a'][n]['mean']:+.4f}  "
          f"<=0 {out['ladder_a'][n]['share_le0']:.3f}  |  pds mean {out['ladder_a'][n]['pds_mean']:.4f}")
    out["thin"] = {}
    P("\n(a') same cells, thinned counts (binomial), 8 draws:")
    for (n, f), arr in sorted(thin_rows.items()):
        x = np.nanmean(np.vstack(arr), 0)
        out["thin"][f"n{n}_f{f}"] = dict(n_cells=n, frac=f, **L.stats(x), boot=L.boot_median(x))
        P(f"  n={n:4d}/side f={f:<6g} cos_resid median {x[np.isfinite(x)].mean()*0 + np.nanmedian(x):+.4f}  "
          f"<=0 {float(np.nanmean(x <= 0)):.3f}")
    # store the per-target arrays for the analysis stage
    arrays = {f"ladder_a|{n}": ladder[n] for n in RUN}
    arrays.update({f"ladder_a_pds|{n}": ladder_pds[n] for n in RUN})
    arrays.update({f"thin|{n}|{f}": np.nanmean(np.vstack(v), 0) for (n, f), v in thin_rows.items()})

    # the GLOBAL centring vector (273's `halves_sub` convention), used for every SUBSET row so that a 4-row
    # matrix is never centred on its own mean: m = the full-depth mean profile over all 85 targets, in delta
    # space against each side's own control.
    prof_full = L.prof(full)[:, scored].mean(0)

    def gm(cv):
        return prof_full - L.prof(cv)[scored]

    def cv_(pm, c):
        """centring vector in delta space; subtracting it leaves prof(S) - pm (control-free).  `pm` is the
        mean profile over ALL 85 targets of a MATCHED-n random sample on the same side, so the centring is
        depth-matched -- centring a n=25 side on the 400-cell panel mean leaves a shared depth artefact in
        both sides and inflates the cosine."""
        return pm - L.prof(c)[scored]

    def allmeans(n, r):
        """matched-depth centring means: one random n/n split of every target."""
        p = {t: r.permutation(400) for t in targets}
        SA = np.stack([dense[t][p[t][:n]].sum(0) for t in targets])
        SB = np.stack([dense[t][p[t][n:2 * n]].sum(0) for t in targets])
        return L.prof(SA)[:, scored].mean(0), L.prof(SB)[:, scored].mean(0)

    # sanity: the two centrings on the FULL 85-target ladder row
    pfix = {t: np.random.default_rng(L.SEED + 3).permutation(400) for t in targets}
    SAx = np.stack([dense[t][pfix[t][:100]].sum(0) for t in targets])
    SBx = np.stack([dense[t][pfix[t][100:200]].sum(0) for t in targets])
    v_glob = L.views(SAx, SBx, ctrl_full, ctrl_full, scored, mA=gm(ctrl_full), mB=gm(ctrl_full))
    v_own = L.views(SAx, SBx, ctrl_full, ctrl_full, scored)
    out["centring_check_n100"] = dict(side_own_median=float(np.nanmedian(v_own["cos_resid"])),
                                      global_panel_mean_median=float(np.nanmedian(v_glob["cos_resid"])),
                                      spearman=float(np.corrcoef(v_own["cos_resid"], v_glob["cos_resid"])[0, 1]))
    P(f"\ncentring check at n=100: side-own-mean median {out['centring_check_n100']['side_own_median']:+.4f} vs "
      f"full-panel-mean median {out['centring_check_n100']['global_panel_mean_median']:+.4f} "
      f"(r {out['centring_check_n100']['spearman']:+.3f})")

    # ------------------------------------------------------------------ (d) lanes at matched n
    P("\n(d) lane halves at matched cells per side (8 draws; targets with < n cells in a parity half are dropped,\n"
      "    and the rows kept are centred on the full-panel mean so the cohort change cannot inflate the cosine):")
    out["lanes"] = {}
    for n in (25, 50, 100, 200):
        acc = np.full((L.NDRAW, targets.size), np.nan)
        for d in range(L.NDRAW):
            r = np.random.default_rng(L.SEED + 31 * d + n)
            pmA, pmB = allmeans(n, r)
            keep, SA, SB = [], [], []
            for j, t in enumerate(targets):
                a = np.flatnonzero(lane_par[idx_t[t]] == 0); b = np.flatnonzero(lane_par[idx_t[t]] == 1)
                if len(a) < n or len(b) < n:
                    continue
                keep.append(j)
                SA.append(dense[t][r.choice(a, n, replace=False)].sum(0))
                SB.append(dense[t][r.choice(b, n, replace=False)].sum(0))
            if len(keep) < 3:
                continue
            v = L.views(np.stack(SA), np.stack(SB), cA, cB, scored, mA=cv_(pmA, cA), mB=cv_(pmB, cB))
            acc[d, np.asarray(keep)] = v["cos_resid"]
        x = np.nanmean(acc, 0)
        arrays[f"lane|{n}"] = x
        out["lanes"][n] = dict(**L.stats(x), boot=L.boot_median(x))
        if out["lanes"][n]["n"]:
            P(f"  n={n:4d}/side  median {L.fmt(out['lanes'][n]['boot'])}  n_targets {out['lanes'][n]['n']}  "
              f"<=0 {out['lanes'][n]['share_le0']:.3f}")

    # the matched (a) row under the SAME centring, for the lane contrast
    out["cells_globalcentre"] = {}
    for n in (25, 50, 100):
        acc = np.full((L.NDRAW, targets.size), np.nan)
        for d in range(L.NDRAW):
            r = np.random.default_rng(L.SEED + 131 * d + n)
            p = {t: r.permutation(400) for t in targets}
            SA = np.stack([dense[t][p[t][:n]].sum(0) for t in targets])
            SB = np.stack([dense[t][p[t][n:2 * n]].sum(0) for t in targets])
            pmA = L.prof(SA)[:, scored].mean(0); pmB = L.prof(SB)[:, scored].mean(0)
            acc[d] = L.views(SA, SB, ctrl_full, ctrl_full, scored,
                             mA=cv_(pmA, ctrl_full), mB=cv_(pmB, ctrl_full))["cos_resid"]
        x = np.nanmean(acc, 0)
        arrays[f"cells_gc|{n}"] = x
        out["cells_globalcentre"][n] = dict(**L.stats(x), boot=L.boot_median(x),
                                            vs_lane=L.boot_median_diff(x, arrays[f"lane|{n}"]))
        P(f"  (a) at n={n:4d} under the same centring: median {L.fmt(out['cells_globalcentre'][n]['boot'])}; "
          f"CELLS − LANE {L.fmt(out['cells_globalcentre'][n]['vs_lane'])}")

    # ------------------------------------------------------------------ (b) guide constructs, matched to (a)
    multi = {}
    for t in targets:
        gs = sorted(set(guide[idx_t[t]].tolist()))
        if len(gs) > 1:
            multi[t] = gs
    gt = sorted(multi)
    info = {}
    for t in gt:
        cs = [np.flatnonzero(guide[idx_t[t]] == g) for g in multi[t]]
        big, small = (cs[0], cs[1]) if len(cs[0]) >= len(cs[1]) else (cs[1], cs[0])
        info[t] = (big, small, int(min(len(big) // 2, len(small))))
    P(f"\n(b) guide construct: {len(multi)} two-construct targets; "
      + "; ".join(f"{t} {len(info[t][0])}/{len(info[t][1])} (n_t {info[t][2]})" for t in gt))
    out["guide"] = {}
    GN = [gn for gn in (37, 62) if sum(info[t][2] >= gn for t in gt) >= 3]
    for gn in GN:
        ts = [t for t in gt if info[t][2] >= gn]
        gacc = np.full((L.NDRAW, len(ts)), np.nan); aacc = np.full((L.NDRAW, len(ts)), np.nan)
        for d in range(L.NDRAW):
            r = np.random.default_rng(L.SEED + 97 * d + gn)
            pmA, pmB = allmeans(gn, r)
            SG1, SG2, SA1, SB1 = [], [], [], []
            for t in ts:
                big, small, _ = info[t]
                pb = r.permutation(big)
                SA1.append(dense[t][pb[:gn]].sum(0)); SB1.append(dense[t][pb[gn:2 * gn]].sum(0))
                SG1.append(dense[t][pb[:gn]].sum(0))
                SG2.append(dense[t][r.choice(small, gn, replace=False)].sum(0))
            gacc[d] = L.views(np.stack(SG1), np.stack(SG2), ctrl_full, ctrl_full, scored,
                              mA=cv_(pmA, ctrl_full), mB=cv_(pmB, ctrl_full))["cos_resid"]
            aacc[d] = L.views(np.stack(SA1), np.stack(SB1), ctrl_full, ctrl_full, scored,
                              mA=cv_(pmA, ctrl_full), mB=cv_(pmB, ctrl_full))["cos_resid"]
        xg, xa = np.nanmean(gacc, 0), np.nanmean(aacc, 0)
        arrays[f"guide|{gn}"] = xg; arrays[f"guide_cellctrl|{gn}"] = xa
        rat = xg / np.where(np.abs(xa) > 1e-9, xa, np.nan)
        out["guide"][gn] = dict(targets=ts, guide=L.stats(xg), cells=L.stats(xa),
                                diff=L.boot_median_diff(xa, xg), ratio_median=float(np.nanmedian(rat)),
                                per_target={t: dict(guide=float(xg[i]), cells=float(xa[i])) for i, t in enumerate(ts)})
        P(f"  n={gn}/side, {len(ts)} targets: GUIDE {L.fmt(L.boot_median(xg))}  CELLS in one construct "
          f"{L.fmt(L.boot_median(xa))}  CELLS - GUIDE {L.fmt(out['guide'][gn]['diff'])}  (guide/cells ratio "
          f"median {out['guide'][gn]['ratio_median']:+.3f})")
        for i, t in enumerate(ts):
            P(f"      {t:8s} guide {xg[i]:+.4f}  cells {xa[i]:+.4f}")

    # ------------------------------------------------------------------ leave-one-lane-out sensitivity on (a)
    P("\nleave-one-lane-out sensitivity of (a) at n=100 (drop each lane's cells, re-draw once):")
    loo = []
    for li, lname in enumerate(lanes):
        keep = {t: np.flatnonzero(lane[idx_t[t]] != lname) for t in targets}
        if min(len(v) for v in keep.values()) < 200:
            continue
        r = np.random.default_rng(L.SEED + 5 + li)
        pk = {t: r.permutation(keep[t]) for t in targets}          # ONE permutation, split into two disjoint sides
        SA = np.stack([dense[t][pk[t][:100]].sum(0) for t in targets])
        SB = np.stack([dense[t][pk[t][100:200]].sum(0) for t in targets])
        loo.append(float(np.nanmedian(score(SA, SB)["cos_resid"])))
    out["loo_lane_n100"] = dict(n=len(loo), min=min(loo) if loo else None, max=max(loo) if loo else None,
                               median=float(np.median(loo)) if loo else None)
    P(f"  {len(loo)} lanes droppable; median-of-medians {out['loo_lane_n100']['median']} "
      f"range [{out['loo_lane_n100']['min']}, {out['loo_lane_n100']['max']}]")

    # ------------------------------------------------------------------ save
    out["ka"] = ka
    out["failures"] = FAIL
    out["targets"] = targets.tolist()
    out["magnitude_decile"] = dec.tolist()
    out["n_lanes"] = len(lanes)
    out["multi_construct"] = {t: multi[t] for t in multi}
    out["seconds"] = round(time.time() - t0, 1)
    out["peak_rss_gib"] = L.peak_gib()
    np.savez_compressed(L.WORK / "h1_285.npz", targets=targets, magnitude=mag, decile=dec, **arrays)
    L.Logger.dump(P, out, L.WORK / "h1_285.json")
    P(f"\ndone: {len(FAIL)} KA failures, {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB")
    P.save()
    if FAIL:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
