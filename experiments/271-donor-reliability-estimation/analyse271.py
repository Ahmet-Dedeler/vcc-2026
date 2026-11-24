"""271 -- analysis on saved arrays: item 1 (donor-row reliability tables), item 3 (the estimators vs cd4svd5 with paired
target bootstraps, per fold and pooled, by 262 group and by reference-reliability tercile, the marginal-band crossing,
the controls), item 5 (the 400-cell extrapolation).  ARITHMETIC ON SAVED ARRAYS ONLY.

Inputs: work/sweep271-<fold>.json, work/pertarget271-<fold>.npz, work/identity271-<fold>.npz, work/half/rows-<line>.npz,
        experiments/267-.../work/half/truth-<line>.npz (the rand-split reliabilities, for the table only).
Outputs: work/analysis271.json, work/tables271.md, work/analyse271.log.
    PYTHONDONTWRITEBYTECODE=1 uv run python analyse271.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/271-donor-reliability-estimation"
WORK = HERE / "work"
E267 = ROOT / "experiments/267-replicate-validated-oracle/work/half"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t")
LINE = {"k562_gwps": "K562-GWPS", "orion_hct116": "HCT116", "orion_hek293t": "HEK293T"}
B = 20_000
SEED = 20260906
rng = np.random.default_rng(SEED)
ARMS = ["relw", "relw_eb", "relw_glob", "relw_wht", "ebrow", "ebrow_eb", "wiener_own", "wiener_gene", "a13:nm", "a13:lit",
        "relw:shuf", "relw:rev", "ebrow:shuf", "wiener_own:shuf", "random:all:s0:nm", "oracle_wht:all:nm", "eqR:nm", "STACK"]
MAIN = ["relw", "relw_eb", "relw_glob", "relw_wht", "ebrow", "ebrow_eb", "wiener_own", "wiener_gene", "a13:nm"]
CONTROL_OF = {"relw": "relw:shuf", "relw_eb": "relw:shuf", "relw_glob": "relw:shuf", "relw_wht": "relw:shuf", "ebrow": "ebrow:shuf",
              "ebrow_eb": "ebrow:shuf", "wiener_own": "wiener_own:shuf", "wiener_gene": "wiener_own:shuf"}
LINES: list[str] = []
MD: list[str] = []
OUT: dict = dict(seed=SEED, B=B, folds=list(FOLDS))


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def boot(diff: np.ndarray) -> dict:
    d = diff[np.isfinite(diff)]
    n = d.size
    idx = rng.integers(0, n, size=(B, n))
    means = d[idx].mean(1)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return dict(point=float(d.mean()), lo=float(lo), hi=float(hi), sd=float(means.std(ddof=1)), p_le0=float((means <= 0).mean()),
                p_ge0=float((means >= 0).mean()), n=int(n), n_nonzero=int((d != 0).sum()))


def boot_pool(diffs: list[np.ndarray]) -> dict:
    acc = np.zeros(B)
    for diff in diffs:
        d = diff[np.isfinite(diff)]
        idx = rng.integers(0, d.size, size=(B, d.size))
        acc += d[idx].mean(1)
    means = acc / len(diffs)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return dict(point=float(np.mean([np.nanmean(d) for d in diffs])), lo=float(lo), hi=float(hi), sd=float(means.std(ddof=1)),
                p_le0=float((means <= 0).mean()), p_ge0=float((means >= 0).mean()), per_fold=[float(np.nanmean(d)) for d in diffs],
                n=[int(np.isfinite(d).sum()) for d in diffs])


def fb(b: dict) -> str:
    return f"{b['point']:+.4f} [{b['lo']:+.4f}, {b['hi']:+.4f}]"


def q(x, p):
    return float(np.nanpercentile(x, p))


# ------------------------------------------------------------------------------------------------- load
D = {}
for f in FOLDS:
    sw = json.loads((WORK / f"sweep271-{f}.json").read_text())
    rows = {r["arm"]: r for r in sw["rows"]}
    per = {k: np.asarray(v) for k, v in np.load(WORK / f"pertarget271-{f}.npz", allow_pickle=True).items()}
    idn = {k: np.asarray(v, dtype=np.float64) for k, v in np.load(WORK / f"identity271-{f}.npz").items()}
    targets = [str(t) for t in sw["meta"]["targets"]]
    base_cos = np.asarray(per["base262_cos_truth"], dtype=np.float64)
    group = np.where(base_cos <= 0, "e_anti", np.where(base_cos <= 0.05, "marginal", "saturated"))
    ref = np.asarray(per["ref_rel"], dtype=np.float64)
    cuts = np.nanpercentile(ref, [100 / 3, 200 / 3])
    terc = np.where(ref <= cuts[0], "T1 low", np.where(ref <= cuts[1], "T2 mid", "T3 high"))
    D[f] = dict(meta=sw["meta"], rows=rows, per=per, idn=idn, targets=targets, group=group, terc=terc, cuts=cuts.tolist(), donors=sw["meta"]["donors"])
    P(f"[{f}] {len(targets)} targets, donors {sw['meta']['donors']}, KA {sw['meta']['ka_pass']}/{sw['meta']['ka_total']} failures {sw['meta']['failures']}, "
      f"{sw['meta']['n_rows']} rows, {sw['meta']['seconds']}s, peak {sw['meta']['peak_rss_gib']:.2f} GiB; est271 {sw['meta']['est271_sha256'][:12]}")
    OUT.setdefault("ka", {})[f] = sw["meta"]["ka"]
    OUT.setdefault("est271_sha256", sw["meta"]["est271_sha256"])
    assert sw["meta"]["est271_sha256"] == OUT["est271_sha256"], "est271.py changed between folds"

# ------------------------------------------------------------------------------------------------- item 1: the donor-row reliability table
P("\n== item 1: reliability of each DONOR's rows at its own depth (batch-parity halves; SB = Spearman-Brown to the full row) ==")
MD.append("## D1. Donor-row reliability: cos(row_A, row_B) per target, un-whitened residual space (`raw`) and after the fold's operator (`wht`)\n")
MD.append("| fold (held) | donor | n | cells/half (median) | cells full (median) | raw mean | raw p10 | raw p25 | **raw median** | raw p75 | raw p90 | **share ≤ 0** | share ≥ 0.05 | **SB full-depth median** | SB p10 / p90 | wht median | wht share ≤ 0 | wht p10 / p90 | ρ(r, cells) | raw by 262 group e / marg / sat | raw by ref-reliability tercile T1 / T2 / T3 | EB κ (r̄, s_null) |")
MD.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---|---:|---|---|---|")
OUT["reliability"] = {}
for f in FOLDS:
    d = D[f]
    per = d["per"]
    for j, line in enumerate(d["donors"]):
        ok = per["avail"][j].astype(bool) & np.isfinite(per["r_raw"][j])
        r = per["r_raw"][j]
        w = per["r_wht"][j]
        R = per["R_raw"][j]
        nc_half = np.minimum(per["ncells_A"][j], per["ncells_B"][j]).astype(float)
        ncf = per["ncells_full"][j].astype(float)
        rho = float(spearmanr(r[ok], ncf[ok]).statistic)
        by_g = {g: float(np.nanmean(r[ok & (d["group"] == g)])) for g in ("e_anti", "marginal", "saturated")}
        by_t = {g: float(np.nanmean(r[ok & (d["terc"] == g)])) for g in ("T1 low", "T2 mid", "T3 high")}
        eb = d["meta"]["eb"][line]
        rec = dict(fold=f, donor=line, n=int(ok.sum()), cells_half_median=float(np.median(nc_half[ok])), cells_full_median=float(np.median(ncf[ok])),
                   raw=dict(mean=float(np.nanmean(r[ok])), p10=q(r[ok], 10), p25=q(r[ok], 25), median=q(r[ok], 50), p75=q(r[ok], 75), p90=q(r[ok], 90),
                            share_le0=float(np.mean(r[ok] <= 0)), share_ge005=float(np.mean(r[ok] >= 0.05))),
                   sb=dict(median=q(R[ok], 50), p10=q(R[ok], 10), p90=q(R[ok], 90), mean=float(np.nanmean(R[ok]))),
                   wht=dict(median=q(w[ok], 50), p10=q(w[ok], 10), p90=q(w[ok], 90), share_le0=float(np.mean(w[ok] <= 0)), mean=float(np.nanmean(w[ok]))),
                   rho_cells=rho, by_group=by_g, by_ref_tercile=by_t, eb=eb,
                   rho_raw_wht=float(spearmanr(r[ok], w[ok]).statistic), rho_raw_ref=float(spearmanr(r[ok], per["ref_rel"][ok], nan_policy="omit").statistic))
        OUT["reliability"][f"{f}|{line}"] = rec
        MD.append(f"| `{f}` | `{line}` | {rec['n']} | {rec['cells_half_median']:.0f} | {rec['cells_full_median']:.0f} | {rec['raw']['mean']:+.4f} | {rec['raw']['p10']:+.4f} | "
                  f"{rec['raw']['p25']:+.4f} | **{rec['raw']['median']:+.4f}** | {rec['raw']['p75']:+.4f} | {rec['raw']['p90']:+.4f} | **{rec['raw']['share_le0']:.3f}** | "
                  f"{rec['raw']['share_ge005']:.3f} | **{rec['sb']['median']:+.4f}** | {rec['sb']['p10']:+.4f} / {rec['sb']['p90']:+.4f} | {rec['wht']['median']:+.4f} | "
                  f"{rec['wht']['share_le0']:.3f} | {rec['wht']['p10']:+.4f} / {rec['wht']['p90']:+.4f} | {rho:+.2f} | {by_g['e_anti']:+.4f} / {by_g['marginal']:+.4f} / {by_g['saturated']:+.4f} | "
                  f"{by_t['T1 low']:+.4f} / {by_t['T2 mid']:+.4f} / {by_t['T3 high']:+.4f} | {eb['kappa']:.2f} ({eb['rbar']:+.3f}, {np.sqrt(eb['s2_null']):.3f}) |")
        P(f"  {f:14s} donor {line:14s} n {rec['n']:3d} cells/half {rec['cells_half_median']:.0f}: raw median {rec['raw']['median']:+.4f} share<=0 {rec['raw']['share_le0']:.3f} "
          f"SB median {rec['sb']['median']:+.4f}; wht median {rec['wht']['median']:+.4f} share<=0 {rec['wht']['share_le0']:.3f} (rho raw~wht {rec['rho_raw_wht']:+.2f}); "
          f"rho(cells) {rho:+.2f}; group e/m/s {by_g['e_anti']:+.4f}/{by_g['marginal']:+.4f}/{by_g['saturated']:+.4f}; rho(r, ref) {rec['rho_raw_ref']:+.2f}; kappa {eb['kappa']:.2f}")
    # between the two donors of the fold: agreement of their reliabilities and whether reliability knows which donor is closer
    ok2 = per["avail"].astype(bool).all(0) & np.isfinite(per["r_raw"]).all(0)
    r0, r1 = per["r_raw"][0][ok2], per["r_raw"][1][ok2]
    R0, R1 = per["R_raw"][0][ok2], per["R_raw"][1][ok2]
    c0, c1 = per["cos_truth_wht"][:, 0][ok2], per["cos_truth_wht"][:, 1][ok2]
    rel_ratio = np.maximum(np.clip(R0, 0.005, 0.99), np.clip(R1, 0.005, 0.99)) / np.minimum(np.clip(R0, 0.005, 0.99), np.clip(R1, 0.005, 0.99))
    W = per["W|relw"][ok2]
    rec2 = dict(n=int(ok2.sum()), rho_r_donor0_donor1=float(spearmanr(r0, r1).statistic),
                rho_Rdiff_truthdiff=float(spearmanr(R0 - R1, c0 - c1).statistic),
                pick_by_R_matches_oracle=float(np.mean(np.sign(R0 - R1) == np.sign(c0 - c1))),
                rho_R_costruth=[float(spearmanr(R0, c0).statistic), float(spearmanr(R1, c1).statistic)],
                share_R_ratio_gt3=float(np.mean(rel_ratio > 3)), median_max_weight=float(np.median(W.max(1))), share_max_weight_gt09=float(np.mean(W.max(1) > 0.9)))
    OUT["reliability"][f"{f}|pair"] = rec2
    P(f"  {f}: rho(r donor0, r donor1) {rec2['rho_r_donor0_donor1']:+.3f}; rho(R0-R1, costruth0-costruth1) {rec2['rho_Rdiff_truthdiff']:+.3f}; "
      f"pick-by-R == oracle pick {100*rec2['pick_by_R_matches_oracle']:.0f}%; rho(R_j, cos_truth_j) {rec2['rho_R_costruth']}; R ratio > 3 on {100*rec2['share_R_ratio_gt3']:.0f}% of "
      f"two-donor targets; relw max weight median {rec2['median_max_weight']:.3f}, > 0.9 on {100*rec2['share_max_weight_gt09']:.0f}%")

# the rand-split reliabilities (267) per line, for the table only
MD.append("\nRandom-cell halves (267's `rel|rand|cos_resid`, the same lines; the batch split is the conservative reading used above):\n")
MD.append("| line | n | rand median | rand share ≤ 0 | batch median (267) | batch share ≤ 0 |")
MD.append("|---|---:|---:|---:|---:|---:|")
for line in FOLDS:
    tr = np.load(E267 / f"truth-{line}.npz", allow_pickle=True)
    rr = np.asarray(tr["rel|rand|cos_resid"], dtype=np.float64)
    rb = np.asarray(tr["rel|batch|cos_resid"], dtype=np.float64)
    MD.append(f"| `{line}` | {int(np.isfinite(rr).sum())} | {np.nanmedian(rr):+.4f} | {np.nanmean(rr <= 0):.3f} | {np.nanmedian(rb):+.4f} | {np.nanmean(rb <= 0):.3f} |")

# ------------------------------------------------------------------------------------------------- item 3: the measurement
P("\n== item 3: the estimators vs cd4svd5 (a23 pds on the cache truth; paired target bootstraps, 20,000, stratified by fold) ==")
MD.append("\n## M1. Every arm vs `cd4svd5`, three cancer folds (Δ`pds` per fold, pooled paired target bootstrap; cosine and identity legs)\n")
MD.append("| arm | Δ`pds` K562 / HCT116 / HEK293T | **pooled [95 % CI]** | P(≤ 0) | Δ`cos_agg` % (mean of folds) | Δ macro cos % | delivered ratio | Δ`nmae` identity (3 folds concat) | Δ top-89 precision | rows changed |")
MD.append("|---|---|---|---:|---:|---:|---:|---:|---:|---|")
OUT["arms"] = {}


def nmae_concat(f_list, arm):
    U = np.concatenate([D[f]["idn"][f"{arm}|U"] for f in f_list])
    M = np.concatenate([D[f]["idn"][f"{arm}|M"] for f in f_list])
    Dd = np.concatenate([D[f]["idn"][f"{arm}|D"] for f in f_list])
    return float(np.mean(1 + (U - 2 * M) / Dd))


nm_base = nmae_concat(FOLDS, "cd4svd5")
for arm in ARMS:
    if not all(arm in D[f]["rows"] for f in FOLDS):
        continue
    diffs = [D[f]["per"][f"{arm}|pds"] - D[f]["per"]["cd4svd5|pds"] for f in FOLDS]
    bp = boot_pool(diffs)
    per_fold = {f: boot(dd) for f, dd in zip(FOLDS, diffs)}
    ra = [D[f]["rows"][arm]["cos_agg"] / D[f]["rows"]["cd4svd5"]["cos_agg"] - 1 for f in FOLDS]
    rm = [D[f]["rows"][arm]["a23_cos"] / D[f]["rows"]["cd4svd5"]["a23_cos"] - 1 for f in FOLDS]
    dn = np.mean([D[f]["rows"][arm]["delivered_norm"] for f in FOLDS]) / np.mean([D[f]["rows"]["cd4svd5"]["delivered_norm"] for f in FOLDS])
    dnm = nmae_concat(FOLDS, arm) - nm_base
    dprec = float(np.mean([D[f]["rows"][arm]["top89_precision"] - D[f]["rows"]["cd4svd5"]["top89_precision"] for f in FOLDS]))
    chg = [D[f]["rows"][arm]["n_rows_changed"] for f in FOLDS]
    rec = dict(pooled=bp, per_fold=per_fold, rel_agg=[float(x) for x in ra], rel_macro=[float(x) for x in rm], delivered_ratio=float(dn), d_nmae=dnm,
               d_prec89=dprec, rows_changed=chg)
    OUT["arms"][arm] = rec
    MD.append(f"| `{arm}` | {' / '.join(f'{x:+.4f}' for x in bp['per_fold'])} | **{fb(bp)}** | {bp['p_le0']:.4f} | {100*np.mean(ra):+.3f} | {100*np.mean(rm):+.3f} | "
              f"{dn:.4f} | {dnm:+.5f} | {dprec:+.4f} | {'/'.join(map(str, chg))} |")
    P(f"  {arm:18s} d pds {' / '.join(f'{x:+.4f}' for x in bp['per_fold'])}  pooled {fb(bp)} P<=0 {bp['p_le0']:.3f}  agg {100*np.mean(ra):+.3f}%  macro {100*np.mean(rm):+.3f}%  "
      f"dn {dn:.4f}  d nmae {dnm:+.5f}  chg {chg}")

# arm minus its controls (paired)
MD.append("\n## M2. Each estimator against its controls (paired target bootstrap of the per-target difference, pooled over folds)\n")
MD.append("| arm | vs its shuffled control | vs `relw:rev` | vs `a13:nm` (the baseline (d)) | vs `random:all:s0:nm` | vs `relw_glob` (scalar weights) |")
MD.append("|---|---|---|---|---|---|")
OUT["vs_controls"] = {}
for arm in MAIN:
    rec = {}
    for other in (CONTROL_OF.get(arm, "relw:shuf"), "relw:rev", "a13:nm", "random:all:s0:nm", "relw_glob"):
        if other == arm:
            rec[other] = None
            continue
        diffs = [D[f]["per"][f"{arm}|pds"] - D[f]["per"][f"{other}|pds"] for f in FOLDS]
        rec[other] = boot_pool(diffs)
    OUT["vs_controls"][arm] = rec
    cells = [fb(rec[o]) if rec[o] else "—" for o in (CONTROL_OF.get(arm, "relw:shuf"), "relw:rev", "a13:nm", "random:all:s0:nm", "relw_glob")]
    MD.append(f"| `{arm}` | " + " | ".join(cells) + " |")
    P(f"  {arm:18s} " + "  ".join(f"vs {o}: {fb(rec[o])}" for o in rec if rec[o]))

# by 262 group and by reference-reliability tercile
MD.append("\n## M3. Δ`pds` vs `cd4svd5` by 262 group (base·truth cosine of the whitened mean) and by the held reference's split-half reliability tercile (267, batch halves); pooled over folds, mean per target [95 % CI]\n")
strata = [("group (e) base·truth ≤ 0", "group", "e_anti"), ("marginal band (0, 0.05]", "group", "marginal"), ("saturated > 0.05", "group", "saturated"),
          ("ref-reliability T1 (low)", "terc", "T1 low"), ("ref-reliability T2 (mid)", "terc", "T2 mid"), ("ref-reliability T3 (high)", "terc", "T3 high")]
MD.append("| arm | " + " | ".join(f"{s[0]} (n {sum(int((D[f][s[1]] == s[2]).sum()) for f in FOLDS)})" for s in strata) + " |")
MD.append("|---|" + "---|" * len(strata))
OUT["strata"] = {}
for arm in MAIN + ["relw:shuf", "random:all:s0:nm", "oracle_wht:all:nm"]:
    rec = {}
    cells = []
    for label, key, val in strata:
        diffs = [np.where(D[f][key] == val, D[f]["per"][f"{arm}|pds"] - D[f]["per"]["cd4svd5|pds"], np.nan) for f in FOLDS]
        b = boot_pool(diffs)
        rec[label] = b
        cells.append(fb(b))
    OUT["strata"][arm] = rec
    MD.append(f"| `{arm}` | " + " | ".join(cells) + " |")
    P(f"  {arm:18s} " + "  ".join(f"{s[2]} {rec[s[0]]['point']:+.4f}" for s in strata))

# the marginal band: crossing 0.02 in own cosine
MD.append("\n## M4. The marginal band (cd4svd5's own cosine in (0, 0.05]): own-cosine change and the count of targets crossing 0.02 (up − down), pooled\n")
MD.append("| arm | band n | mean Δ own cosine on the band | up across 0.02 | down across 0.02 | net share | mean Δ`pds` on the band [95 % CI] | Δ own cosine, all targets |")
MD.append("|---|---:|---:|---:|---:|---:|---|---:|")
OUT["band"] = {}
for arm in MAIN + ["relw:shuf", "relw:rev", "random:all:s0:nm", "oracle_wht:all:nm"]:
    nb = up = down = 0
    dc, dpds, dc_all = [], [], []
    for f in FOLDS:
        c0 = D[f]["per"]["cd4svd5|cos"]
        c1 = D[f]["per"][f"{arm}|cos"]
        band = (c0 > 0) & (c0 <= 0.05)
        nb += int(band.sum())
        up += int((band & (c0 < 0.02) & (c1 >= 0.02)).sum())
        down += int((band & (c0 >= 0.02) & (c1 < 0.02)).sum())
        dc.append(np.where(band, c1 - c0, np.nan))
        dc_all.append(c1 - c0)
        dpds.append(np.where(band, D[f]["per"][f"{arm}|pds"] - D[f]["per"]["cd4svd5|pds"], np.nan))
    b = boot_pool(dpds)
    rec = dict(n=nb, up=up, down=down, net_share=(up - down) / nb, d_cos_band=float(np.nanmean(np.concatenate(dc))), d_pds_band=b,
               d_cos_all=float(np.nanmean(np.concatenate(dc_all))))
    OUT["band"][arm] = rec
    MD.append(f"| `{arm}` | {nb} | {rec['d_cos_band']:+.5f} | {up} | {down} | {100*rec['net_share']:+.1f} % | {fb(b)} | {rec['d_cos_all']:+.5f} |")
    P(f"  {arm:18s} band n {nb}: d own cos {rec['d_cos_band']:+.5f}, up {up} down {down} net {100*rec['net_share']:+.1f}%, d pds on band {fb(b)}")

# leave-top-3-out retention (261's statistic) for the main arms
MD.append("\n## M5. Leave-top-3-out retention of the per-target Δ`pds` (261's statistic; each fold on its own ordering) and the top-3 share\n")
MD.append("| arm | full | leave-top-3-out | retained | per fold | top-3 share |")
MD.append("|---|---:|---:|---:|---|---:|")
OUT["retention"] = {}
for arm in MAIN + ["relw:shuf", "random:all:s0:nm", "oracle_wht:all:nm"]:
    full, loo, share = [], [], []
    for f in FOLDS:
        d = D[f]["per"][f"{arm}|pds"] - D[f]["per"]["cd4svd5|pds"]
        order = np.argsort(-d)
        keep = np.ones(d.size, bool)
        keep[order[:3]] = False
        full.append(d.mean())
        loo.append(d[keep].mean())
        share.append(float(d[order[:3]].sum() / d.sum()) if d.sum() != 0 else float("nan"))
    fm, lm = float(np.mean(full)), float(np.mean(loo))
    rec = dict(full=fm, leave_top_out=lm, retained=(lm / fm if fm else float("nan")), per_fold=[float(l / g) if g else float("nan") for l, g in zip(loo, full)],
               share=float(np.nanmean(share)))
    OUT["retention"][arm] = rec
    MD.append(f"| `{arm}` | {fm:+.5f} | {lm:+.5f} | {100*rec['retained']:+.0f} % | {' / '.join(f'{100*x:+.0f} %' for x in rec['per_fold'])} | {100*rec['share']:+.0f} % |")

# ------------------------------------------------------------------------------------------------- item 5: the 400-cell question
P("\n== item 5: the 400-cell extrapolation (159's laws on 267's half-depth numbers) ==")
MD.append("\n## Q. The 400-cell question: reliability of the target-specific component vs cells per condition (batch halves; 159 §1.3's laws)\n")
MD.append("| line | median cells / half | r at half depth (median) | SB → full donor depth (2× cells) | full depth (cells) | power law p = 1.0 → 400 cells | p = 0.70 → 400 | √ρ attenuation at donor depth (SB) | √ρ at 400 (p = 1.0 / 0.70) |")
MD.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
OUT["cells400"] = {}
for line in FOLDS:
    tr = np.load(E267 / f"truth-{line}.npz", allow_pickle=True)
    r = np.asarray(tr["rel|batch|cos_resid"], dtype=np.float64)
    ncA, ncB = np.asarray(tr["rel|batch|ncells_A"], float), np.asarray(tr["rel|batch|ncells_B"], float)
    m = float(np.nanmedian(np.minimum(ncA, ncB)))
    rm = float(np.nanmedian(r))
    sb = 2 * rm / (1 + rm)
    def law(r_m, m_from, m_to, p):
        x = (1.0 / r_m - 1.0) * (m_from / m_to) ** p
        return 1.0 / (1.0 + x)
    r400_p1, r400_p07 = law(rm, m, 400, 1.0), law(rm, m, 400, 0.70)
    rec = dict(cells_half=m, r_half_median=rm, sb_full=sb, cells_full=2 * m, r400_p1=r400_p1, r400_p07=r400_p07, r_full_law_p1=law(rm, m, 2 * m, 1.0),
               att_full=float(np.sqrt(max(sb, 0))), att400_p1=float(np.sqrt(r400_p1)), att400_p07=float(np.sqrt(r400_p07)))
    OUT["cells400"][line] = rec
    MD.append(f"| `{line}` | {m:.0f} | {rm:+.4f} | {sb:+.4f} | {2*m:.0f} | {r400_p1:.3f} | {r400_p07:.3f} | {rec['att_full']:.3f} | {rec['att400_p1']:.3f} / {rec['att400_p07']:.3f} |")
    P(f"  {line}: half {m:.0f} cells r {rm:+.4f} -> SB full {sb:+.4f} (law p=1 {rec['r_full_law_p1']:+.4f}); 400 cells: {r400_p1:.3f} (p=1) / {r400_p07:.3f} (p=0.7); "
      f"sqrt-rho {rec['att_full']:.3f} -> {rec['att400_p1']:.3f} / {rec['att400_p07']:.3f}")

(WORK / "analysis271.json").write_text(json.dumps(OUT, indent=1, default=float))
(WORK / "tables271.md").write_text("\n".join(MD) + "\n")
(WORK / "analyse271.log").write_text("\n".join(LINES) + "\n")
P("\nwritten work/analysis271.json, work/tables271.md, work/analyse271.log")
