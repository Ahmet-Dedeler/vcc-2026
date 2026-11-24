"""262 -- the analysis: distribution, deciles, fold agreement, covariate ranking, deficit decomposition,
and the per-target h1union -> cd4svd5 step.  Reads `work/diag262-<fold>.{npz,json}` and the fold caches
(for the cross-line truth agreement covariate); writes `work/analysis262.json` and `work/tables262.md`.
Nothing is fitted; every number is a rank correlation, an AUC, a median, or a share."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
CANCER = FOLDS[:3]
LINE = {"k562_gwps": "K562", "orion_hct116": "HCT116", "orion_hek293t": "HEK293T", "kolf2": "KOLF2"}
ARM, BASEARM = "cd4svd5", "h1union"

D = {f: np.load(WORK / f"diag262-{f}.npz") for f in FOLDS}
M = {f: json.loads((WORK / f"diag262-{f}.json").read_text()) for f in FOLDS}
CACHE = {f: np.load(WORK / "plolo" / f"fold-{f}.npz", allow_pickle=True) for f in FOLDS}

out: dict = {}
md: list[str] = []


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    return f"{x:.{nd}f}"


def unit_rows(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


# ------------------------------------------------------------------ cross-line agreement of the TRUTHS
# For every pair of lines, the cosine between the two lines' raw truths for the same target on the
# intersection of their scored columns (panel columns are zero on both).  This covariate is defined
# entirely outside the pipeline: it is how much the target's response agrees between two real lines.
truths = {f: np.asarray(CACHE[f]["truth"], dtype=np.float64) for f in FOLDS}
scored = {f: np.asarray(CACHE[f]["scored"]) for f in FOLDS}
tlist = {f: [str(x) for x in CACHE[f]["targets"]] for f in FOLDS}
tpos = {f: {t: i for i, t in enumerate(tlist[f])} for f in FOLDS}
pair_cos: dict[tuple, dict] = {}
for i, a in enumerate(FOLDS):
    for b in FOLDS[i + 1:]:
        cols = np.intersect1d(scored[a], scored[b])
        common = [t for t in tlist[a] if t in tpos[b]]
        A = truths[a][[tpos[a][t] for t in common]][:, cols]
        B = truths[b][[tpos[b][t] for t in common]][:, cols]
        c = (unit_rows(A) * unit_rows(B)).sum(1)
        pair_cos[(a, b)] = dict(zip(common, c.tolist()))
        pair_cos[(b, a)] = pair_cos[(a, b)]
        out.setdefault("cross_line_truth_cosine", {})[f"{LINE[a]}-{LINE[b]}"] = dict(
            n=len(common), mean=float(c.mean()), median=float(np.median(c)),
            q10=float(np.quantile(c, 0.1)), q90=float(np.quantile(c, 0.9)),
            frac_below_0=float((c <= 0).mean()), frac_below_0_05=float((c <= 0.05).mean()))

# assemble per-fold covariate tables
COV_STRUCT = [  # (key, label, higher-means)
    ("cov|base_cos_truth", "base·truth cos (whitened donor base vs the truth, effect level)", "base points at truth"),
    ("agree_held_vs_donors", "cross-line truth agreement: held line vs the fold's donor lines (mean cos)", "transferable response"),
    ("agree_donor_donor", "cross-line truth agreement between the fold's donor lines (cos; observable at build time)", "donors agree"),
    ("cov|resid_cos_truth", "residual·truth cos (fr·CD4_svd5 + g·union vs the truth)", "residual points at truth"),
    ("cov|cd4_deflated_cos_truth", "deflated CD4 row·truth cos", "CD4 points at truth"),
    ("cov|union_cos_truth", "union row·truth cos (carried rows only)", "union points at truth"),
    ("cov|resid_cos_base", "residual·base cos (does the residual agree with the base?)", "agree"),
    ("cov|resid_norm", "residual row norm (scored cols)", "bigger residual"),
    ("cov|resid_share", "residual share of the row, ‖R‖/‖W+R‖", "more residual"),
    ("cov|cd4_norm", "CD4 row norm (complete rows; 0 if absent)", "bigger CD4 row"),
    ("cov|cd4_energy_removed", "CD4 row energy removed by the k=5 deflation", "more deflated"),
    ("cov|union_norm", "union row norm (0 if not carried)", "bigger union row"),
    ("cov|union_count", "union carriers (0/1/2)", "more H1 carriers"),
    ("cov|n_whitened", "whitened sources carrying the target (0/1/2; 3 on kolf2)", "more donors"),
    ("cov|pm_carried", "PerturbMatch-blended row (K562 source only)", "PM"),
    ("cov|cd4_complete", "CD4 complete (3/3 conditions)", "CD4 complete"),
    ("cov|cd4_count", "CD4 condition count (0–3)", "more CD4 conditions"),
    ("cov|base_norm", "whitened base row norm", "bigger base row"),
    ("cov|essential", "k562_essential flag", "essential"),
    (f"{ARM}|truth_norm", "truth ‖Δ‖ (scored cols)", "bigger true effect"),
    (f"{ARM}|truth_nn", "truth's nearest-neighbour cosine to another truth", "truth confusable"),
    (f"{ARM}|truth_mean", "truth's mean cosine to the other truths", "truth generic"),
    (f"{ARM}|pred_nn", "prediction's nearest-neighbour cosine to another prediction", "prediction confusable"),
    (f"{ARM}|pred_mean", "prediction's mean cosine to the other predictions (hubness)", "hub prediction"),
    (f"{ARM}|comp_med", "median cosine of the prediction to the OTHER truths", "prediction near everything"),
    (f"{ARM}|comp_max", "max cosine of the prediction to another truth", "a close wrong truth"),
    ("amp_ratio", "delivered ‖pred‖ / ‖truth‖ (scored cols)", "over-scaled"),
    ("render_delta", "render effect on pds: delivered − effect-level", "render helps"),
]
COV_GEOM = [
    (f"{ARM}|own", "own cosine (delivered prediction vs its truth)"),
    (f"{ARM}|own_vs_nn", "own cosine minus the best competitor cosine"),
    ("BASE|pds", "pds of the whitened base alone (no residual)"),
    ("CD4ONLY|pds", "pds of base + deflated CD4 alone (no union)"),
    ("UNIONONLY|pds", "pds of base + union alone (no CD4)"),
    (f"{BASEARM}|pds", "pds of h1union (the un-deflated stack)"),
    (f"{ARM}|effect_pds", "pds at effect level (un-rendered)"),
]


def fold_table(f):
    z = D[f]
    t = np.asarray([str(x) for x in z["targets"]])
    n = len(t)
    cov = {k: z[k].astype(np.float64) for k in z.files if k.startswith("cov|") or k.split("|")[-1] in
           ("pds", "own", "own_vs_nn", "truth_norm", "truth_nn", "truth_mean", "pred_nn", "pred_mean", "comp_med",
            "comp_max", "comp_mean", "pred_norm", "effect_pds")}
    cov["amp_ratio"] = z[f"{ARM}|pred_norm"] / np.maximum(z[f"{ARM}|truth_norm"], 1e-300)
    cov["render_delta"] = z[f"{ARM}|pds"] - z[f"{ARM}|effect_pds"]
    donors = M[f]["meta"]["sources"]
    held_vs = np.array([np.mean([pair_cos[(f, d)].get(x, np.nan) for d in donors]) for x in t])
    if len(donors) == 2:
        dd = np.array([pair_cos[(donors[0], donors[1])].get(x, np.nan) for x in t])
    else:
        dd = np.array([np.mean([pair_cos[(donors[i], donors[j])].get(x, np.nan)
                                for i in range(len(donors)) for j in range(i + 1, len(donors))]) for x in t])
    cov["agree_held_vs_donors"] = held_vs
    cov["agree_donor_donor"] = dd
    return t, cov


TAB = {f: fold_table(f) for f in FOLDS}


def spearman(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 10 or np.std(x[ok]) == 0 or np.std(y[ok]) == 0:
        return float("nan")
    return float(spearmanr(x[ok], y[ok]).correlation)


def auc(v, mask):
    ok = np.isfinite(v)
    a, b = v[mask & ok], v[~mask & ok]
    if a.size == 0 or b.size == 0:
        return float("nan")
    return float((a[:, None] > b[None, :]).mean() + 0.5 * (a[:, None] == b[None, :]).mean())


def ranks(v):
    ok = np.isfinite(v)
    r = np.full(v.shape, np.nan)
    r[ok] = (np.argsort(np.argsort(v[ok])) + 1) / ok.sum()
    return r


# ================================================================== 1. distribution, deciles, agreement
md.append("## T1. Per-target `pds` of `cd4svd5`, per fold\n")
md.append("| fold (held line) | n | mean `pds` | deficit 1−mean | q10 | q25 | median | q75 | at 1.0 | < 0.5 | < 0.9 | bottom-decile share of deficit | bottom-quartile share |")
md.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
dist = {}
for f in FOLDS:
    t, cov = TAB[f]
    p = cov[f"{ARM}|pds"]
    n = len(p)
    dec = int(np.ceil(0.1 * n))
    order = np.argsort(p, kind="stable")
    bottom, top = order[:dec], order[::-1][:dec]
    dfc = 1 - p
    q = np.quantile(p, [0.1, 0.25, 0.5, 0.75])
    dist[f] = dict(n=n, mean=float(p.mean()), deficit=float(dfc.mean()), q10=float(q[0]), q25=float(q[1]),
                   median=float(q[2]), q75=float(q[3]), at_one=int((p == 1).sum()), below_half=int((p < 0.5).sum()),
                   below_0_9=int((p < 0.9).sum()), bottom_decile=[str(x) for x in t[bottom]],
                   bottom_decile_pds=[float(x) for x in p[bottom]],
                   top_decile=[str(x) for x in t[top]], n_decile=dec,
                   bottom_decile_deficit_share=float(dfc[bottom].sum() / dfc.sum()),
                   bottom_quartile_deficit_share=float(dfc[order[:n // 4]].sum() / dfc.sum()),
                   bottom_half_deficit_share=float(dfc[order[:n // 2]].sum() / dfc.sum()))
    d = dist[f]
    md.append(f"| `{f}` ({LINE[f]}){' — apart' if f == 'kolf2' else ''} | {n} | {d['mean']:.4f} | {d['deficit']:.4f} | "
              f"{d['q10']:.3f} | {d['q25']:.3f} | {d['median']:.3f} | {d['q75']:.3f} | {d['at_one']} | {d['below_half']} | "
              f"{d['below_0_9']} | {100*d['bottom_decile_deficit_share']:.1f} % | {100*d['bottom_quartile_deficit_share']:.1f} % |")
out["distribution"] = dist

# fold agreement
agree = {}
common3 = set.intersection(*[set(TAB[f][0]) for f in CANCER])
sets = {f: set(dist[f]["bottom_decile"]) for f in FOLDS}
tops = {f: set(dist[f]["top_decile"]) for f in FOLDS}
md.append("\n## T2. Fold-to-fold agreement of the bottom decile (cancer folds; targets restricted to the 259 common to all three)\n")
md.append("| pair | \\|A\\| | \\|B\\| | overlap | expected by chance | Jaccard | Spearman of per-target `pds` (all common targets) |")
md.append("|---|---:|---:|---:|---:|---:|---:|")
for i, a in enumerate(CANCER):
    for b in CANCER[i + 1:]:
        A, B = sets[a] & common3, sets[b] & common3
        ta, ca = TAB[a]; tb, cb = TAB[b]
        cm = sorted(set(ta) & set(tb))
        ia = {x: j for j, x in enumerate(ta)}; ib = {x: j for j, x in enumerate(tb)}
        rho = spearman(ca[f"{ARM}|pds"][[ia[x] for x in cm]], cb[f"{ARM}|pds"][[ib[x] for x in cm]])
        agree[f"{a}|{b}"] = dict(A=len(A), B=len(B), overlap=len(A & B), expected=len(A) * len(B) / len(common3),
                                 jaccard=len(A & B) / len(A | B), spearman_pds=rho, n_common=len(cm))
        md.append(f"| {LINE[a]} vs {LINE[b]} | {len(A)} | {len(B)} | **{len(A & B)}** | {len(A)*len(B)/len(common3):.1f} | "
                  f"{len(A & B)/len(A | B):.3f} | {rho:+.3f} (n = {len(cm)}) |")
all3 = sorted(sets[CANCER[0]] & sets[CANCER[1]] & sets[CANCER[2]])
ge2 = sorted(x for x in common3 if sum(x in sets[f] for f in CANCER) >= 2)
top_all3 = sorted(tops[CANCER[0]] & tops[CANCER[1]] & tops[CANCER[2]])
agree["bottom_all_three"] = all3
agree["bottom_two_or_more"] = ge2
agree["kolf2_bottom_overlap_with_ge2"] = sorted(set(ge2) & sets["kolf2"])
md.append(f"\nBottom decile on **all three** cancer folds ({len(all3)}): {', '.join(f'`{x}`' for x in all3)}.  "
          f"On **two or more** ({len(ge2)}): {', '.join(f'`{x}`' for x in ge2)}.  Of those, in `kolf2`'s bottom decile: "
          f"{', '.join(f'`{x}`' for x in agree['kolf2_bottom_overlap_with_ge2']) or 'none'}.")
# per-target pds across folds for the >=2 set
md.append("\n| target (bottom decile on ≥ 2 cancer folds) | " + " | ".join(f"`pds` {LINE[f]}" for f in FOLDS) +
          " | base·truth cos (mean, cancer) | held-vs-donor truth agreement (mean, cancer) | donor-donor agreement (mean) | CD4 | union carriers | whitened sources (cancer folds) |")
md.append("|---|" + "---:|" * (len(FOLDS) + 6) + "---|")
rows_ge2 = {}
for x in ge2:
    vals, bc, hv, ddv = [], [], [], []
    cd4c = unc = None
    nw = []
    for f in FOLDS:
        t, cov = TAB[f]
        pos = {s: j for j, s in enumerate(t)}
        if x in pos:
            j = pos[x]
            vals.append(cov[f"{ARM}|pds"][j])
            if f != "kolf2":
                bc.append(cov["cov|base_cos_truth"][j]); hv.append(cov["agree_held_vs_donors"][j]); ddv.append(cov["agree_donor_donor"][j])
                nw.append(int(cov["cov|n_whitened"][j]))
            cd4c = int(cov["cov|cd4_count"][j]); unc = int(cov["cov|union_count"][j])
        else:
            vals.append(np.nan)
    rows_ge2[x] = dict(pds=[float(v) for v in vals], base_cos=float(np.nanmean(bc)), held_vs_donors=float(np.nanmean(hv)),
                       donor_donor=float(np.nanmean(ddv)), cd4_count=cd4c, union_count=unc, n_whitened=nw)
    md.append(f"| `{x}` | " + " | ".join(fmt(v, 2) if np.isfinite(v) else "—" for v in vals) +
              f" | {np.nanmean(bc):+.3f} | {np.nanmean(hv):+.3f} | {np.nanmean(ddv):+.3f} | {cd4c}/3 | {unc} | {'/'.join(map(str, nw))} |")
agree["rows_two_or_more"] = rows_ge2
agree["top_decile_all_three"] = top_all3
out["agreement"] = agree
md.append(f"\nTop decile (all at `pds` = 1.0 on every fold; the 1.0 block is {', '.join(str(dist[f]['at_one']) for f in CANCER)} targets on the "
          f"three cancer folds) shared by all three cancer folds: {len(top_all3)} of the {dist[CANCER[0]]['n_decile']}-target deciles — "
          f"{', '.join(f'`{x}`' for x in top_all3[:40])}{' …' if len(top_all3) > 40 else ''}.")

# ================================================================== 2. bottom / top decile tables with covariates
COLS = [("pds", f"{ARM}|pds", 2), ("Δ vs h1union", "delta", 2), ("base pds", "BASE|pds", 2), ("own cos", f"{ARM}|own", 3),
        ("base·truth", "cov|base_cos_truth", 3), ("resid·truth", "cov|resid_cos_truth", 3), ("held-vs-donor agree", "agree_held_vs_donors", 3),
        ("donor-donor agree", "agree_donor_donor", 3), ("‖truth‖", f"{ARM}|truth_norm", 1), ("truth NN", f"{ARM}|truth_nn", 2),
        ("pred NN", f"{ARM}|pred_nn", 2), ("comp med", f"{ARM}|comp_med", 3), ("src", "cov|n_whitened", 0), ("CD4", "cov|cd4_count", 0),
        ("union", "cov|union_count", 0), ("‖CD4‖", "cov|cd4_norm", 1), ("CD4 E rm", "cov|cd4_energy_removed", 2)]
for f in FOLDS:
    t, cov = TAB[f]
    cov["delta"] = cov[f"{ARM}|pds"] - cov[f"{BASEARM}|pds"]
    p = cov[f"{ARM}|pds"]
    n = len(p); dec = int(np.ceil(0.1 * n))
    order = np.argsort(p, kind="stable")
    for which, idx in (("bottom", order[:dec]), ("top", order[::-1][:dec])):
        md.append(f"\n## T3{'b' if which == 'bottom' else 't'}. `{f}` ({LINE[f]}) — {which} decile ({len(idx)} targets){' — apart' if f == 'kolf2' else ''}\n")
        md.append("| target | " + " | ".join(c[0] for c in COLS) + " |")
        md.append("|---|" + "---:|" * len(COLS))
        if which == "top":
            # the top decile is the 1.0 block: list names, then medians only
            md.append(f"| *{len(idx)} targets, all at 1.0:* {', '.join(f'`{x}`' for x in t[idx])} | " + " | ".join("" for _ in COLS) + " |")
            med = [np.nanmedian(cov[c[1]][idx]) for c in COLS]
            md.append("| **median, top decile** | " + " | ".join(fmt(v, c[2]) for v, c in zip(med, COLS)) + " |")
            rest = np.ones(n, bool); rest[idx] = False
            medr = [np.nanmedian(cov[c[1]][rest]) for c in COLS]
            md.append("| median, the rest | " + " | ".join(fmt(v, c[2]) for v, c in zip(medr, COLS)) + " |")
            continue
        for j in idx:
            md.append(f"| `{t[j]}` | " + " | ".join(fmt(float(cov[c[1]][j]), c[2]) for c in COLS) + " |")
        med = [np.nanmedian(cov[c[1]][idx]) for c in COLS]
        rest = np.ones(n, bool); rest[idx] = False
        medr = [np.nanmedian(cov[c[1]][rest]) for c in COLS]
        md.append("| **median, bottom decile** | " + " | ".join(fmt(v, c[2]) for v, c in zip(med, COLS)) + " |")
        md.append("| median, the rest | " + " | ".join(fmt(v, c[2]) for v, c in zip(medr, COLS)) + " |")

# ================================================================== 3. covariate ranking
rank = {}
md.append("\n## T4. Covariate ranking — structural covariates vs per-target `pds` of `cd4svd5` (three cancer folds; `kolf2` apart)\n")
md.append("Spearman ρ with `pds` per fold and pooled by within-fold rank; AUC = P(a bottom-decile target has the HIGHER covariate value than a random other target) — 0.5 is no separation, far from 0.5 is separation (below 0.5: the bottom decile has LOWER values); pooled median of the covariate in the bottom decile vs the rest.\n")
md.append("| covariate | ρ K562 | ρ HCT116 | ρ HEK293T | **ρ pooled** | AUC K562 | AUC HCT116 | AUC HEK293T | **AUC pooled** | bottom-decile median | rest median | ρ / AUC `kolf2` |")
md.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
pooled_p, pooled_bot = [], []
for f in CANCER:
    t, cov = TAB[f]
    p = cov[f"{ARM}|pds"]; dec = int(np.ceil(0.1 * len(p)))
    bot = np.zeros(len(p), bool); bot[np.argsort(p, kind="stable")[:dec]] = True
    pooled_p.append(ranks(p)); pooled_bot.append(bot)
pooled_p = np.concatenate(pooled_p); pooled_bot = np.concatenate(pooled_bot)
for key, label, hi in COV_STRUCT + [(k, l, "") for k, l in COV_GEOM]:
    rhos, aucs, vals = [], [], []
    for f in CANCER:
        t, cov = TAB[f]
        v = cov[key]; p = cov[f"{ARM}|pds"]
        dec = int(np.ceil(0.1 * len(p))); bot = np.zeros(len(p), bool); bot[np.argsort(p, kind="stable")[:dec]] = True
        rhos.append(spearman(v, p)); aucs.append(auc(v, bot)); vals.append(v)
    vr = np.concatenate([ranks(v) for v in vals]); vraw = np.concatenate(vals)
    rho_pool = spearman(vr, pooled_p); auc_pool = auc(vr, pooled_bot)
    med_b = float(np.nanmedian(vraw[pooled_bot])) if np.isfinite(vraw[pooled_bot]).any() else float("nan")
    med_r = float(np.nanmedian(vraw[~pooled_bot])) if np.isfinite(vraw[~pooled_bot]).any() else float("nan")
    tk, ck = TAB["kolf2"]
    pk = ck[f"{ARM}|pds"]; deck = int(np.ceil(0.1 * len(pk))); botk = np.zeros(len(pk), bool); botk[np.argsort(pk, kind="stable")[:deck]] = True
    rk, ak = spearman(ck[key], pk), auc(ck[key], botk)
    rank[key] = dict(label=label, rho=rhos, rho_pooled=rho_pool, auc=aucs, auc_pooled=auc_pool, median_bottom=med_b,
                     median_rest=med_r, kolf2=dict(rho=rk, auc=ak))
    nd = 1 if "norm" in key and "cos" not in key else 3
    md.append(f"| {label} | " + " | ".join(fmt(r, 3) if np.isfinite(r) else "—" for r in rhos) + f" | **{fmt(rho_pool, 3)}** | " +
              " | ".join(fmt(a, 3) if np.isfinite(a) else "—" for a in aucs) + f" | **{fmt(auc_pool, 3)}** | {fmt(med_b, nd)} | {fmt(med_r, nd)} | "
              f"{fmt(rk, 2)} / {fmt(ak, 2)} |")
out["covariate_ranking"] = rank

# ================================================================== 4. deficit decomposition
md.append("\n## T5. Deficit decomposition (ceiling: `pds` = 1 on every target, i.e. every prediction nearest its own truth; the fold deficit is 1 − mean `pds`)\n")
decomp = {}
md.append("### T5a. Where the beats are: prediction side vs truth side\n")
md.append("| fold | beats (prediction, truth, frame) | share from the 10 most-beaten predictions | bottom-decile predictions | bottom quartile | share from the 10 most-beating (hub) truths | top-30 hub truths | truths beating nothing | top-10 hub truths |")
md.append("|---|---:|---:|---:|---:|---:|---:|---:|---|")
for f in FOLDS:
    z = D[f]; t = np.asarray([str(x) for x in z["targets"]])
    B = z[f"{ARM}|beats"].astype(np.int64)
    n = len(t); dec = int(np.ceil(0.1 * n))
    hub = B.sum(0); rowb = B.sum(1)
    hs = np.sort(hub)[::-1]; rs = np.sort(rowb)[::-1]
    tot = B.sum()
    top_hubs = np.argsort(hub)[::-1][:10]
    decomp.setdefault(f, {})["beats"] = dict(total=int(tot), pred_top10=float(rs[:10].sum() / tot), pred_bottom_decile=float(rs[:dec].sum() / tot),
                                             pred_bottom_quartile=float(rs[:n // 4].sum() / tot), truth_top10=float(hs[:10].sum() / tot),
                                             truth_top30=float(hs[:30].sum() / tot), truth_zero=int((hub == 0).sum()),
                                             hub_truths=[(str(t[i]), int(hub[i])) for i in top_hubs])
    md.append(f"| `{f}` | {tot} | {100*rs[:10].sum()/tot:.1f} % | {100*rs[:dec].sum()/tot:.1f} % | {100*rs[:n//4].sum()/tot:.1f} % | "
              f"**{100*hs[:10].sum()/tot:.1f} %** | {100*hs[:30].sum()/tot:.1f} % | {(hub == 0).sum()} | "
              + ", ".join(f"`{t[i]}` {hub[i]}" for i in top_hubs[:6]) + " |")
md.append("\nA hub truth would beat many predictions; a uniformly spread loss puts every truth at ≈ 1/n of the beats (0.3–0.4 %). The top-10 hub truths take 4–5 % against a uniform 3.4–3.7 %: there is no truth-side hub left. The loss sits on the prediction side.\n")

# hub remaining: h1union vs cd4svd5
md.append("### T5b. How much of the residual-hub story remains (h1union → cd4svd5)\n")
md.append("| fold | mean pairwise prediction cosine (hubness) h1union → cd4svd5 | ρ(prediction hubness, `pds`) h1union → cd4svd5 | ρ(median cos to other truths, `pds`) h1union → cd4svd5 | top-10 hub-truth beat share h1union → cd4svd5 |")
md.append("|---|---:|---:|---:|---:|")
for f in FOLDS:
    z = D[f]
    hh = float(z[f"{BASEARM}|pred_mean"].mean()); ha = float(z[f"{ARM}|pred_mean"].mean())
    r1 = spearman(z[f"{BASEARM}|pred_mean"], z[f"{BASEARM}|pds"]); r2 = spearman(z[f"{ARM}|pred_mean"], z[f"{ARM}|pds"])
    c1 = spearman(z[f"{BASEARM}|comp_med"], z[f"{BASEARM}|pds"]); c2 = spearman(z[f"{ARM}|comp_med"], z[f"{ARM}|pds"])
    Bh = z[f"{BASEARM}|beats"].astype(np.int64); Ba = z[f"{ARM}|beats"].astype(np.int64)
    sh = float(np.sort(Bh.sum(0))[::-1][:10].sum() / Bh.sum()); sa = float(np.sort(Ba.sum(0))[::-1][:10].sum() / Ba.sum())
    decomp[f]["hub_remaining"] = dict(pred_mean_h1union=hh, pred_mean_cd4svd5=ha, rho_predmean_h1union=r1, rho_predmean_cd4svd5=r2,
                                      rho_compmed_h1union=c1, rho_compmed_cd4svd5=c2, hub10_h1union=sh, hub10_cd4svd5=sa)
    md.append(f"| `{f}` | {hh:.4f} → {ha:.4f} | {r1:+.3f} → {r2:+.3f} | {c1:+.3f} → {c2:+.3f} | {100*sh:.1f} % → {100*sa:.1f} % |")

# the partition
md.append("\n### T5c. The partition of the deficit by mechanism (each target assigned to the FIRST group it matches, in this order)\n")
md.append("Groups: **(b) no base support** — carried by none of the fold's whitened sources (the base row is the mean-response filler); "
          "**(e) non-transferable direction** — the whitened donor base's cosine to the truth ≤ 0 (the donors' consensus points away from this line's response); "
          "**(r) residual breaks it** — base·truth > 0 but `pds` fell below the base-only `pds` by more than 0.1 after the residual; "
          "**(a) hub / confusable** — the prediction's median cosine to the other truths is in the fold's top quintile; "
          "**(c) small true effect** — truth ‖Δ‖ in the fold's bottom quintile; **(d) render / magnitude** — the render moved `pds` by more than 0.05; "
          "**(o) other** — none of the above (the prediction points the right way but a few specific truths sit closer).\n")
md.append("| fold | deficit | (b) no base | (e) base points away | (r) residual breaks it | (a) hub | (c) small truth | (d) render | (o) other: specific confusion | n targets per group |")
md.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---|")
GROUPS = ["b_no_base", "e_base_away", "r_resid_breaks", "a_hub", "c_small_truth", "d_render", "o_specific"]
for f in FOLDS:
    t, cov = TAB[f]
    p = cov[f"{ARM}|pds"]; dfc = 1 - p; n = len(p)
    bc = cov["cov|base_cos_truth"]; bp = cov["BASE|pds"]
    cm = cov[f"{ARM}|comp_med"]; tn = cov[f"{ARM}|truth_norm"]; rd = cov["render_delta"]
    g = np.full(n, "o_specific", dtype=object)
    conds = [
        ("b_no_base", cov["cov|n_whitened"] == 0),
        ("e_base_away", bc <= 0),
        ("r_resid_breaks", (bc > 0) & (p < bp - 0.1)),
        ("a_hub", cm >= np.quantile(cm, 0.8)),
        ("c_small_truth", tn <= np.quantile(tn, 0.2)),
        ("d_render", np.abs(rd) > 0.05),
    ]
    assigned = np.zeros(n, bool)
    for name, c in conds:
        pick = c & ~assigned & (dfc > 0)
        g[pick] = name; assigned |= pick
    shares = {k: float(dfc[g == k].sum() / dfc.sum()) for k in GROUPS}
    counts = {k: int((g == k).sum()) for k in GROUPS}
    counts["zero_deficit"] = int((dfc == 0).sum())
    # marginal (overlapping) shares too
    marginal = {name: float(dfc[c].sum() / dfc.sum()) for name, c in conds}
    marginal_n = {name: int(c.sum()) for name, c in conds}
    decomp[f]["partition"] = dict(shares=shares, counts=counts, marginal_shares=marginal, marginal_counts=marginal_n,
                                  deficit=float(dfc.mean()), group=g.tolist())
    md.append(f"| `{f}` | {dfc.mean():.4f} | " + " | ".join(f"{100*shares[k]:.1f} %" for k in GROUPS) + " | " +
              ", ".join(f"{k.split('_')[0]} {counts[k]}" for k in GROUPS) + f", zero-deficit {counts['zero_deficit']} |")
md.append("\nMarginal (overlapping) shares of the deficit, i.e. each condition on its own:\n")
md.append("| fold | (b) no base | (e) base·truth ≤ 0 | (r) residual breaks | (a) hub quintile | (c) small-truth quintile | (d) render > 0.05 |")
md.append("|---|---:|---:|---:|---:|---:|---:|")
for f in FOLDS:
    mg = decomp[f]["partition"]["marginal_shares"]; mn = decomp[f]["partition"]["marginal_counts"]
    md.append(f"| `{f}` | " + " | ".join(f"{100*mg[k]:.1f} % (n {mn[k]})" for k, _ in conds) + " |")

# finer: base-cos bins
md.append("\n### T5d. The deficit by bins of base·truth cosine and of held-vs-donor truth agreement (cancer folds pooled)\n")
for key, label in (("cov|base_cos_truth", "base·truth cosine"), ("agree_held_vs_donors", "held-vs-donor truth agreement"),
                   ("agree_donor_donor", "donor-donor truth agreement (observable)")):
    v = np.concatenate([TAB[f][1][key] for f in CANCER]); p = np.concatenate([TAB[f][1][f"{ARM}|pds"] for f in CANCER])
    dfc = 1 - p
    edges = [-1, 0.0, 0.05, 0.1, 0.2, 1.01]
    md.append(f"\n| {label} bin | n targets | share of targets | mean `pds` | share of the pooled deficit | deficit per target |")
    md.append("|---|---:|---:|---:|---:|---:|")
    bins = {}
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (v > lo) & (v <= hi) & np.isfinite(v)
        if m.sum() == 0:
            continue
        bins[f"({lo}, {hi}]"] = dict(n=int(m.sum()), mean_pds=float(p[m].mean()), deficit_share=float(dfc[m].sum() / dfc.sum()),
                                     deficit_per_target=float(dfc[m].mean()))
        md.append(f"| ({lo:g}, {hi:g}] | {m.sum()} | {100*m.mean():.1f} % | {p[m].mean():.3f} | **{100*dfc[m].sum()/dfc.sum():.1f} %** | {dfc[m].mean():.3f} |")
    decomp.setdefault("bins", {})[key] = bins
out["decomposition"] = decomp

# ================================================================== 5. the cd4svd5 step per target
step = {}
md.append("\n## T6. The `h1union` → `cd4svd5` step, per target\n")
md.append("| fold | mean Δ`pds` | gainers / losers / unchanged | sum of gains | sum of losses | top-3 share of the net gain | ρ(Δ, base·truth) | ρ(Δ, CD4 row energy removed) | ρ(Δ, CD4·truth before) | ρ(Δ, CD4 row norm) | ρ(Δ, `pds` before) | AUC losers-higher: E removed / ‖CD4‖ / base·truth / `pds` before |")
md.append("|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|")
for f in FOLDS:
    t, cov = TAB[f]
    dlt = cov["delta"]; p0 = cov[f"{BASEARM}|pds"]
    gain, loss = dlt[dlt > 0].sum(), dlt[dlt < 0].sum()
    losers = dlt < 0
    o = np.argsort(dlt)[::-1]
    top3 = float(dlt[o[:3]].sum() / dlt.sum()) if dlt.sum() != 0 else float("nan")
    r = {k: spearman(cov[k], dlt) for k in ("cov|base_cos_truth", "cov|cd4_energy_removed", "cov|cd4_cos_truth", "cov|cd4_norm", f"{BASEARM}|pds",
                                             "cov|cd4_deflated_cos_truth", "agree_held_vs_donors", "cov|union_count", f"{BASEARM}|pred_mean")}
    a = {k: auc(cov[k], losers) for k in ("cov|cd4_energy_removed", "cov|cd4_norm", "cov|base_cos_truth", f"{BASEARM}|pds")}
    step[f] = dict(mean=float(dlt.mean()), gainers=int((dlt > 0).sum()), losers=int(losers.sum()), same=int((dlt == 0).sum()),
                   sum_gain=float(gain), sum_loss=float(loss), top3_share=top3, rho=r, auc_losers=a,
                   top_losers=[(str(t[i]), float(dlt[i])) for i in np.argsort(dlt)[:10]],
                   top_gainers=[(str(t[i]), float(dlt[i])) for i in o[:10]])
    md.append(f"| `{f}` | {dlt.mean():+.5f} | {int((dlt > 0).sum())} / {int(losers.sum())} / {int((dlt == 0).sum())} | {gain:+.3f} | {loss:+.3f} | "
              f"{100*top3:.0f} % | {r['cov|base_cos_truth']:+.3f} | {r['cov|cd4_energy_removed']:+.3f} | {r['cov|cd4_cos_truth']:+.3f} | "
              f"{r['cov|cd4_norm']:+.3f} | {r[f'{BASEARM}|pds']:+.3f} | {a['cov|cd4_energy_removed']:.2f} / {a['cov|cd4_norm']:.2f} / "
              f"{a['cov|base_cos_truth']:.2f} / {a[f'{BASEARM}|pds']:.2f} |")
# fold consistency of the step
md.append("\nFold-to-fold consistency of the per-target step (Spearman of Δ`pds` across cancer folds, common targets): ")
cons = {}
for i, a_ in enumerate(CANCER):
    for b_ in CANCER[i + 1:]:
        ta, ca = TAB[a_]; tb, cb = TAB[b_]
        cm = sorted(set(ta) & set(tb)); ia = {x: j for j, x in enumerate(ta)}; ib = {x: j for j, x in enumerate(tb)}
        cons[f"{a_}|{b_}"] = spearman(ca["delta"][[ia[x] for x in cm]], cb["delta"][[ib[x] for x in cm]])
md[-1] += "; ".join(f"{LINE[k.split('|')[0]]}–{LINE[k.split('|')[1]]} {v:+.3f}" for k, v in cons.items()) + "."
step["consistency"] = cons
# biggest movers with per-fold values
movers = {}
for f in CANCER:
    for nm, v in step[f]["top_losers"][:5] + step[f]["top_gainers"][:5]:
        movers.setdefault(nm, {})
for nm in movers:
    for f in FOLDS:
        t, cov = TAB[f]; pos = {s: j for j, s in enumerate(t)}
        movers[nm][f] = float(cov["delta"][pos[nm]]) if nm in pos else float("nan")
step["movers"] = movers
md.append("\n| biggest movers (top-5 losers and gainers of each cancer fold) | " + " | ".join(f"Δ {LINE[f]}" for f in FOLDS) + " | same sign on all cancer folds? |")
md.append("|---|" + "---:|" * len(FOLDS) + "---|")
for nm in sorted(movers, key=lambda x: -np.nanmean([abs(movers[x][f]) for f in CANCER])):
    vv = [movers[nm][f] for f in FOLDS]
    s = [np.sign(movers[nm][f]) for f in CANCER if np.isfinite(movers[nm][f])]
    md.append(f"| `{nm}` | " + " | ".join(f"{v:+.3f}" if np.isfinite(v) else "—" for v in vv) + f" | {'yes' if len(set(s)) == 1 else 'no'} |")
# losers vs gainers covariate medians (pooled cancer)
md.append("\n| covariate (cancer folds pooled) | median, losers (Δ < 0) | median, gainers (Δ > 0) | median, unchanged | AUC losers-higher |")
md.append("|---|---:|---:|---:|---:|")
lg = {}
allc = {k: np.concatenate([TAB[f][1][k] for f in CANCER]) for k in ("delta", "cov|base_cos_truth", "cov|cd4_energy_removed", "cov|cd4_norm",
                                                                   "cov|cd4_cos_truth", "cov|cd4_deflated_cos_truth", "cov|cd4_cos_base",
                                                                   f"{BASEARM}|pds", f"{BASEARM}|pred_mean", "agree_held_vs_donors",
                                                                   "agree_donor_donor", "cov|union_count", "cov|n_whitened", f"{ARM}|truth_norm",
                                                                   "cov|resid_share")}
dl = allc["delta"]
for k in list(allc)[1:]:
    v = allc[k]
    lg[k] = dict(losers=float(np.nanmedian(v[dl < 0])), gainers=float(np.nanmedian(v[dl > 0])), same=float(np.nanmedian(v[dl == 0])),
                 auc_losers=auc(v, dl < 0))
    nd = 1 if ("norm" in k and "cos" not in k) else 3
    md.append(f"| {k} | {fmt(lg[k]['losers'], nd)} | {fmt(lg[k]['gainers'], nd)} | {fmt(lg[k]['same'], nd)} | {lg[k]['auc_losers']:.3f} |")
step["losers_vs_gainers"] = lg
out["step"] = step

# ================================================================== 6. component attribution per target
md.append("\n## T7. Component attribution of the loss (per fold): the same targets, four objects\n")
md.append("| fold | mean `pds`: BASE / CD4ONLY / UNIONONLY / h1union / cd4svd5 | ρ(BASE `pds`, cd4svd5 `pds`) | targets with BASE `pds` < 0.5 that cd4svd5 lifts above 0.5 | targets with BASE `pds` ≥ 0.9 that cd4svd5 drops below 0.5 | share of cd4svd5's deficit on targets already < 0.5 in BASE | render: mean `pds` delivered − effect-level (BASE / cd4svd5) |")
md.append("|---|---|---:|---:|---:|---:|---|")
attr = {}
for f in FOLDS:
    z = D[f]; t, cov = TAB[f]
    pb, pc, pu, ph, pa = z["BASE|pds"], z["CD4ONLY|pds"], z["UNIONONLY|pds"], z[f"{BASEARM}|pds"], z[f"{ARM}|pds"]
    lifted = int(((pb < 0.5) & (pa >= 0.5)).sum()); dropped = int(((pb >= 0.9) & (pa < 0.5)).sum())
    share = float((1 - pa)[pb < 0.5].sum() / (1 - pa).sum())
    rb = float(z["BASE|pds"].mean() - z["BASE|effect_pds"].mean()); ra = float(pa.mean() - z[f"{ARM}|effect_pds"].mean())
    attr[f] = dict(means=dict(BASE=float(pb.mean()), CD4ONLY=float(pc.mean()), UNIONONLY=float(pu.mean()), h1union=float(ph.mean()), cd4svd5=float(pa.mean())),
                   rho_base_arm=spearman(pb, pa), lifted=lifted, dropped=dropped, deficit_share_on_base_below_half=share,
                   render_base=rb, render_arm=ra)
    md.append(f"| `{f}` | {pb.mean():.4f} / {pc.mean():.4f} / {pu.mean():.4f} / {ph.mean():.4f} / {pa.mean():.4f} | {spearman(pb, pa):+.3f} | "
              f"{lifted} of {int((pb < 0.5).sum())} | {dropped} of {int((pb >= 0.9).sum())} | {100*share:.1f} % | {rb:+.4f} / {ra:+.4f} |")
out["attribution"] = attr

(WORK / "analysis262.json").write_text(json.dumps(out, indent=1, default=float))
(WORK / "tables262.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
