"""267 -- analysis: the oracle table with paired target bootstraps, per fold and pooled, by 262 group and by
reliability tercile; the reliability distribution; item 5's sanity rows.  ARITHMETIC ON SAVED ARRAYS ONLY.

Inputs: work/sweep267-<fold>.json, work/pertarget267-<fold>.npz, work/half/truth-<fold>.npz, work/inv/<fold>-targets.csv,
        experiments/263-.../work/cov263-<fold>.npz (262's base.truth cosine -> groups), 263's tables for the reference.
Outputs: work/analysis267.json, work/tables267.md, work/analyse267.log.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/267-replicate-validated-oracle"
WORK = HERE / "work"
E263 = ROOT / "experiments/263-donor-selection-arbiters/work"
FOLDS = [f for f in ("k562_gwps", "orion_hct116", "orion_hek293t") if (WORK / f"sweep267-{f}.json").exists()]
B = 20_000
SEED = 20260905
rng = np.random.default_rng(SEED)
LINES: list[str] = []
MD: list[str] = []
OUT: dict = dict(seed=SEED, B=B, folds=FOLDS)


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def boot(diff: np.ndarray) -> dict:
    """Paired target bootstrap of the fold mean of a per-target difference (NaN targets excluded)."""
    d = diff[np.isfinite(diff)]
    n = d.size
    idx = rng.integers(0, n, size=(B, n))
    means = d[idx].mean(1)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return dict(point=float(d.mean()), lo=float(lo), hi=float(hi), sd=float(means.std(ddof=1)), p_le0=float((means <= 0).mean()),
                p_ge0=float((means >= 0).mean()), n=int(n), n_nonzero=int((d != 0).sum()))


def boot_pool(diffs: list[np.ndarray]) -> dict:
    """The K-fold mean of fold means, each fold resampled within itself."""
    acc = np.zeros(B)
    for diff in diffs:
        d = diff[np.isfinite(diff)]
        idx = rng.integers(0, d.size, size=(B, d.size))
        acc += d[idx].mean(1)
    means = acc / len(diffs)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return dict(point=float(np.mean([np.nanmean(d) for d in diffs])), lo=float(lo), hi=float(hi), sd=float(means.std(ddof=1)),
                p_le0=float((means <= 0).mean()), p_ge0=float((means >= 0).mean()), per_fold=[float(np.nanmean(d)) for d in diffs])


def boot_ratio(num: list[np.ndarray], den: list[np.ndarray]) -> dict:
    """Retained fraction = (pooled mean of num) / (pooled mean of den) under the SAME paired resamples."""
    an = np.zeros(B)
    ad = np.zeros(B)
    for x, y in zip(num, den):
        ok = np.isfinite(x) & np.isfinite(y)
        x, y = x[ok], y[ok]
        idx = rng.integers(0, x.size, size=(B, x.size))
        an += x[idx].mean(1)
        ad += y[idx].mean(1)
    r = an / ad
    lo, hi = np.percentile(r, [2.5, 97.5])
    pt = np.mean([np.nanmean(x) for x in num]) / np.mean([np.nanmean(y) for y in den])
    return dict(point=float(pt), lo=float(lo), hi=float(hi))


def fb(b: dict) -> str:
    return f"{b['point']:+.4f} [{b['lo']:+.4f}, {b['hi']:+.4f}]"


# ---------------------------------------------------------------------------------------------- load
D = {}
for f in FOLDS:
    sw = json.loads((WORK / f"sweep267-{f}.json").read_text())
    rows = {r["arm"]: r for r in sw["rows"]}
    per = {k: np.asarray(v, dtype=np.float64) for k, v in np.load(WORK / f"pertarget267-{f}.npz").items()}
    tr = np.load(WORK / f"half/truth-{f}.npz", allow_pickle=True)
    cov = np.load(E263 / f"cov263-{f}.npz", allow_pickle=True)
    inv = pd.read_csv(WORK / f"inv/{f}-targets.csv").set_index("target")
    targets = [str(t) for t in sw["meta"]["targets"]]
    assert [str(t) for t in tr["targets"]] == targets and [str(t) for t in cov["targets"]] == targets
    base_cos = np.asarray(cov["base262_cos_truth"], dtype=np.float64)
    group = np.where(base_cos <= 0, "e_anti", np.where(base_cos <= 0.05, "marginal", "saturated"))
    cells = inv.loc[targets, "cells"].to_numpy()
    D[f] = dict(meta=sw["meta"], rows=rows, per=per, tr=tr, group=group, base_cos=base_cos, cells=cells, targets=targets,
                guide_targets=[str(t) for t in sw["meta"].get("guide_targets", [])])
    P(f"[{f}] {len(targets)} targets; KA {sw['meta']['ka_pass']}/{sw['meta']['ka_total']} failures {sw['meta']['failures']}; "
      f"groups e/marg/sat {int((group == 'e_anti').sum())}/{int((group == 'marginal').sum())}/{int((group == 'saturated').sum())}")
    OUT.setdefault("ka", {})[f] = sw["meta"]["ka"]

# ---------------------------------------------------------------------------------------------- item 3: reliability
P("\n== item 3: reliability of the reference, cos(truth_A, truth_B) in the instrument's residual space ==")
MD.append("## R. Reliability of the held reference: cos(truth_A, truth_B) per target\n")
MD.append("| fold | split | n | cells/half (median) | mean | p10 | p25 | median | p75 | p90 | share ≤ 0 | share ≥ 0.05 | full-space median | panel-centred median | ρ(cells) | (e) | marginal | saturated |")
MD.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
OUT["reliability"] = {}
for f in FOLDS:
    d = D[f]
    for split in ("batch", "rand", "guide"):
        k = f"rel|{split}|cos_resid"
        if k not in d["tr"].files:
            continue
        c = np.asarray(d["tr"][k], dtype=np.float64)
        ok = np.isfinite(c)
        full = np.asarray(d["tr"][f"rel|{split}|cos_full"], dtype=np.float64)
        pan = np.asarray(d["tr"][f"rel|{split}|cos_panel"], dtype=np.float64)
        ncA = np.asarray(d["tr"][f"rel|{split}|ncells_A"], dtype=np.float64)
        ncB = np.asarray(d["tr"][f"rel|{split}|ncells_B"], dtype=np.float64)
        half_cells = np.minimum(ncA, ncB)
        from scipy.stats import spearmanr
        rho = float(spearmanr(c[ok], d["cells"][ok]).statistic) if ok.sum() > 3 else float("nan")
        g = d["group"]
        by_g = {gg: float(np.nanmean(c[ok & (g == gg)])) if (ok & (g == gg)).any() else float("nan") for gg in ("e_anti", "marginal", "saturated")}
        q = np.percentile(c[ok], [10, 25, 50, 75, 90])
        rec = dict(n=int(ok.sum()), median_cells_half=float(np.median(half_cells[ok])), mean=float(c[ok].mean()), p10=float(q[0]), p25=float(q[1]),
                   median=float(q[2]), p75=float(q[3]), p90=float(q[4]), share_le0=float((c[ok] <= 0).mean()), share_ge005=float((c[ok] >= 0.05).mean()),
                   full_median=float(np.nanmedian(full[ok])), panel_median=float(np.nanmedian(pan[ok])), rho_cells=rho, by_group=by_g)
        OUT["reliability"][f"{f}|{split}"] = rec
        MD.append(f"| `{f}` | {split} | {rec['n']} | {rec['median_cells_half']:.0f} | {rec['mean']:+.4f} | {rec['p10']:+.4f} | {rec['p25']:+.4f} | "
                  f"{rec['median']:+.4f} | {rec['p75']:+.4f} | {rec['p90']:+.4f} | {rec['share_le0']:.3f} | {rec['share_ge005']:.3f} | "
                  f"{rec['full_median']:+.4f} | {rec['panel_median']:+.4f} | {rho:+.2f} | {by_g['e_anti']:+.4f} | {by_g['marginal']:+.4f} | {by_g['saturated']:+.4f} |")
        P(f"  {f} {split}: n {rec['n']} median {rec['median']:+.4f} (p10 {rec['p10']:+.4f}, p90 {rec['p90']:+.4f}) share<=0 {rec['share_le0']:.3f} "
          f"rho(cells) {rho:+.2f}  by group e/marg/sat {by_g['e_anti']:+.4f}/{by_g['marginal']:+.4f}/{by_g['saturated']:+.4f}")
    # reliability vs 262's per-target pds (cd4svd5) and base.truth
    c = np.asarray(d["tr"]["rel|batch|cos_resid"], dtype=np.float64)
    pds262 = np.asarray(np.load(E263 / f"cov263-{f}.npz")["cd4svd5_262_pds"], dtype=np.float64)
    from scipy.stats import spearmanr
    OUT["reliability"][f"{f}|batch"]["rho_pds262"] = float(spearmanr(c, pds262, nan_policy="omit").statistic)
    OUT["reliability"][f"{f}|batch"]["rho_basecos"] = float(spearmanr(c, d["base_cos"], nan_policy="omit").statistic)
    P(f"  {f}: Spearman(reliability, cd4svd5 per-target pds) {OUT['reliability'][f + '|batch']['rho_pds262']:+.3f}; "
      f"(reliability, base.truth cosine) {OUT['reliability'][f + '|batch']['rho_basecos']:+.3f}")

# ---------------------------------------------------------------------------------------------- item 4: the oracle table
P("\n== item 4: the replicate-validated oracle (norm-matched; a23 pds; Δ vs the equal-weight mean scored on the SAME half) ==")
OUT["oracle"] = {}


def get(f, arm, truth):
    return D[f]["per"][f"{arm}|a23|{truth}|pds"]


def arm_table(split: str, gate: str, only=None):
    """Rows: the arms; for each direction (A->B, B->A) the Δ vs the mean on the scoring half; pooled over folds."""
    res = {}
    A, Bn = f"{split}A", f"{split}B"
    if split == "guide":
        A, Bn = "guidemixA", "guidemixB"
    spec = [
        ("(i) mean (cd4svd5's base)", lambda f, x, y: get(f, "mean", y), None),
        ("(ii) in-sample oracle (chosen on the scoring half)", lambda f, x, y: get(f, f"oracle:{y}:{gate}", y), None),
        ("(iii) CROSS-VALIDATED oracle (chosen on the other half)", lambda f, x, y: get(f, f"oracle:{x}:{gate}", y), None),
        ("(iii') oracle chosen on the FULL (deep) truth, scored on the half", lambda f, x, y: get(f, f"oracle:cache:{gate}", y), None),
        ("(iv) random pick, 3 seeds", lambda f, x, y: np.mean([get(f, f"random:{gate}:s{k}", y) for k in range(3)], axis=0), None),
        ("(iii)−(iv) cross-validated oracle MINUS random pick (information above a coin flip)",
         lambda f, x, y: get(f, f"oracle:{x}:{gate}", y) - np.mean([get(f, f"random:{gate}:s{k}", y) for k in range(3)], axis=0) + get(f, "mean", y), None),
        ("(v) 263 arbiter union:wht:tau0:nm", lambda f, x, y: get(f, "arb:union:wht:tau0:nm", y), None),
        ("(v) 263 arbiter basal:raw:tau0:nm", lambda f, x, y: get(f, "arb:basal:raw:tau0:nm", y), None),
        ("(v) 263 arbiter cd4d:wht:tau0:nm", lambda f, x, y: get(f, "arb:cd4d:wht:tau0:nm", y), None),
    ]
    if gate == "all":
        spec += [
            ("(ii-3) in-sample ABSTAINING oracle (best of donor 1 / donor 2 / mean)", lambda f, x, y: get(f, f"oracle3:{y}:all", y), None),
            ("(iii-3) cross-validated abstaining oracle", lambda f, x, y: get(f, f"oracle3:{x}:all", y), None),
        ]
    for label, fn, _ in spec:
        per_dir = {}
        for x, y in ((A, Bn), (Bn, A)):
            diffs = []
            absol = []
            for f in FOLDS:
                v = fn(f, x, y)
                m = get(f, "mean", y)
                mask = np.ones(len(v), bool) if only is None else only(f)
                dd = np.where(mask, v - m, np.nan)
                diffs.append(dd)
                absol.append(float(np.nanmean(np.where(mask, v, np.nan))))
            per_dir[f"{x}->{y}"] = dict(pooled=boot_pool(diffs), per_fold=[boot(dd) for dd in diffs], abs=absol)
        # both directions averaged (each target contributes its A->B and B->A Δ)
        diffs_sym = []
        for f in FOLDS:
            mask = np.ones(len(get(f, "mean", A)), bool) if only is None else only(f)
            d1 = fn(f, A, Bn) - get(f, "mean", Bn)
            d2 = fn(f, Bn, A) - get(f, "mean", A)
            diffs_sym.append(np.where(mask, 0.5 * (d1 + d2), np.nan))
        per_dir["sym"] = dict(pooled=boot_pool(diffs_sym), per_fold=[boot(dd) for dd in diffs_sym])
        res[label] = per_dir
    return res, (A, Bn)


def emit_table(title: str, res: dict, AB, retained_key=("(ii) in-sample oracle (chosen on the scoring half)", "(iii) CROSS-VALIDATED oracle (chosen on the other half)")):
    A, Bn = AB
    MD.append(f"\n### {title}\n")
    MD.append(f"| arm | Δ`pds` {A}→{Bn} pooled [95 % CI] | per fold | Δ`pds` {Bn}→{A} pooled | per fold | symmetric pooled | P(≤ 0) sym |")
    MD.append("|---|---|---|---|---|---|---:|")
    for label, pd_ in res.items():
        a = pd_[f"{A}->{Bn}"]
        b = pd_[f"{Bn}->{A}"]
        s = pd_["sym"]
        MD.append(f"| {label} | {fb(a['pooled'])} | {' / '.join(f'{x:+.4f}' for x in a['pooled']['per_fold'])} | {fb(b['pooled'])} | "
                  f"{' / '.join(f'{x:+.4f}' for x in b['pooled']['per_fold'])} | {fb(s['pooled'])} | {s['pooled']['p_le0']:.4f} |")
        P(f"  {label:70s} {A}->{Bn} {fb(a['pooled'])}  {Bn}->{A} {fb(b['pooled'])}  sym {fb(s['pooled'])} P<=0 {s['pooled']['p_le0']:.4f}")


for split in ("batch", "rand"):
    for gate in ("all", "tau0"):
        res, AB = arm_table(split, gate)
        OUT["oracle"][f"{split}|{gate}"] = res
        emit_table(f"Oracle table — {split} split, gate {gate} (norm-matched; Δ`pds` vs the mean scored on the same half; paired target bootstrap, {B:,} draws, stratified by fold)", res, AB)
        # retained fraction
        A, Bn = AB
        ins, crs = [], []
        for f in FOLDS:
            m1, m2 = get(f, "mean", Bn), get(f, "mean", A)
            ins.append(0.5 * ((get(f, f"oracle:{Bn}:{gate}", Bn) - m1) + (get(f, f"oracle:{A}:{gate}", A) - m2)))
            crs.append(0.5 * ((get(f, f"oracle:{A}:{gate}", Bn) - m1) + (get(f, f"oracle:{Bn}:{gate}", A) - m2)))
        rf = boot_ratio(crs, ins)
        OUT["oracle"][f"{split}|{gate}"]["retained_fraction"] = rf
        MD.append(f"\n**Retained fraction** (symmetric cross-validated Δ / symmetric in-sample Δ, same paired resamples): **{rf['point']:+.3f}** [{rf['lo']:+.3f}, {rf['hi']:+.3f}]\n")
        P(f"  RETAINED FRACTION {split}/{gate}: {rf['point']:+.3f} [{rf['lo']:+.3f}, {rf['hi']:+.3f}]")

# ---------------------------------------------------------------------------------------------- by 262 group and by reliability tercile
P("\n== by 262 group and by reliability tercile (batch split, gate all; symmetric Δ per target; sums = contribution to the fold Δ) ==")
MD.append("\n### By 262 group and by reliability tercile (batch split, gate all; per-target symmetric Δ`pds` vs the mean, i.e. the mean of A→B and B→A; pooled over the three folds)\n")
MD.append("| stratum | n | in-sample oracle: mean Δ per target | share of the in-sample gain | cross-validated: mean Δ per target [95 % CI] | share | retained | random pick | cross − random [95 % CI] | 263 union arbiter |")
MD.append("|---|---:|---:|---:|---|---:|---:|---:|---|---:|")
OUT["strata"] = {}
for split in ("batch", "rand"):
    A, Bn = f"{split}A", f"{split}B"
    per_t = {}
    for f in FOLDS:
        m1, m2 = get(f, "mean", Bn), get(f, "mean", A)
        per_t[f] = dict(
            ins=0.5 * ((get(f, f"oracle:{Bn}:all", Bn) - m1) + (get(f, f"oracle:{A}:all", A) - m2)),
            crs=0.5 * ((get(f, f"oracle:{A}:all", Bn) - m1) + (get(f, f"oracle:{Bn}:all", A) - m2)),
            rnd=0.5 * ((np.mean([get(f, f"random:all:s{k}", Bn) for k in range(3)], axis=0) - m1) + (np.mean([get(f, f"random:all:s{k}", A) for k in range(3)], axis=0) - m2)),
            arb=0.5 * ((get(f, "arb:union:wht:tau0:nm", Bn) - m1) + (get(f, "arb:union:wht:tau0:nm", A) - m2)),
            rel=np.asarray(D[f]["tr"][f"rel|{split}|cos_resid"], dtype=np.float64), group=D[f]["group"])
        rel = per_t[f]["rel"]
        q = np.nanpercentile(rel, [100 / 3, 200 / 3])
        per_t[f]["terc"] = np.where(rel <= q[0], "T1 low", np.where(rel <= q[1], "T2 mid", "T3 high"))
        per_t[f]["terc_q"] = q.tolist()
    strata = [("group (e) base·truth ≤ 0", "group", "e_anti"), ("marginal band (0, 0.05]", "group", "marginal"), ("saturated > 0.05", "group", "saturated"),
              ("reliability tercile T1 (low)", "terc", "T1 low"), ("reliability tercile T2 (mid)", "terc", "T2 mid"), ("reliability tercile T3 (high)", "terc", "T3 high")]
    tot_ins = sum(float(np.nansum(per_t[f]["ins"])) for f in FOLDS)
    for label, key, val in strata:
        masks = [per_t[f][key] == val for f in FOLDS]
        n = sum(int(m.sum()) for m in masks)
        ins = [np.where(m, per_t[f]["ins"], np.nan) for f, m in zip(FOLDS, masks)]
        crs = [np.where(m, per_t[f]["crs"], np.nan) for f, m in zip(FOLDS, masks)]
        rnd = [np.where(m, per_t[f]["rnd"], np.nan) for f, m in zip(FOLDS, masks)]
        arb = [np.where(m, per_t[f]["arb"], np.nan) for f, m in zip(FOLDS, masks)]
        s_ins = sum(float(np.nansum(x)) for x in ins)
        s_crs = sum(float(np.nansum(x)) for x in crs)
        b_crs = boot_pool(crs)
        b_cmr = boot_pool([x - y for x, y in zip(crs, rnd)])          # cross-validated oracle minus the random pick
        mean_ins = float(np.nanmean(np.concatenate(ins)))
        rf = s_crs / s_ins if s_ins != 0 else float("nan")
        rec = dict(n=n, ins_mean=mean_ins, ins_sum=s_ins, ins_share=s_ins / tot_ins if tot_ins else float("nan"), crs=b_crs, crs_sum=s_crs, retained=rf,
                   rnd_mean=float(np.nanmean(np.concatenate(rnd))), arb_mean=float(np.nanmean(np.concatenate(arb))), cross_minus_random=b_cmr)
        OUT["strata"][f"{split}|{label}"] = rec
        if split == "batch":
            MD.append(f"| {label} | {n} | {mean_ins:+.4f} | {rec['ins_share'] * 100:.0f} % | {fb(b_crs)} | {(s_crs / tot_ins) * 100 if tot_ins else float('nan'):.0f} % | "
                      f"{rf:+.2f} | {rec['rnd_mean']:+.4f} | {fb(b_cmr)} | {rec['arb_mean']:+.4f} |")
        P(f"  {split} {label:32s} n {n:4d}  in-sample {mean_ins:+.4f}/target ({rec['ins_share'] * 100:.0f} % of gain)  cross {fb(b_crs)}  retained {rf:+.2f}  random {rec['rnd_mean']:+.4f}  cross-random {fb(b_cmr)}")
    OUT["strata"][f"{split}|tercile_cuts"] = {f: per_t[f]["terc_q"] for f in FOLDS}

# ---------------------------------------------------------------------------------------------- pick agreement and the accuracy framing
P("\n== pick agreement (the oracle's choice on A vs on B; 263's break-even is 57-59 % at tau0, 66-70 % at all) ==")
MD.append("\n### Pick agreement between the halves' oracles (263's break-even accuracy: 57–59 % at τ0, 66–70 % at gate all)\n")
MD.append("| fold | batchA~batchB all | tau0 | randA~randB all | tau0 | batchA~full all | batchB~full all |")
MD.append("|---|---:|---:|---:|---:|---:|---:|")
OUT["agree"] = {}
for f in FOLDS:
    ag = D[f]["meta"]["agree"]
    OUT["agree"][f] = ag
    g = lambda k: f"{ag[k]['agree'] * 100:.0f} % (n {ag[k]['n']})" if k in ag else "—"
    MD.append(f"| `{f}` | {g('batchA~batchB:all')} | {g('batchA~batchB:tau0')} | {g('randA~randB:all')} | {g('randA~randB:tau0')} | {g('batchA~cache:all')} | {g('batchB~cache:all')} |")
    P(f"  {f}: " + ", ".join(f"{k} {v['agree'] * 100:.0f}% (n {v['n']})" for k, v in ag.items()))

# ---------------------------------------------------------------------------------------------- item 5: sanity on the split
P("\n== item 5: sanity ==")
MD.append("\n### Item 5 — sanity on the split\n")
MD.append("| fold | mean `pds` on the full truth | on batchA | on batchB | on randA | on randB | in-sample oracle on full / batchA / batchB |")
MD.append("|---|---:|---:|---:|---:|---:|---|")
OUT["sanity"] = {}
for f in FOLDS:
    r = D[f]["rows"]
    m = r["mean"]
    rec = dict(mean=dict(cache=m["pds|cache"], batchA=m["pds|batchA"], batchB=m["pds|batchB"], randA=m["pds|randA"], randB=m["pds|randB"]),
               oracle_in=dict(cache=r["oracle:cache:all"]["pds|cache"], batchA=r["oracle:batchA:all"]["pds|batchA"], batchB=r["oracle:batchB:all"]["pds|batchB"]))
    dAB = get(f, "mean", "batchA") - get(f, "mean", "batchB")
    rec["mean_A_minus_B"] = boot(dAB)
    OUT["sanity"][f] = rec
    MD.append(f"| `{f}` | {m['pds|cache']:.4f} | {m['pds|batchA']:.4f} | {m['pds|batchB']:.4f} | {m['pds|randA']:.4f} | {m['pds|randB']:.4f} | "
              f"{rec['oracle_in']['cache']:.4f} / {rec['oracle_in']['batchA']:.4f} / {rec['oracle_in']['batchB']:.4f} |")
    P(f"  {f}: mean pds full {m['pds|cache']:.4f}  batchA {m['pds|batchA']:.4f}  batchB {m['pds|batchB']:.4f} (A−B {fb(rec['mean_A_minus_B'])})  "
      f"randA {m['pds|randA']:.4f} randB {m['pds|randB']:.4f}")

# guide split vs random split on the two-construct targets
MD.append("\n**Guide split vs random-cell split on the two-construct targets** (cross-validated symmetric Δ`pds` vs the mean on those targets only; the guide-mix truths carry the guide halves on these targets and the batch halves elsewhere):\n")
MD.append("| fold | n guide targets | cross-validated Δ, guide split | cross-validated Δ, batch split (same targets) | cross-validated Δ, rand split (same targets) | in-sample Δ guide / batch / rand |")
MD.append("|---|---:|---|---|---|---|")
OUT["guide"] = {}
for f in FOLDS:
    gt = D[f]["guide_targets"]
    if not gt or "oracle:guidemixA:all" not in D[f]["rows"]:
        continue
    mask = np.isin(D[f]["targets"], gt)
    rec = {}
    for split, A, Bn in (("guide", "guidemixA", "guidemixB"), ("batch", "batchA", "batchB"), ("rand", "randA", "randB")):
        m1, m2 = get(f, "mean", Bn), get(f, "mean", A)
        crs = 0.5 * ((get(f, f"oracle:{A}:all", Bn) - m1) + (get(f, f"oracle:{Bn}:all", A) - m2))
        ins = 0.5 * ((get(f, f"oracle:{Bn}:all", Bn) - m1) + (get(f, f"oracle:{A}:all", A) - m2))
        rec[split] = dict(cross=boot(np.where(mask, crs, np.nan)), ins=float(np.nanmean(np.where(mask, ins, np.nan))))
    OUT["guide"][f] = rec
    MD.append(f"| `{f}` | {int(mask.sum())} | {fb(rec['guide']['cross'])} | {fb(rec['batch']['cross'])} | {fb(rec['rand']['cross'])} | "
              f"{rec['guide']['ins']:+.4f} / {rec['batch']['ins']:+.4f} / {rec['rand']['ins']:+.4f} |")
    P(f"  {f} guide targets n {int(mask.sum())}: cross Δ guide {fb(rec['guide']['cross'])}, batch {fb(rec['batch']['cross'])}, rand {fb(rec['rand']['cross'])}; "
      f"in-sample guide/batch/rand {rec['guide']['ins']:+.4f}/{rec['batch']['ins']:+.4f}/{rec['rand']['ins']:+.4f}")

# ---------------------------------------------------------------------------------------------- the in-sample oracle on the full truth vs 263 (reference)
MD.append("\n### Reference: 263's oracle on the full (deep) truth, reproduced here (KA), and the same oracle scored on the halves\n")
MD.append("| fold | Δ`pds` oracle_wht:all:nm on the full truth (263) | here, chosen on full, scored on full | chosen on full, scored on batchA / batchB | chosen on batchA, scored on full (contains A) |")
MD.append("|---|---:|---:|---|---:|")
for f in FOLDS:
    r = D[f]["rows"]
    j263 = json.loads((E263 / f"sweep263-{f}.json").read_text())
    r263 = {x["arm"]: x for x in j263["rows"]}
    d263 = r263["oracle_wht:all:nm"]["pds"] - r263["cd4svd5"]["pds"]
    here = r["oracle:cache:all"]["pds|cache"] - r["mean"]["pds|cache"]
    hA = r["oracle:cache:all"]["pds|batchA"] - r["mean"]["pds|batchA"]
    hB = r["oracle:cache:all"]["pds|batchB"] - r["mean"]["pds|batchB"]
    cA = r["oracle:batchA:all"]["pds|cache"] - r["mean"]["pds|cache"]
    MD.append(f"| `{f}` | {d263:+.4f} | {here:+.4f} | {hA:+.4f} / {hB:+.4f} | {cA:+.4f} |")
    P(f"  {f}: 263 oracle Δ {d263:+.4f}; here {here:+.4f}; chosen on full scored on halves {hA:+.4f}/{hB:+.4f}; chosen on batchA scored on full {cA:+.4f}")

(WORK / "analysis267.json").write_text(json.dumps(OUT, indent=1, default=float))
(WORK / "tables267.md").write_text("\n".join(MD) + "\n")
(WORK / "analyse267.log").write_text("\n".join(LINES) + "\n")
P("\nwritten work/analysis267.json, work/tables267.md, work/analyse267.log")
