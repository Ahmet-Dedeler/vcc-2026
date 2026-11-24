"""272 -- analysis on saved arrays: the spectra (A on three folds, C on production), the anti-correlation tables (A's
KA-0 and B's cross-fit test), every arm vs cd4svd5 with paired target bootstraps (per fold and pooled, 20,000 draws,
seed 20260906), by 262 group and reference-reliability tercile, the marginal band, leave-top-3-out retention, and the
arm-minus-control paired intervals.  ARITHMETIC ON SAVED ARRAYS ONLY.

Inputs: work/sweep272-<fold>.{json}, pertarget272-<fold>.npz, identity272-<fold>.npz, spectrum272-<fold>.npz,
        work/prod272.json, work/refit272-orion_hct116.json + pertarget272-refit-*.npz (if present).
Outputs: work/analysis272.json, work/tables272.md, work/analyse272.log.
    PYTHONDONTWRITEBYTECODE=1 uv run python analyse272.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/272-noise-corrected-operator"
WORK = HERE / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t")
LINE = {"k562_gwps": "K562-GWPS", "orion_hct116": "HCT116", "orion_hek293t": "HEK293T"}
B = 20_000
SEED = 20260906
rng = np.random.default_rng(SEED)
ARMS_A = [f"wiener:k{k:g}" for k in (0.25, 0.5, 1, 2, 4)] + [f"floor:c{c:g}" for c in (0.25, 0.5, 1, 2, 4)] + ["hard:0.5"]
CONTROLS_A = [f"wiener:k{k:g}:{c}" for k in (0.5, 1, 2) for c in ("shuf", "mirror")] + ["floor:c1:shuf", "floor:c1:mirror", "hard:0.5:mirror"]
CONTROL_OF = {**{f"wiener:k{k:g}": [f"wiener:k{k:g}:shuf", f"wiener:k{k:g}:mirror"] for k in (0.5, 1, 2)},
              "wiener:k0.25": ["wiener:k0.5:shuf", "wiener:k0.5:mirror"], "wiener:k4": ["wiener:k2:shuf", "wiener:k2:mirror"],
              **{f"floor:c{c:g}": ["floor:c1:shuf", "floor:c1:mirror"] for c in (0.25, 0.5, 1, 2, 4)}, "hard:0.5": ["hard:0.5:mirror"]}
ARMS_B = ["refit_full", "op_A", "op_B", "xfit_avg", "half:inc", "half:full", "half:avgop", "half:xfit", "half:infit"]
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
                p_ge0=float((means >= 0).mean()), n=int(n))


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


# ------------------------------------------------------------------------------------------------- load
D = {}
for f in FOLDS:
    sw = json.loads((WORK / f"sweep272-{f}.json").read_text())
    rows = {r["arm"]: r for r in sw["rows"]}
    per = {k: np.asarray(v) for k, v in np.load(WORK / f"pertarget272-{f}.npz", allow_pickle=True).items()}
    idn = {k: np.asarray(v, dtype=np.float64) for k, v in np.load(WORK / f"identity272-{f}.npz").items()}
    base_cos = np.asarray(per["base262_cos_truth"], dtype=np.float64)
    group = np.where(base_cos <= 0, "e_anti", np.where(base_cos <= 0.05, "marginal", "saturated"))
    ref = np.asarray(per["ref_rel"], dtype=np.float64)
    cuts = np.nanpercentile(ref, [100 / 3, 200 / 3])
    terc = np.where(ref <= cuts[0], "T1 low", np.where(ref <= cuts[1], "T2 mid", "T3 high"))
    D[f] = dict(meta=sw["meta"], rows=rows, per=per, idn=idn, targets=[str(t) for t in sw["meta"]["targets"]], group=group, terc=terc)
    m = sw["meta"]
    P(f"[{f}] {m['n_targets']} targets, donors {m['donors']}, KA {m['ka_pass']}/{m['ka_total']} failures {m['failures']}, {m['n_rows']} rows, "
      f"{m['seconds']}s, peak {m['peak_rss_gib']:.2f} GiB; op272 {m['op272_sha256'][:12]}")
    OUT.setdefault("ka", {})[f] = m["ka"]
    OUT.setdefault("op272_sha256", m["op272_sha256"])
    assert m["op272_sha256"] == OUT["op272_sha256"], "op272.py changed between folds"
    assert not m["failures"], f"{f}: KA failures {m['failures']}"
prod = json.loads((WORK / "prod272.json").read_text())
assert prod["op272_sha256"] == OUT["op272_sha256"]
REFITS = {}
for rf in FOLDS:
    refit_path = WORK / f"refit272-{rf}.json"
    if not refit_path.exists():
        continue
    rj = json.loads(refit_path.read_text())
    REFITS[rf] = dict(meta=rj["meta"], rows={r["arm"]: r for r in rj["rows"]},
                      per={k: np.asarray(v) for k, v in np.load(WORK / f"pertarget272-refit-{rf}.npz", allow_pickle=True).items()},
                      idn={k: np.asarray(v, dtype=np.float64) for k, v in np.load(WORK / f"identity272-refit-{rf}.npz").items()})
    assert REFITS[rf]["meta"]["op272_sha256"] == OUT["op272_sha256"]
    m = REFITS[rf]["meta"]
    P(f"[refit {rf}] KA {m['ka_pass']}/{m['ka_total']} failures {m['failures']}; fit KA {m['fit']['ka_pass']}/{m['fit']['ka_total']}")
R = REFITS.get("orion_hct116")

# ------------------------------------------------------------------------------------------------- the spectra
P("\n== the s_i spectra (A: three fold operators; C: production) ==")
MD.append("## S1. The noise spectrum in the fitted eigenbasis: s_i = clip((λ_i − n_i)/λ_i, 0, 1) by eigenvalue-rank decile (decile 0 = the largest eigenvalues); "
          "`amp` = median λ^(−γ/2) relative to the top eigenvalue; `n/λ` = median noise-over-eigenvalue; energy shares of the prediction rows before / after whitening\n")
MD.append("| operator | decile | λ range | amp | mean s | share s < 0.5 | share s = 0 | n/λ median | raw energy | whitened energy |")
MD.append("|---|---:|---|---:|---:|---:|---:|---:|---:|---:|")
OUT["spectra"] = {}
specs = {f: D[f]["meta"]["spectrum"] for f in FOLDS}
specs["production"] = prod["spectrum"]
for name, rep in specs.items():
    OUT["spectra"][name] = {k: v for k, v in rep.items() if k != "deciles"}
    OUT["spectra"][name]["deciles"] = rep["deciles"]
    for r in rep["deciles"]:
        MD.append(f"| `{name}` | {r['decile']} | {r['lambda_min']:.1e} – {r['lambda_max']:.1e} | {r['amplification_median']:.0f}× | {r['mean_s']:.3f} | {r['share_s_lt_half']:.2f} | "
                  f"{r['share_s_zero']:.2f} | {r['median_noise_over_lambda']:.2f} | {100*r['raw_energy_share']:.1f} % | {100*r['whitened_energy_share']:.1f} % |")
MD.append("\n## S2. The noise-dominated set per operator\n")
MD.append("| operator (basis rows) | p | floored λ | s < 0.5 | s = 0 | rank of the s < 0.5 set: first / p10 / median / p90 | whitened energy in s < 0.5 | in s = 0 | raw energy in s < 0.5 | Σn / Σλ | rank corr(log n, log λ) | half-factor ratio per donor | signal share of the donor rows |")
MD.append("|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---|---|")
for name, rep in specs.items():
    bw = rep["basis_weights"]
    hf = rep["half_factor_check"]
    MD.append(f"| `{name}` ({', '.join(f'{k} {int(v)}' for k, v in bw.items())}) | {rep['p']} | {rep['n_lambda_floored']} | {rep['n_noise_dominated']} ({100*rep['share_noise_dominated']:.1f} %) | "
              f"{rep['n_s_zero']} ({100*rep['share_s_zero']:.1f} %) | {rep['rank_of_noise_dominated']['first']} / {rep['rank_of_noise_dominated']['p10']:.0f} / {rep['rank_of_noise_dominated']['median']:.0f} / {rep['rank_of_noise_dominated']['p90']:.0f} | "
              f"**{100*rep['whitened_energy_share_noise_dominated']:.1f} %** | {100*rep['whitened_energy_share_s_zero']:.1f} % | {100*rep['raw_energy_share_noise_dominated']:.1f} % | {rep['total_noise_over_total_lambda']:.3f} | "
              f"{rep.get('spearman_log_n_log_lambda', float('nan')):+.3f} | {' / '.join(f'{k}: {v['ratio_d_over_implied']:.3f}' for k, v in hf.items())} | {' / '.join(f'{k}: {v['signal_share_of_real_energy']:.3f}' for k, v in hf.items())} |")
    P(f"  {name:14s} s<0.5 {100*rep['share_noise_dominated']:.1f} % (s=0 {100*rep['share_s_zero']:.1f} %) carrying {100*rep['whitened_energy_share_noise_dominated']:.1f} % of the whitened energy "
      f"(s=0: {100*rep['whitened_energy_share_s_zero']:.1f} %); rank corr {rep.get('spearman_log_n_log_lambda', float('nan')):+.3f}; sum n / sum lambda {rep['total_noise_over_total_lambda']:.3f}")

# ------------------------------------------------------------------------------------------------- the anti-correlation tables
P("\n== anti-correlation (KA-0 on the fold operators; the production operator; B's cross-fit test) ==")
MD.append("\n## AC1. Whitened split-half cosine cos(WA_t, WB_t), median over the carried panel rows (271's statistic, KA-0 at 0.0e+00 against `diagwht271.json`)\n")
MD.append("| operator | line | role | n | raw median | whitened median | share ≤ 0 | null median (sd) |")
MD.append("|---|---|---|---:|---:|---:|---:|---|")
OUT["anti"] = {}
diag271 = json.loads((ROOT / "experiments/271-donor-reliability-estimation/work/diagwht271.json").read_text())
for f in FOLDS:
    for line, a in D[f]["meta"]["anti_correlation"].items():
        role = a.get("role", "donor (in-sample)")
        OUT["anti"][f"{f}|{line}"] = a
        raw = a.get("raw_median", diag271[f][line]["raw"]["median"])          # the held line's raw median is 271's (not recomputed in the sweep)
        null = f"{a['wht_null_median']:+.4f} ({a['wht_null_sd']:.3f})" if "wht_null_median" in a else f"{diag271[f][line]['wht']['null_median']:+.4f} ({diag271[f][line]['wht']['null_sd']:.3f})"
        MD.append(f"| `held_{f}` | `{line}` | {role} | {a['n']} | {raw:+.4f} | **{a['wht_median']:+.4f}** | {a['wht_share_le0']:.3f} | {null} |")
for line, a in prod["anti_correlation"].items():
    OUT["anti"][f"production|{line}"] = a
    MD.append(f"| `union_pruned` (production) | `{line}` | donor (in-sample for the basis rows drawn) | {a['n']} | {a['raw_median']:+.4f} | **{a['wht_median']:+.4f}** | {a['wht_share_le0']:.3f} | {a['wht_null_median']:+.4f} ({a['wht_null_sd']:.3f}) |")
if REFITS:
    MD.append("\n## AC2. Construction B: the same statistic through the refitted operators — cached (both halves in-sample), the full-row refit, W_A (half A in-sample, half B out), W_B (the mirror)\n")
    MD.append("| fold | line | role | n | raw | cached | full refit | **W_A** (A in / B out) | share ≤ 0 | **W_B** (B in / A out) | share ≤ 0 | null sd (W_A) |")
    MD.append("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for rf, Rf in REFITS.items():
        for line, a in Rf["meta"]["anti_correlation"].items():
            OUT["anti"][f"refit|{rf}|{line}"] = a
            if "role" in a:
                MD.append(f"| `{rf}` | `{line}` | {a['role']} | {a['inc']['n']} | — | {a['inc']['median']:+.4f} | {a['full']['median']:+.4f} | **{a['A']['median']:+.4f}** | {a['A']['share_le0']:.3f} | **{a['B']['median']:+.4f}** | {a['B']['share_le0']:.3f} | — |")
            else:
                MD.append(f"| `{rf}` | `{line}` | donor | {a['inc']['n']} | {a['raw']['median']:+.4f} | {a['inc']['median']:+.4f} | {a['full']['median']:+.4f} | **{a['A']['median']:+.4f}** | {a['A']['share_le0']:.3f} | "
                          f"**{a['B']['median']:+.4f}** | {a['B']['share_le0']:.3f} | {a['A']['null_sd']:.3f} |")
            P(f"  refit {rf} {line:14s}: cached {a['inc']['median']:+.4f} full {a['full']['median']:+.4f} W_A {a['A']['median']:+.4f} W_B {a['B']['median']:+.4f}")


# ------------------------------------------------------------------------------------------------- the measurement (A)
def nmae_concat(Dd, f_list, arm):
    U = np.concatenate([Dd[f]["idn"][f"{arm}|U"] for f in f_list])
    M = np.concatenate([Dd[f]["idn"][f"{arm}|M"] for f in f_list])
    Den = np.concatenate([Dd[f]["idn"][f"{arm}|D"] for f in f_list])
    return float(np.mean(1 + (U - 2 * M) / Den))


P("\n== construction A: every arm vs cd4svd5 (a23 pds on the cache truth; paired target bootstraps, pooled over three folds) ==")
MD.append("\n## M1. Construction A: every arm vs `cd4svd5`, three cancer folds (Δ`pds` per fold K562 / HCT116 / HEK293T, pooled paired target bootstrap; cosine, delivered amplitude, the `nmae` identity)\n")
MD.append("| arm | whitened energy kept | Δ`pds` K562 / HCT116 / HEK293T | **pooled [95 % CI]** | P(≤ 0) | Δ`cos_agg` % | Δ macro % | delivered ratio | Δ`nmae` identity | Δ top-89 precision | base cosine to the incumbent base |")
MD.append("|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|")
OUT["arms"] = {}
nm_base = nmae_concat(D, FOLDS, "cd4svd5")
for arm in ["STACK"] + ARMS_A + CONTROLS_A:
    diffs = [D[f]["per"][f"{arm}|pds"] - D[f]["per"]["cd4svd5|pds"] for f in FOLDS]
    bp = boot_pool(diffs)
    per_fold = {f: boot(dd) for f, dd in zip(FOLDS, diffs)}
    ra = [D[f]["rows"][arm]["cos_agg"] / D[f]["rows"]["cd4svd5"]["cos_agg"] - 1 for f in FOLDS]
    rm = [D[f]["rows"][arm]["a23_cos"] / D[f]["rows"]["cd4svd5"]["a23_cos"] - 1 for f in FOLDS]
    dn = np.mean([D[f]["rows"][arm]["delivered_norm"] for f in FOLDS]) / np.mean([D[f]["rows"]["cd4svd5"]["delivered_norm"] for f in FOLDS])
    dnm = nmae_concat(D, FOLDS, arm) - nm_base
    dprec = float(np.mean([D[f]["rows"][arm]["top89_precision"] - D[f]["rows"]["cd4svd5"]["top89_precision"] for f in FOLDS]))
    kept = [D[f]["rows"][arm]["audit"].get("whitened_energy_kept", 1.0) for f in FOLDS]
    bc = [D[f]["rows"][arm]["audit"]["base_cos_to_incumbent_base_median"] for f in FOLDS]
    rec = dict(pooled=bp, per_fold=per_fold, rel_agg=[float(x) for x in ra], rel_macro=[float(x) for x in rm], delivered_ratio=float(dn), d_nmae=dnm, d_prec89=dprec,
               whitened_energy_kept=[float(x) for x in kept], base_cos=[float(x) for x in bc])
    OUT["arms"][arm] = rec
    MD.append(f"| `{arm}` | {100*np.mean(kept):.0f} % | {' / '.join(f'{x:+.4f}' for x in bp['per_fold'])} | **{fb(bp)}** | {bp['p_le0']:.3f} | {100*np.mean(ra):+.2f} | {100*np.mean(rm):+.2f} | "
              f"{dn:.4f} | {dnm:+.5f} | {dprec:+.4f} | {np.mean(bc):.3f} |")
    P(f"  {arm:20s} kept {100*np.mean(kept):5.1f}%  d pds {' / '.join(f'{x:+.4f}' for x in bp['per_fold'])}  pooled {fb(bp)}  P<=0 {bp['p_le0']:.3f}  agg {100*np.mean(ra):+.2f}%  macro {100*np.mean(rm):+.2f}%  dn {dn:.4f}  d nmae {dnm:+.5f}")

# arm minus its controls
MD.append("\n## M2. Each arm against its controls (paired target bootstrap of the per-target Δ`pds` difference, pooled; positive = the arm loses less than the control)\n")
MD.append("| arm | control | arm − control |")
MD.append("|---|---|---|")
OUT["vs_controls"] = {}
for arm in ARMS_A:
    OUT["vs_controls"][arm] = {}
    for ctl in CONTROL_OF[arm]:
        diffs = [D[f]["per"][f"{arm}|pds"] - D[f]["per"][f"{ctl}|pds"] for f in FOLDS]
        b = boot_pool(diffs)
        OUT["vs_controls"][arm][ctl] = b
        MD.append(f"| `{arm}` | `{ctl}` | {fb(b)} |")

# by 262 group and by reference-reliability tercile
MD.append("\n## M3. Δ`pds` vs `cd4svd5` by 262 group (base·truth cosine of the whitened mean) and by the held reference's split-half reliability tercile (267, batch halves); pooled, mean per target [95 % CI]\n")
strata = [("group (e) base·truth ≤ 0", "group", "e_anti"), ("marginal band (0, 0.05]", "group", "marginal"), ("saturated > 0.05", "group", "saturated"),
          ("ref-reliability T1 (low)", "terc", "T1 low"), ("ref-reliability T2 (mid)", "terc", "T2 mid"), ("ref-reliability T3 (high)", "terc", "T3 high")]
MD.append("| arm | " + " | ".join(f"{s[0]} (n {sum(int((D[f][s[1]] == s[2]).sum()) for f in FOLDS)})" for s in strata) + " |")
MD.append("|---|" + "---|" * len(strata))
OUT["strata"] = {}
for arm in ARMS_A + ["wiener:k1:shuf", "wiener:k1:mirror", "floor:c1:shuf"]:
    rec = {}
    cells = []
    for label, key, val in strata:
        diffs = [np.where(D[f][key] == val, D[f]["per"][f"{arm}|pds"] - D[f]["per"]["cd4svd5|pds"], np.nan) for f in FOLDS]
        b = boot_pool(diffs)
        rec[label] = b
        cells.append(fb(b))
    OUT["strata"][arm] = rec
    MD.append(f"| `{arm}` | " + " | ".join(cells) + " |")

# the marginal band
MD.append("\n## M4. The marginal band (cd4svd5's own cosine in (0, 0.05]): own-cosine change and the count of targets crossing 0.02, pooled\n")
MD.append("| arm | band n | mean Δ own cosine on the band | up across 0.02 | down across 0.02 | net share | mean Δ`pds` on the band [95 % CI] | Δ own cosine, all targets |")
MD.append("|---|---:|---:|---:|---:|---:|---|---:|")
OUT["band"] = {}
for arm in ARMS_A + ["wiener:k1:shuf", "wiener:k1:mirror", "floor:c1:shuf", "floor:c1:mirror"]:
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
    rec = dict(n=nb, up=up, down=down, net_share=(up - down) / nb, d_cos_band=float(np.nanmean(np.concatenate(dc))), d_pds_band=b, d_cos_all=float(np.nanmean(np.concatenate(dc_all))))
    OUT["band"][arm] = rec
    MD.append(f"| `{arm}` | {nb} | {rec['d_cos_band']:+.5f} | {up} | {down} | {100*rec['net_share']:+.1f} % | {fb(b)} | {rec['d_cos_all']:+.5f} |")

# leave-top-3-out retention
MD.append("\n## M5. Leave-top-3-out retention of the per-target Δ`pds` (261's statistic; each fold on its own ordering)\n")
MD.append("| arm | full | leave-top-3-out | retained | per fold |")
MD.append("|---|---:|---:|---:|---|")
OUT["retention"] = {}
for arm in ARMS_A:
    full, loo = [], []
    for f in FOLDS:
        d = D[f]["per"][f"{arm}|pds"] - D[f]["per"]["cd4svd5|pds"]
        keep = np.ones(d.size, bool)
        keep[np.argsort(-d)[:3]] = False
        full.append(d.mean())
        loo.append(d[keep].mean())
    fm, lm = float(np.mean(full)), float(np.mean(loo))
    rec = dict(full=fm, leave_top_out=lm, retained=(lm / fm if fm else float("nan")), per_fold=[float(l / g) if g else float("nan") for l, g in zip(loo, full)])
    OUT["retention"][arm] = rec
    MD.append(f"| `{arm}` | {fm:+.5f} | {lm:+.5f} | {100*rec['retained']:+.0f} % | {' / '.join(f'{100*x:+.0f} %' for x in rec['per_fold'])} |")

# ------------------------------------------------------------------------------------------------- construction B
OUT["refit"] = {}
for rf, R in REFITS.items():
    P(f"\n== construction B ({rf}): the cross-fitted arms vs cd4svd5 and vs the matched half-row baseline ==")
    MD.append(f"\n## B1 ({rf}). Construction B on `{rf}` (donors {' + '.join(R['meta']['donors'])}): the refitted / cross-fitted operators through cd4svd5's construction (paired target bootstrap on this fold)\n")
    MD.append("| arm | vs | Δ`pds` [95 % CI] | P(≤ 0) | Δ`cos_agg` % | Δ macro % | delivered ratio | Δ`nmae` identity | base cosine to the incumbent base |")
    MD.append("|---|---|---|---:|---:|---:|---:|---:|---:|")
    OUT["refit"][rf] = dict(meta={k: v for k, v in R["meta"].items() if k not in ("targets",)}, arms={})
    Rr = {rf: R}

    def nm1(arm, Rr=Rr, rf=rf):
        return nmae_concat(Rr, [rf], arm)

    for arm, base in [("refit_full", "cd4svd5"), ("op_A", "cd4svd5"), ("op_B", "cd4svd5"), ("xfit_avg", "cd4svd5"), ("half:inc", "cd4svd5"),
                      ("half:full", "half:inc"), ("half:avgop", "half:inc"), ("half:xfit", "half:inc"), ("half:infit", "half:inc"), ("half:xfit", "half:infit")]:
        d = R["per"][f"{arm}|pds"] - R["per"][f"{base}|pds"]
        b = boot(d)
        ra = R["rows"][arm]["cos_agg"] / R["rows"][base]["cos_agg"] - 1
        rm = R["rows"][arm]["a23_cos"] / R["rows"][base]["a23_cos"] - 1
        dn = R["rows"][arm]["delivered_norm"] / R["rows"][base]["delivered_norm"]
        dnm = nm1(arm) - nm1(base)
        OUT["refit"][rf]["arms"][f"{arm}|{base}"] = dict(boot=b, rel_agg=float(ra), rel_macro=float(rm), delivered_ratio=float(dn), d_nmae=float(dnm))
        MD.append(f"| `{arm}` | `{base}` | **{fb(b)}** | {b['p_le0']:.3f} | {100*ra:+.2f} | {100*rm:+.2f} | {dn:.4f} | {dnm:+.5f} | {R['rows'][arm]['audit']['base_cos_to_incumbent_base_median']:.3f} |")
        P(f"  {arm:12s} vs {base:10s} d pds {fb(b)} P<=0 {b['p_le0']:.3f} agg {100*ra:+.2f}% macro {100*rm:+.2f}% dn {dn:.4f} d nmae {dnm:+.5f}")

(WORK / "analysis272.json").write_text(json.dumps(OUT, indent=1, default=float))
(WORK / "tables272.md").write_text("\n".join(MD) + "\n")
(WORK / "analyse272.log").write_text("\n".join(LINES) + "\n")
P("\nwritten work/analysis272.json, work/tables272.md, work/analyse272.log")
