"""264 -- the analysis: known-answer rows (R1-R3), the per-fold table for (a) the three-donor mean, (b) the majority
rule, (c) the incumbent and the controls; the two increments (b) - (a) and (b) - (c) with retention and the
three-scenario price (261's `price()` / `scenarios()`, 257's constants); the confound (a) - (c); the verdict.
Reads work/sweep264-<fold>.json and the npz siblings.  Writes work/analysis264.json, work/tables264.md, work/analyse264.log.

    PYTHONDONTWRITEBYTECODE=1 uv run python analyse264.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E261 = ROOT / "experiments/261-union-deflation"
E257 = ROOT / "experiments/257-soft-deflation-price"
sys.path.insert(0, str(E261))
sys.path.insert(0, str(E257))
import price257 as P  # noqa: E402  (constants, raw_for, scaled_terms; read-only at import)
import price261 as Q  # noqa: E402  (price, factors, scenarios, retention; read-only at import)

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
FOLDS = CANCER + ("kolf2",)
LINE = {"k562_gwps": "K562", "orion_hct116": "HCT116", "orion_hek293t": "HEK293T", "kolf2": "KOLF2"}
C, A, B = "cd4svd5", "mean3", "majority:raw:all:nm"
SIGMA = P.SIGMA_PDS
GATE_SHAPE, GATE_AMPL = P.GATE_SHAPE, P.GATE_AMPL
md: list[str] = []
out: dict = {}


def Pr(s: str = "") -> None:
    print(s, flush=True)
    md.append(s)


def f(x, nd=5, sign=True):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    return f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"


def pct(x, nd=2):
    return "—" if x is None or not np.isfinite(x) else f"{100 * x:+.{nd}f} %"


# ------------------------------------------------------------------------------------------ load
ap = argparse.ArgumentParser()
ap.add_argument("--tag", default="", choices=("", "m"), help="'' = the mandate's equal three-donor mean; 'm' = the matched-kolf2 sensitivity (cancer folds)")
TAG = ap.parse_args().tag
rows, pt, idn, picks, cov, meta = {}, {}, {}, {}, {}, {}
for h in FOLDS:
    t = TAG if h != "kolf2" else ""                       # the kolf2 fold has no third donor; always the letter run
    d = json.loads((WORK / f"sweep264{t}-{h}.json").read_text())
    meta[h] = d["meta"]
    for r in d["rows"]:
        rows[(h, r["arm"])] = r
    z = np.load(WORK / f"pertarget264{t}-{h}.npz")
    for k in z.files:
        pt[(h, k)] = z[k]
    z = np.load(WORK / f"identity264{t}-{h}.npz")
    for k in z.files:
        idn[(h, k)] = z[k]
    picks[h] = dict(np.load(WORK / f"picks264{t}-{h}.npz"))
    cov[h] = dict(np.load(WORK / f"cov264{t}-{h}.npz"))
Pr(f"# 264 analysis — {'the MATCHED-kolf2 sensitivity (kolf2 scaled to the production donors\' whitened amplitude)' if TAG else 'the mandate\'s construction (equal three-donor mean)'}\n")
targets = {h: [str(x) for x in meta[h]["targets"]] for h in FOLDS}
arms_cancer = [a for a in dict.fromkeys(a for (_, a) in rows) if all((h, a) in rows for h in CANCER)]
arms_kolf = [a for (h, a) in rows if h == "kolf2"]


# ------------------------------------------------------------------------------------------ helpers
def nmae_concat(arm, folds):
    U = np.concatenate([idn[(h, f"{arm}|U")] for h in folds])
    M = np.concatenate([idn[(h, f"{arm}|M")] for h in folds])
    D = np.concatenate([idn[(h, f"{arm}|D")] for h in folds])
    return float(np.mean(1 + (U - 2 * M) / D))


def compare(arm, base, folds):
    ra, rm, dp, dn_a, dn_b, ns, dnm = [], [], [], [], [], [], []
    for h in folds:
        a, b = rows[(h, arm)], rows[(h, base)]
        ra.append(a["cos_agg"] / b["cos_agg"] - 1.0)
        rm.append(a["a23_cos"] / b["a23_cos"] - 1.0)
        dp.append(a["pds"] - b["pds"])
        dn_a.append(a["delivered_norm"])
        dn_b.append(b["delivered_norm"])
        ns.append(a["n_selected"])
        dnm.append(a["nmae_wald89"] - b["nmae_wald89"])
    nm_a, nm_b = nmae_concat(arm, folds), nmae_concat(base, folds)
    return dict(arm=arm, base=base, n=len(folds), rel_agg=float(np.mean(ra)), rel_macro=float(np.mean(rm)),
                d_pds=float(np.mean(dp)), per_fold_d_pds=[float(x) for x in dp], per_fold_rel_agg=[float(x) for x in ra],
                per_fold_rel_macro=[float(x) for x in rm], per_fold_d_nmae89=[float(x) for x in dnm],
                pds_pos=int(sum(x > 0 for x in dp)), sigma=float(np.mean(dp)) / SIGMA,
                fold_se=float(np.std(dp, ddof=1) / np.sqrt(len(dp))) if len(dp) > 1 else float("nan"),
                delivered_ratio=float(np.mean(dn_a) / np.mean(dn_b)), per_fold_delivered_ratio=[float(x / y) for x, y in zip(dn_a, dn_b)],
                nmae_arm=nm_a, nmae_base=nm_b, d_nmae=nm_a - nm_b,
                d_prec89=float(np.mean([rows[(h, arm)]["top89_precision"] - rows[(h, base)]["top89_precision"] for h in folds])),
                n_selected=float(np.mean(ns)), per_fold_n_selected=[int(x) for x in ns])


def price_cell(c):
    sc = Q.scenarios(c)
    amp = abs(c["delivered_ratio"] - 1.0)
    cls = "shape-only" if amp <= P.AMP_SHAPE_BOUND else "amplitude-changing"
    gate = GATE_SHAPE if cls == "shape-only" else GATE_AMPL
    return dict(pessimistic=sc["PESSIMISTIC"]["total"], point=sc["POINT"]["total"], optimistic=sc["OPTIMISTIC"]["total"],
                legs_pess={k: sc["PESSIMISTIC"][k] for k in ("pds", "mse", "de3", "nmae", "draw_amp", "dcos_board")},
                legs_point={k: sc["POINT"][k] for k in ("pds", "mse", "de3", "nmae", "draw_amp", "dcos_board")},
                legs_opt={k: sc["OPTIMISTIC"][k] for k in ("pds", "mse", "de3", "nmae")},
                factors_pess={k: sc["PESSIMISTIC"][k] for k in ("k_pds", "k_cos", "de3_rate", "x_loss", "x_gain")},
                amp_leg_applied=sc["POINT"]["amp_leg_applied"], cls=cls, gate=gate, gate_x=sc["PESSIMISTIC"]["total"] / gate)


def retention(arm, base, folds):
    return Q.retention(pt, arm, base, folds, targets, n=3)


def group_split(arm, base, folds, key="base_c_cos_truth"):
    """Per-target d pds by 262's groups on the cached base: (e) base·truth <= 0, marginal (0, 0.05], saturated > 0.05;
    and by the selected / unselected split."""
    res = {}
    acc = {g: [] for g in ("e", "band", "sat", "selected", "unselected")}
    for h in folds:
        d = pt[(h, f"{arm}|pds")] - pt[(h, f"{base}|pds")]
        bc = cov[h][key]
        p = picks[h].get(arm)
        sel = (p >= 0) if p is not None else np.zeros(d.size, bool)
        acc["e"].append(d[bc <= 0]); acc["band"].append(d[(bc > 0) & (bc <= 0.05)]); acc["sat"].append(d[bc > 0.05])
        acc["selected"].append(d[sel]); acc["unselected"].append(d[~sel])
    for g, v in acc.items():
        x = np.concatenate(v) if v else np.zeros(0)
        res[g] = dict(n=int(x.size), sum=float(x.sum()), mean=float(x.mean()) if x.size else float("nan"),
                      gainers=int((x > 1e-12).sum()), losers=int((x < -1e-12).sum()))
    return res


def dropped_hist(arm, folds):
    hist = {}
    for h in folds:
        p = picks[h].get(arm)
        if p is None:
            continue
        dn = meta[h]["donors"]
        for j in range(len(dn)):
            name = "kolf2" if dn[j] == "kolf2" else dn[j]
            hist[name] = hist.get(name, 0) + int((p == j).sum())
    return hist


def by_dropped(arm, base, folds):
    res = {}
    for h in folds:
        d = pt[(h, f"{arm}|pds")] - pt[(h, f"{base}|pds")]
        p = picks[h].get(arm)
        if p is None:
            continue
        for j, dn in enumerate(meta[h]["donors"]):
            m = p == j
            r = res.setdefault(dn, dict(n=0, sum=0.0))
            r["n"] += int(m.sum()); r["sum"] += float(d[m].sum())
    for r in res.values():
        r["mean"] = r["sum"] / r["n"] if r["n"] else float("nan")
    return res


def agree_with_oracle(arm, folds, oracle="oracle_drop:wht:all:nm"):
    ag = tot = 0
    for h in folds:
        p, o = picks[h].get(arm), picks[h].get(oracle)
        if p is None or o is None:
            continue
        m = (p >= 0) & (o >= 0)
        ag += int((p[m] == o[m]).sum()); tot += int(m.sum())
    return (ag / tot) if tot else float("nan"), tot


# ================================================================================ KA
Pr("## KA. Known-answer rows per fold\n")
Pr("| fold | KA rows | row KAs | failures | rows | seconds | peak RSS GiB | donors | three-donor targets (all / tau0) | (a) rows changed vs cached pre | kolf2 cols outside operator (energy) |")
Pr("|---|---:|---:|---:|---:|---:|---:|---|---|---:|---|")
ka_all_ok = True
for h in FOLDS:
    m = meta[h]
    rk = [r for r in rows.values() if r["fold"] == h and "ka" in r]
    ok = m["ka_pass"] == m["ka_total"] and not m["failures"] and all(r["ka"]["ok"] for r in rk)
    ka_all_ok &= ok
    a3 = m["aud3"]
    extra = "—" if h == "kolf2" else f"{a3['n_cols_outside_operator']}/{a3['n_carried_cols']} ({100*a3['energy_frac_outside_operator']:.3f} %)"
    Pr(f"| `{h}` | {m['ka_pass']}/{m['ka_total']} | {sum(r['ka']['ok'] for r in rk)}/{len(rk)} (max\\|d\\| {max((r['ka']['max_abs_diff'] for r in rk), default=0.0):.1e}) | "
       f"{len(m['failures'])} | {m['n_rows']} | {m['seconds']} | {m['peak_rss_gib']:.2f} | {', '.join(m['donors'])} | {m['gates']['all']} / {m['gates']['tau0']} | "
       f"{a3['n_changed_vs_cached_pre']} | {extra} |")
out["ka_all_ok"] = ka_all_ok
out["ka_rows"] = {h: meta[h]["ka_rows"] for h in FOLDS}
Pr(f"\nAll known-answer rows pass: **{ka_all_ok}**.\n")
Pr("R3 per cancer fold (kolf2 prepared as a DONOR for the fold's targets vs the kolf2 fold's cached TRUTH rows, float32):\n")
Pr("| fold | targets kolf2 carries | in common with the kolf2 fold | rows bit-identical | max\\|d\\| | sha256 mine / kolf2-fold truth | kolf2 weight | raw median norm |")
Pr("|---|---:|---:|---:|---:|---|---:|---:|")
for h in CANCER:
    a3 = meta[h]["aud3"]
    Pr(f"| `{h}` | {a3['n_kolf_carried']} | {a3['n_common']} | {a3['rows_bit_identical']} | {a3['max_abs_diff_float32']:.1e} | `{a3['sha256_mine'][:16]}…` / `{a3['sha256_kolf2_fold_truth'][:16]}…` | {a3['kolf2_weight']:.6f} | {a3['raw_median_norm']:.6f} |")
out["r3"] = {h: {k: v for k, v in meta[h]["aud3"].items()} for h in CANCER}
Pr()
Pr("## W. Whitening audit: how the fold's operator treats each donor's rows (raw median norm on the operator columns; amplification = un-preserved whitened median / raw median; energy on frame A's CPM gate before / after whitening; energy in the operator's bottom eigen-decile after whitening; delivered / effect norm of the rows alone)\n")
Pr("| fold | rows | n | raw median | amplification | energy on gate A raw / whitened | bottom eigen-decile (whitened) | delivered / effect |")
Pr("|---|---|---:|---:|---:|---|---:|---:|")
for h in FOLDS:
    for nm, v in meta[h]["whitening_audit"].items():
        Pr(f"| `{h}` | {nm} | {v['n']} | {v['raw_median']:.3f} | ×{v['amplification']:.2f} | {100*v['energy_on_gateA_raw']:.1f} % / {100*v['energy_on_gateA_wht']:.1f} % | {100*v['energy_bottom_eigen_decile_wht']:.1f} % | {v['delivered_over_effect']:.3f} |")
    if h != "kolf2":
        Pr(f"| `{h}` | kolf2 weight used | | | {meta[h]['aud3']['kolf2_weight']:.6f} ({meta[h]['aud3'].get('kolf_weight_mode')}) | | | |")
out["whitening_audit"] = {h: meta[h]["whitening_audit"] for h in FOLDS}
Pr()

# ================================================================================ per-fold table
Pr("## T1. Every row per fold: vs (c) `cd4svd5` and vs (a) `mean3` (Δ`pds`, rel `cos_agg`, rel macro, Δ`nmae89`, delivered ratio)\n")
Pr("| fold | arm | n sel | chg vs a / c | `pds` | Δ vs c | Δ vs a | rel agg vs c | rel agg vs a | rel macro vs c | Δ`nmae89` vs c | deliv vs c | deliv vs a |")
Pr("|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
per_fold_rows = {}
for h in FOLDS:
    rc, ra = rows[(h, C)], rows[(h, A)]
    for arm in [a for (hh, a) in rows if hh == h]:
        r = rows[(h, arm)]
        rec = dict(pds=r["pds"], d_pds_c=r["pds"] - rc["pds"], d_pds_a=r["pds"] - ra["pds"], rel_agg_c=r["cos_agg"] / rc["cos_agg"] - 1,
                   rel_agg_a=r["cos_agg"] / ra["cos_agg"] - 1, rel_macro_c=r["a23_cos"] / rc["a23_cos"] - 1, d_nmae_c=r["nmae_wald89"] - rc["nmae_wald89"],
                   d_nmae_a=r["nmae_wald89"] - ra["nmae_wald89"], deliv_c=r["delivered_norm"] / rc["delivered_norm"], deliv_a=r["delivered_norm"] / ra["delivered_norm"],
                   n_selected=r["n_selected"], chg_a=r["n_rows_changed_vs_a"], chg_c=r["n_rows_changed_vs_c"],
                   cd4_fraction=r["audit"]["cd4_fraction"], union_fraction=r["audit"]["union_fraction"], base_median_norm=r["base_median_norm_pre_wcols"])
        per_fold_rows[(h, arm)] = rec
        Pr(f"| `{h}` | `{arm}` | {r['n_selected']} | {r['n_rows_changed_vs_a']} / {r['n_rows_changed_vs_c']} | {r['pds']:.6f} | {f(rec['d_pds_c'])} | {f(rec['d_pds_a'])} | "
           f"{pct(rec['rel_agg_c'])} | {pct(rec['rel_agg_a'])} | {pct(rec['rel_macro_c'])} | {f(rec['d_nmae_c'])} | {rec['deliv_c']:.4f} | {rec['deliv_a']:.4f} |")
out["per_fold_rows"] = {f"{h}|{a}": v for (h, a), v in per_fold_rows.items()}
Pr()

# ================================================================================ the increments on the cancer folds
Pr("## T2. The increments on the THREE CANCER FOLDS (kolf2 the third donor): every arm vs (c) and vs (a), priced\n")
Pr("| arm | base | n sel | Δ`pds` (3-fold) | per fold | +/3 | σ | rel `cos_agg` | rel macro | Δ`nmae89` (identity, concat) | deliv ratio | class | PESSIMISTIC | POINT | OPTIMISTIC | LT3O retained | top-3 share |")
Pr("|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|")
inc = {}
for arm in arms_cancer:
    for base in (C, A):
        if arm == base:
            continue
        c = compare(arm, base, CANCER)
        pr = price_cell(c)
        ret = retention(arm, base, CANCER)
        inc[(arm, base)] = dict(compare=c, price=pr, retention=ret)
        Pr(f"| `{arm}` | {base} | {c['n_selected']:.0f} | **{f(c['d_pds'])}** | {', '.join(f(x, 4) for x in c['per_fold_d_pds'])} | {c['pds_pos']}/3 | {c['sigma']:+.1f} | "
           f"{pct(c['rel_agg'])} | {pct(c['rel_macro'])} | {f(c['d_nmae'])} | {c['delivered_ratio']:.4f} | {pr['cls']} | **{f(pr['pessimistic'])}** | {f(pr['point'])} | {f(pr['optimistic'])} | "
           f"{100*ret['retained']:+.0f} % | {100*ret['share']:+.0f} % |")
out["increments_cancer"] = {f"{a}|vs|{b}": v for (a, b), v in inc.items()}
Pr()

# ================================================================================ kolf2 apart
Pr("## T3. `kolf2` apart (the fold's own three production donors, the production operator; R2 = 263's row at 0.0e+00)\n")
Pr("| arm | n sel | Δ`pds` vs c | rel `cos_agg` | rel macro | Δ`nmae89` | deliv ratio | LT3O retained | top-3 share | KA vs 263 |")
Pr("|---|---:|---:|---:|---:|---:|---:|---:|---:|---|")
kolf = {}
for arm in arms_kolf:
    if arm == C:
        continue
    c = compare(arm, C, ("kolf2",))
    ret = retention(arm, C, ("kolf2",))
    r = rows[("kolf2", arm)]
    ka = r.get("ka")
    kolf[arm] = dict(compare=c, retention=ret, ka=ka)
    Pr(f"| `{arm}` | {r['n_selected']} | **{f(c['d_pds'])}** | {pct(c['rel_agg'])} | {pct(c['rel_macro'])} | {f(c['d_nmae'])} | {c['delivered_ratio']:.4f} | {100*ret['retained']:+.0f} % | {100*ret['share']:+.0f} % | "
       f"{('PASS ' + format(ka['max_abs_diff'], '.1e')) if ka and ka['ok'] else ('FAIL ' + format(ka['max_abs_diff'], '.1e') if ka else '—')} |")
out["kolf2"] = kolf
Pr()

# ================================================================================ the headline cells in full
Pr("## T4. The headline cells in full (cancer folds)\n")
head = {}
for arm in (B, "majority:raw:all", "majority:raw:all:nmc", "majority:raw:tau0:nm", "majority:wht:all:nm", "mean3", "mean3:nm", "mean3:unrestricted"):
    if arm not in arms_cancer:
        continue
    for base in (C, A):
        if arm == base or (arm.startswith("mean3") and base == A):
            continue
        c = inc[(arm, base)]["compare"]
        pr = inc[(arm, base)]["price"]
        ret = inc[(arm, base)]["retention"]
        gs = group_split(arm, base, CANCER)
        hist = dropped_hist(arm, CANCER)
        byd = by_dropped(arm, base, CANCER)
        agr, ntot = agree_with_oracle(arm, CANCER)
        head[f"{arm}|vs|{base}"] = dict(compare=c, price=pr, retention=ret, groups=gs, dropped_hist=hist, by_dropped=byd, oracle_agreement=agr, oracle_n=ntot)
        Pr(f"### `{arm}` vs `{base}`\n")
        Pr(f"- Δ`pds` **{f(c['d_pds'])}** (per fold {', '.join(f(x, 4) for x in c['per_fold_d_pds'])}; {c['pds_pos']}/3; {c['sigma']:+.1f} LOO-σ; fold SE {c['fold_se']:.5f}); "
           f"rel `cos_agg` {pct(c['rel_agg'])} (per fold {', '.join(pct(x) for x in c['per_fold_rel_agg'])}); rel macro {pct(c['rel_macro'])}; Δ`nmae89` {f(c['d_nmae'])} (per fold {', '.join(f(x, 4) for x in c['per_fold_d_nmae89'])}); "
           f"delivered ratio {c['delivered_ratio']:.4f} (per fold {', '.join(format(x, '.4f') for x in c['per_fold_delivered_ratio'])}); top-89 precision {f(c['d_prec89'], 4)}.")
        Pr(f"- Price ({pr['cls']}, gate {pr['gate']:+.4f}): PESSIMISTIC **{f(pr['pessimistic'])}** ({pr['gate_x']:+.2f}× gate; legs pds {f(pr['legs_pess']['pds'])} / mse {f(pr['legs_pess']['mse'])} / DE-3 {f(pr['legs_pess']['de3'])} / nmae {f(pr['legs_pess']['nmae'])}; "
           f"k_pds {pr['factors_pess']['k_pds']}, k_cos {pr['factors_pess']['k_cos']}, Δcos_board {pr['legs_pess']['dcos_board']:+.5f}, amp leg {'applied' if pr['amp_leg_applied'] else 'not applied'}), "
           f"POINT {f(pr['point'])} (pds {f(pr['legs_point']['pds'])} / mse {f(pr['legs_point']['mse'])} / DE-3 {f(pr['legs_point']['de3'])} / nmae {f(pr['legs_point']['nmae'])}), OPTIMISTIC {f(pr['optimistic'])}.")
        Pr(f"- Leave-top-3-out retention **{100*ret['retained']:+.0f} %** (per fold {', '.join(format(100*x, '+.0f') + ' %' for x in ret['per_fold_retained'])}); top-3 share {100*ret['share']:+.0f} %; "
           f"top-3 per fold: " + "; ".join(f"{LINE[h]}: {', '.join(ret['top_names'][h])}" for h in CANCER) + ".")
        Pr(f"- By 262's group on the cached base (base·truth): (e) ≤ 0 n {gs['e']['n']} sum {f(gs['e']['sum'], 3)} mean {f(gs['e']['mean'], 4)}; marginal (0, 0.05] n {gs['band']['n']} sum {f(gs['band']['sum'], 3)} mean {f(gs['band']['mean'], 4)}; "
           f"saturated > 0.05 n {gs['sat']['n']} sum {f(gs['sat']['sum'], 3)}; selected n {gs['selected']['n']} sum {f(gs['selected']['sum'], 3)} gainers {gs['selected']['gainers']} / losers {gs['selected']['losers']}; "
           f"unselected n {gs['unselected']['n']} sum {f(gs['unselected']['sum'], 3)}.")
        if hist:
            Pr(f"- Dropped donor counts {hist}; Δ`pds` sum by dropped donor: " + "; ".join(f"{k}: n {v['n']} sum {f(v['sum'], 3)} mean {f(v['mean'], 4)}" for k, v in byd.items()) +
               f"; drops the same donor as the whitened-truth oracle on {100*agr:.0f} % of {ntot} (chance 33 %).")
        Pr()
out["headline"] = head

# ================================================================================ controls below half?
Pr("## T5. Controls against the arm (cancer folds; Δ`pds` vs the same base)\n")
ctrl_names = {"reverse": ["reverse:raw:all:nm"], "random": ["random:all:s0:nm", "random:all:s1:nm", "random:all:s2:nm"],
              "kolf2 ALONE on the gated targets (information, not a verdict control)": ["only:kolf2:all:nm"],
              "drop:kolf2 (= two production donors, nm)": ["drop:kolf2:all:nm"], "drop:kolf2 (as-is = cd4svd5 rows on the gated targets)": ["drop:kolf2:all"],
              "drop:production donor (nm)": [a for a in arms_cancer if a.startswith("drop:") and a.endswith(":nm") and "kolf2" not in a]}
ctrl_tau = {"reverse tau0": ["reverse:raw:tau0:nm"], "random tau0": ["random:tau0:s0:nm", "random:tau0:s1:nm", "random:tau0:s2:nm"]}
verdict_ctrl = {}
for base in (C, A):
    Pr(f"### vs `{base}`\n")
    Pr("| control | cells | Δ`pds` (mean over cells) | per cell | arm Δ`pds` | below half the arm? |")
    Pr("|---|---:|---:|---|---:|---|")
    arm_d = inc[(B, base)]["compare"]["d_pds"]
    allok = True
    for name, cells in ctrl_names.items():
        cells = [c for c in cells if c in arms_cancer]
        if not cells:
            continue
        vals = [inc[(c, base)]["compare"]["d_pds"] for c in cells]
        ok = all(v < 0.5 * arm_d for v in vals) if arm_d > 0 else False
        if "information" not in name:
            allok &= ok
        Pr(f"| {name} | {len(cells)} | {f(float(np.mean(vals)))} | {', '.join(f(v, 4) for v in vals)} | {f(arm_d)} | {'yes' if ok else 'NO'} |")
    arm_t = inc[("majority:raw:tau0:nm", base)]["compare"]["d_pds"] if ("majority:raw:tau0:nm", base) in inc else float("nan")
    for name, cells in ctrl_tau.items():
        cells = [c for c in cells if c in arms_cancer]
        if not cells:
            continue
        vals = [inc[(c, base)]["compare"]["d_pds"] for c in cells]
        ok = all(v < 0.5 * arm_t for v in vals) if arm_t > 0 else False
        Pr(f"| {name} (vs the tau0 arm) | {len(cells)} | {f(float(np.mean(vals)))} | {', '.join(f(v, 4) for v in vals)} | {f(arm_t)} | {'yes' if ok else 'NO'} |")
    verdict_ctrl[base] = allok
    Pr()
out["controls_below_half"] = verdict_ctrl

# ================================================================================ the oracle references
Pr("## T6. Ceiling references (drop the donor whose row has the lowest cosine to the held truth; never an arm)\n")
for arm in ("oracle_drop:raw:all:nm", "oracle_drop:wht:all:nm"):
    if arm in arms_cancer:
        for base in (C, A):
            c = inc[(arm, base)]["compare"]
            Pr(f"- `{arm}` vs `{base}`: Δ`pds` **{f(c['d_pds'])}** ({', '.join(f(x, 4) for x in c['per_fold_d_pds'])}), rel `cos_agg` {pct(c['rel_agg'])}, delivered {c['delivered_ratio']:.4f}; "
               f"the majority rule captures {100 * inc[(B, base)]['compare']['d_pds'] / c['d_pds']:+.0f} % of it.")
Pr()

# ================================================================================ verdict
cb, ca = inc[(B, C)], inc[(B, A)]
# 12:50Z fix: the gate is class-dependent (shape-only +0.0019 / amplitude-changing +0.0089); the first cut compared against GATE_SHAPE.
cond = dict(pess_ge_gate=cb["price"]["pessimistic"] >= cb["price"]["gate"], pos_3of3=cb["compare"]["pds_pos"] == 3, two_sigma=cb["compare"]["d_pds"] >= 2 * SIGMA,
            retention_50=cb["retention"]["retained"] >= 0.5 and cb["compare"]["d_pds"] > 0, controls_half=verdict_ctrl[C])
cond_a = dict(pess_ge_gate=ca["price"]["pessimistic"] >= ca["price"]["gate"], pos_3of3=ca["compare"]["pds_pos"] == 3, two_sigma=ca["compare"]["d_pds"] >= 2 * SIGMA,
              retention_50=ca["retention"]["retained"] >= 0.5 and ca["compare"]["d_pds"] > 0, controls_half=verdict_ctrl[A])
clears = all(cond.values())
out["conditions_b_minus_c"] = cond
out["conditions_b_minus_a"] = cond_a
out["clears"] = clears
Pr("## Verdict\n")
Pr(f"(b) − (c) `majority:raw:all:nm` vs `cd4svd5`: Δ`pds` {f(cb['compare']['d_pds'])} ({cb['compare']['pds_pos']}/3, {cb['compare']['sigma']:+.1f} σ), PESSIMISTIC {f(cb['price']['pessimistic'])} vs gate {cb['price']['gate']:+.4f} "
   f"({cb['price']['cls']}), retention {100*cb['retention']['retained']:+.0f} %, controls below half: {verdict_ctrl[C]} → conditions {cond}")
Pr(f"(b) − (a) `majority:raw:all:nm` vs `mean3`: Δ`pds` {f(ca['compare']['d_pds'])} ({ca['compare']['pds_pos']}/3, {ca['compare']['sigma']:+.1f} σ), PESSIMISTIC {f(ca['price']['pessimistic'])} "
   f"({ca['price']['cls']}), retention {100*ca['retention']['retained']:+.0f} %, controls below half: {verdict_ctrl[A]} → conditions {cond_a}")
Pr(f"\n**{'CLEARS' if clears else 'DOES NOT CLEAR'}**\n")
(WORK / f"analysis264{TAG}.json").write_text(json.dumps(out, indent=1, default=float))
(WORK / f"tables264{TAG}.md").write_text("\n".join(md) + "\n")
(WORK / f"analyse264{TAG}.log").write_text("\n".join(md) + "\n")
