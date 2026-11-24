"""263 -- the analysis: R1-R3, the arbiter x gating table with its controls, the oracle-headroom fractions, the
marginal-band crossings, leave-top-3-out retention, and the three-scenario price (261's `price()` / `factors()`,
257's constants) of every cell against `cd4svd5`.  Reads work/sweep263-<fold>.json and the npz siblings, 262's
diag262-<fold>.npz (read-only) and pilot262.json (R3); the fold caches under work/plolo (byte copies).
Writes work/analysis263.json and work/tables263.md.

    uv run python analyse263.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E262 = ROOT / "experiments/262-where-pds-is-lost"
E261 = ROOT / "experiments/261-union-deflation"
E257 = ROOT / "experiments/257-soft-deflation-price"
sys.path.insert(0, str(E261))
sys.path.insert(0, str(E257))
import price257 as P  # noqa: E402  (constants, raw_for, scaled_terms; read-only at import)
import price261 as Q  # noqa: E402  (price, factors, retention; read-only at import)

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
FOLDS = CANCER + ("kolf2",)
LINE = {"k562_gwps": "K562", "orion_hct116": "HCT116", "orion_hek293t": "HEK293T", "kolf2": "KOLF2"}
BASE = "cd4svd5"
GATES = ("tau0", "tau0.1", "all")
SHUF_SEED = 20726
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


# ------------------------------------------------------------------------------------------ load
rows, pt, idn, picks, cov, meta = {}, {}, {}, {}, {}, {}
for h in FOLDS:
    d = json.loads((WORK / f"sweep263-{h}.json").read_text())
    meta[h] = d["meta"]
    for r in d["rows"]:
        rows[(h, r["arm"])] = r
    z = np.load(WORK / f"pertarget263-{h}.npz")
    for k in z.files:
        pt[(h, k)] = z[k]
    z = np.load(WORK / f"identity263-{h}.npz")
    for k in z.files:
        idn[(h, k)] = z[k]
    picks[h] = dict(np.load(WORK / f"picks263-{h}.npz"))
    cov[h] = dict(np.load(WORK / f"cov263-{h}.npz"))
arms_all = sorted({a for (_, a) in rows})
arms_cancer = [a for a in arms_all if all((h, a) in rows for h in CANCER)]
targets = {h: [str(x) for x in meta[h]["targets"]] for h in FOLDS}

# ================================================================================ 0. known-answer rows
Pr("## KA. Known-answer rows per fold\n")
Pr("| fold | KA rows | row KAs | failures | rows | seconds | peak RSS GiB | gates tau0 / tau0.1 / all | donor-donor raw cos median (whitened) |")
Pr("|---|---:|---:|---:|---:|---:|---:|---|---|")
ka_all_ok = True
for h in FOLDS:
    m = meta[h]
    rk = [r for r in rows.values() if r["fold"] == h and "ka" in r]
    ok = m["ka_pass"] == m["ka_total"] and not m["failures"] and all(r["ka"]["ok"] for r in rk)
    ka_all_ok &= ok
    g = m["gates"]
    Pr(f"| `{h}` | {m['ka_pass']}/{m['ka_total']} | {sum(r['ka']['ok'] for r in rk)}/{len(rk)} (max\\|d\\| {max(r['ka']['max_abs_diff'] for r in rk):.1e}) | "
       f"{len(m['failures'])} | {m['n_rows']} | {m['seconds']} | {m['peak_rss_gib']:.2f} | {g['tau0']} / {g['tau0.1']} / {g['all']} | "
       f"{m['dd_raw_median']:+.4f} ({m['dd_wht_median']:+.4f}) |")
out["ka_all_ok"] = ka_all_ok
Pr(f"\nAll known-answer rows pass: **{ka_all_ok}**.\n")

# ================================================================================ R3. 262 section 7.2 reproduced
Pr("## R3. 262 §7.2 (CD4 as arbiter, cross-line truth cosines) reproduced from the fold caches, then read in `pds`\n")
cache = {h: np.load(WORK / "plolo" / f"fold-{h}.npz", allow_pickle=True) for h in FOLDS}
truths = {h: np.asarray(cache[h]["truth"], dtype=np.float64) for h in FOLDS}
scored = {h: np.asarray(cache[h]["scored"]) for h in FOLDS}
tlist = {h: [str(x) for x in cache[h]["targets"]] for h in FOLDS}
tpos = {h: {t: i for i, t in enumerate(tlist[h])} for h in FOLDS}
d262 = {h: np.load(E262 / f"work/diag262-{h}.npz") for h in FOLDS}


def unit_rows(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


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
pilot = json.loads((E262 / "work/pilot262.json").read_text())["b"]
r3 = {}
r3_max = 0.0
Pr("| fold, rows | n | mean of two | best single (oracle) | CD4-arbiter pick | random pick | > 0.05: mean / best / arbiter | arbiter = best | max\\|d\\| vs pilot262.json |")
Pr("|---|---:|---:|---:|---:|---:|---|---:|---:|")
for h in CANCER:
    t = np.asarray(tlist[h])
    donors = meta[h]["sources"]
    A = np.array([[pair_cos[(h, s)].get(x, np.nan) for s in donors] for x in t])
    C = np.array([[d262[s]["cov|cd4_cos_truth"][tpos[s][x]] if x in tpos[s] else np.nan for s in donors] for x in t])
    base_ct = d262[h]["cov|base_cos_truth"]
    ok = np.isfinite(A).all(1) & np.isfinite(C).all(1)
    mean_ag, best = A.mean(1), A.max(1)
    pick = np.argmax(C, axis=1)
    arb = A[np.arange(len(t)), pick]
    rng = np.random.default_rng(SHUF_SEED)
    rnd = A[np.arange(len(t)), rng.integers(0, 2, len(t))]
    e = base_ct <= 0
    band = (base_ct > 0) & (base_ct <= 0.05)
    rho = spearmanr(mean_ag[ok], base_ct[ok]).correlation
    res = {}
    for gname, gm in (("all", ok), ("group e (base·truth <= 0)", ok & e), ("marginal band", ok & band)):
        res[gname] = dict(n=int(gm.sum()), mean_of_two=float(mean_ag[gm].mean()), best_single=float(best[gm].mean()),
                          cd4_arbiter=float(arb[gm].mean()), random_pick=float(rnd[gm].mean()),
                          frac_best_above_0_05=float((best[gm] > 0.05).mean()), frac_mean_above_0_05=float((mean_ag[gm] > 0.05).mean()),
                          frac_arbiter_above_0_05=float((arb[gm] > 0.05).mean()),
                          arbiter_picks_best=float((pick[gm] == np.argmax(A, axis=1)[gm]).mean()))
        ref = pilot[h]["groups"][gname]
        dmax = max(abs(res[gname][k] - ref[k]) for k in res[gname])
        r3_max = max(r3_max, dmax)
        Pr(f"| {LINE[h]}, {gname} | {res[gname]['n']} | {res[gname]['mean_of_two']:+.4f} | {res[gname]['best_single']:+.4f} | "
           f"{res[gname]['cd4_arbiter']:+.4f} | {res[gname]['random_pick']:+.4f} | {100*res[gname]['frac_mean_above_0_05']:.0f} % / "
           f"{100*res[gname]['frac_best_above_0_05']:.0f} % / {100*res[gname]['frac_arbiter_above_0_05']:.0f} % | "
           f"{100*res[gname]['arbiter_picks_best']:.0f} % | {dmax:.1e} |")
    r3_max = max(r3_max, abs(rho - pilot[h]["rho_mean_agreement_vs_base_cos_truth"]))
    r3[h] = dict(groups=res, rho=float(rho))
out["r3"] = dict(per_fold=r3, max_abs_diff_vs_pilot262=r3_max, ok=bool(r3_max < 1e-9))
Pr(f"\nR3 reproduces 262 §7.2 at max\\|d\\| **{r3_max:.1e}** → **{'PASS' if r3_max < 1e-9 else 'FAIL'}**.\n")


# ================================================================================ helpers
def nmae_concat(arm, folds):
    U = np.concatenate([idn[(h, f"{arm}|U")] for h in folds])
    M = np.concatenate([idn[(h, f"{arm}|M")] for h in folds])
    D = np.concatenate([idn[(h, f"{arm}|D")] for h in folds])
    return float(np.mean(1 + (U - 2 * M) / D))


def compare(arm, folds):
    ra, rm, dp, dq, dn_a, dn_b, ns = [], [], [], [], [], [], []
    for h in folds:
        a, b = rows[(h, arm)], rows[(h, BASE)]
        ra.append(a["cos_agg"] / b["cos_agg"] - 1.0)
        rm.append(a["a23_cos"] / b["a23_cos"] - 1.0)
        dp.append(a["pds"] - b["pds"])
        dq.append(a["a15_sgn_w"] - b["a15_sgn_w"])
        dn_a.append(a["delivered_norm"])
        dn_b.append(b["delivered_norm"])
        ns.append(a["n_selected"])
    nm_a, nm_b = nmae_concat(arm, folds), nmae_concat(BASE, folds)
    return dict(arm=arm, base=BASE, n=len(folds), rel_agg=float(np.mean(ra)), rel_macro=float(np.mean(rm)),
                d_pds=float(np.mean(dp)), per_fold_d_pds=[float(x) for x in dp], per_fold_rel_agg=[float(x) for x in ra],
                pds_pos=int(sum(x > 0 for x in dp)), sigma=float(np.mean(dp)) / SIGMA,
                fold_se=float(np.std(dp, ddof=1) / np.sqrt(len(dp))) if len(dp) > 1 else float("nan"),
                d_qw=float(np.mean(dq)), delivered_ratio=float(np.mean(dn_a) / np.mean(dn_b)),
                per_fold_delivered_ratio=[float(x / y) for x, y in zip(dn_a, dn_b)],
                nmae_arm=nm_a, nmae_base=nm_b, d_nmae=nm_a - nm_b,
                d_prec89=float(np.mean([rows[(h, arm)]["top89_precision"] - rows[(h, BASE)]["top89_precision"] for h in folds])),
                n_selected=float(np.mean(ns)), per_fold_n_selected=[int(x) for x in ns])


def band_stats(arm, folds):
    """Marginal-band crossings (262's base·truth in (0, 0.05]): own cosine crossing 0.02 upward, falling below 0;
    group-(e) rescues (base·truth <= 0 and the arm's own cosine > 0.02); gainers / losers among the selected targets."""
    up = down = resc = gain = lose = nsel = 0
    for h in folds:
        c = cov[h]
        b262 = c["base262_cos_truth"]
        band = (b262 > 0) & (b262 <= 0.05)
        ge = b262 <= 0
        cb, ca = pt[(h, f"{BASE}|cos")], pt[(h, f"{arm}|cos")]
        pb, pa = pt[(h, f"{BASE}|pds")], pt[(h, f"{arm}|pds")]
        up += int((band & (cb < 0.02) & (ca >= 0.02)).sum())
        down += int((band & (cb >= 0) & (ca < 0)).sum())
        resc += int((ge & (ca > 0.02)).sum())
        p = picks[h].get(arm)
        if p is not None:
            sel = p >= 0
            nsel += int(sel.sum())
            gain += int((sel & (pa > pb + 1e-12)).sum())
            lose += int((sel & (pa < pb - 1e-12)).sum())
    return dict(band_up=up, band_down=down, e_rescued=resc, sel_gainers=gain, sel_losers=lose, n_selected_pt=nsel)


def oracle_agreement(arm, folds):
    """Fraction of the arm's selected targets whose pick equals the whitened oracle's pick at the same gate."""
    gate = arm.split(":")[2] if arm.count(":") >= 2 else ("all" if ":" not in arm else arm.split(":")[1])
    agree = tot = 0
    for h in folds:
        p = picks[h].get(arm)
        o = picks[h].get(f"oracle_wht:{gate}")
        if p is None or o is None:
            continue
        m = (p >= 0) & (o >= 0)
        agree += int((p[m] == o[m]).sum())
        tot += int(m.sum())
    return float(agree / tot) if tot else float("nan")


def gate_of(arm):
    parts = arm.split(":")
    for g in ("tau0.1", "tau0", "all"):
        if g in parts:
            return g
    return "all"


def is_nm(arm):
    return arm.endswith(":nm")


def price_cell(c):
    sc = Q.scenarios(c)
    amp = abs(c["delivered_ratio"] - 1.0)
    cls = "shape-only" if amp <= P.AMP_SHAPE_BOUND else "amplitude-changing"
    gate = GATE_SHAPE if cls == "shape-only" else GATE_AMPL
    return dict(pessimistic=sc["PESSIMISTIC"]["total"], point=sc["POINT"]["total"], optimistic=sc["OPTIMISTIC"]["total"],
                legs_pess={k: sc["PESSIMISTIC"][k] for k in ("pds", "mse", "de3", "nmae", "draw_amp", "dcos_board")},
                legs_point={k: sc["POINT"][k] for k in ("pds", "mse", "de3", "nmae", "draw_amp", "dcos_board")},
                amp_leg_applied=sc["POINT"]["amp_leg_applied"], cls=cls, gate=gate, gate_x=sc["PESSIMISTIC"]["total"] / gate)


# ================================================================================ R2. the oracle ceiling
Pr("## R2. The oracle ceiling (per target the donor whose WHITENED row has the higher cosine to the held truth), vs `cd4svd5`\n")
Pr("| oracle | gate | norm | n selected (mean) | Δ`pds` (3-fold; per fold) | rel `cos_agg` | rel macro | Δ`nmae89` | delivered ratio | band ↑0.02 / ↓0 | (e) rescued | selected gainers / losers |")
Pr("|---|---|---|---:|---|---:|---:|---:|---:|---|---:|---|")
orc = {}
for name in ("oracle_wht", "oracle_raw", "oracle3_wht"):
    for g in GATES:
        for nm in ("", ":nm"):
            a = f"{name}:{g}{nm}"
            if a not in arms_cancer:
                continue
            c = compare(a, CANCER)
            b = band_stats(a, CANCER)
            orc[a] = dict(**c, **b)
            Pr(f"| `{name}` | {g} | {'matched' if nm else 'as-is'} | {c['n_selected']:.0f} | **{c['d_pds']:+.5f}** ({', '.join('%+.4f' % x for x in c['per_fold_d_pds'])}) | "
               f"{100*c['rel_agg']:+.2f} % | {100*c['rel_macro']:+.2f} % | {c['d_nmae']:+.5f} | {c['delivered_ratio']:.4f} | {b['band_up']} / {b['band_down']} | "
               f"{b['e_rescued']} | {b['sel_gainers']} / {b['sel_losers']} |")
out["oracle"] = orc
Pr("")

# ================================================================================ the arbiter x gating table
ARBS = ["cd4", "cd4d", "union", "kolf2", "norm", "selfcons", "basal"]
Pr("## A. Arbiters × gating (three cancer folds, vs `cd4svd5`; norm-matched = the selected row rescaled to the mean row's norm)\n")
Pr("| arbiter | space | gate | norm | n sel | Δ`pds` (3-fold) | per fold | +/3 | σ (0.001494) | rel `cos_agg` | rel macro | Δ`nmae89` | deliv. ratio | oracle headroom captured (same gate / all) | = oracle pick | band ↑0.02 / ↓0 | LT3O retention | REVERSE Δ`pds` | RANDOM Δ`pds` (3 seeds) |")
Pr("|---|---|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---|---:|---|---:|---:|---:|")
cells = {}
for arb in ARBS:
    for space in ("raw", "wht"):
        for g in GATES:
            for nm in ("", ":nm"):
                a = f"{arb}:{space}:{g}{nm}"
                if a not in arms_cancer:
                    continue
                c = compare(a, CANCER)
                b = band_stats(a, CANCER)
                rev = compare(f"{arb}:{space}:{g}:rev{nm}", CANCER)["d_pds"]
                rnd = [compare(f"random:{g}:s{k}{nm}", CANCER)["d_pds"] for k in range(3)]
                o_same = orc[f"oracle_wht:{g}{nm}"]["d_pds"]
                o_all = orc[f"oracle_wht:all{nm}"]["d_pds"]
                ret = Q.retention(pt, a, BASE, CANCER, targets, n=3)
                ag = oracle_agreement(a, CANCER)
                pr = price_cell(c)
                cells[a] = dict(**c, **b, rev_d_pds=rev, random_d_pds=rnd, random_mean=float(np.mean(rnd)), headroom_same=c["d_pds"] / o_same if o_same else float("nan"),
                                headroom_all=c["d_pds"] / o_all if o_all else float("nan"), oracle_agree=ag, retention=ret["retained"], top3_share=ret["share"],
                                top_names=ret["top_names"], price=pr)
                Pr(f"| `{arb}` | {space} | {g} | {'matched' if nm else 'as-is'} | {c['n_selected']:.0f} | **{c['d_pds']:+.5f}** | {', '.join('%+.4f' % x for x in c['per_fold_d_pds'])} | "
                   f"{c['pds_pos']}/3 | {c['sigma']:+.1f} | {100*c['rel_agg']:+.2f} % | {100*c['rel_macro']:+.2f} % | {c['d_nmae']:+.5f} | {c['delivered_ratio']:.4f} | "
                   f"{100*cells[a]['headroom_same']:+.0f} % / {100*cells[a]['headroom_all']:+.0f} % | {100*ag:.0f} % | {b['band_up']} / {b['band_down']} | "
                   f"{100*ret['retained']:+.0f} % | {rev:+.5f} | {np.mean(rnd):+.5f} |")
out["cells"] = cells
Pr("")

# softmax weighting
Pr("## B. Continuous weighting, w ∝ softmax(β · standardised arbiter score), gate all (three cancer folds)\n")
Pr("| arbiter | space | β | norm | Δ`pds` (3-fold) | per fold | +/3 | σ | rel `cos_agg` | rel macro | Δ`nmae89` | deliv. ratio | headroom (oracle all) | LT3O retention |")
Pr("|---|---|---:|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|")
sm = {}
for arb in ARBS:
    for space in ("raw", "wht"):
        for beta in (1, 3):
            for nm in ("", ":nm"):
                a = f"{arb}:{space}:all:sm{beta}{nm}"
                if a not in arms_cancer:
                    continue
                c = compare(a, CANCER)
                ret = Q.retention(pt, a, BASE, CANCER, targets, n=3)
                o_all = orc[f"oracle_wht:all{nm}"]["d_pds"]
                sm[a] = dict(**c, retention=ret["retained"], headroom_all=c["d_pds"] / o_all, price=price_cell(c))
                Pr(f"| `{arb}` | {space} | {beta} | {'matched' if nm else 'as-is'} | **{c['d_pds']:+.5f}** | {', '.join('%+.4f' % x for x in c['per_fold_d_pds'])} | {c['pds_pos']}/3 | "
                   f"{c['sigma']:+.1f} | {100*c['rel_agg']:+.2f} % | {100*c['rel_macro']:+.2f} % | {c['d_nmae']:+.5f} | {c['delivered_ratio']:.4f} | {100*sm[a]['headroom_all']:+.0f} % | {100*ret['retained']:+.0f} % |")
out["softmax"] = sm
Pr("")

# the controls on their own
Pr("## C. The controls on their own (three cancer folds, vs `cd4svd5`): what a per-target donor CHOICE costs when it is uninformed\n")
Pr("| control | gate | norm | Δ`pds` (3-fold) | per fold | rel `cos_agg` | rel macro | Δ`nmae89` | deliv. ratio |")
Pr("|---|---|---|---:|---|---:|---:|---:|---:|")
ctrl = {}
for g in GATES:
    for nm in ("", ":nm"):
        for k in range(3):
            a = f"random:{g}:s{k}{nm}"
            c = compare(a, CANCER)
            ctrl[a] = c
            Pr(f"| random s{k} | {g} | {'matched' if nm else 'as-is'} | {c['d_pds']:+.5f} | {', '.join('%+.4f' % x for x in c['per_fold_d_pds'])} | {100*c['rel_agg']:+.2f} % | "
               f"{100*c['rel_macro']:+.2f} % | {c['d_nmae']:+.5f} | {c['delivered_ratio']:.4f} |")
out["controls"] = ctrl
Pr("")

# ================================================================================ kolf2 apart
Pr("## D. `kolf2` apart (three production donors, the production operator): the majority rule and every arbiter, vs `cd4svd5` on that fold\n")
Pr("| arm | n sel | Δ`pds` | rel `cos_agg` | rel macro | Δ`nmae89` | deliv. ratio | REVERSE | RANDOM (3 seeds) | oracle (same gate) | headroom captured |")
Pr("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
kol = {}
for a in sorted(a for (h, a) in rows if h == "kolf2"):
    if a in (BASE, "mean_explicit") or ":rev" in a or a.startswith("random"):
        continue
    c = compare(a, ("kolf2",))
    g = gate_of(a)
    nm = ":nm" if is_nm(a) else ""
    rv = f"{a[:-3] if nm else a}:rev{nm}"
    rev = compare(rv, ("kolf2",))["d_pds"] if ("kolf2", rv) in rows else float("nan")
    rnd = float(np.mean([compare(f"random:{g}:s{k}{nm}", ("kolf2",))["d_pds"] for k in range(3)]))
    o = compare(f"oracle_wht:{g}{nm}", ("kolf2",))["d_pds"]
    kol[a] = dict(**c, rev=rev, random=rnd, oracle=o, headroom=c["d_pds"] / o if o else float("nan"))
    if a.startswith("oracle"):
        continue
    Pr(f"| `{a}` | {c['n_selected']:.0f} | **{c['d_pds']:+.5f}** | {100*c['rel_agg']:+.2f} % | {100*c['rel_macro']:+.2f} % | {c['d_nmae']:+.5f} | {c['delivered_ratio']:.4f} | "
       f"{rev:+.5f} | {rnd:+.5f} | {o:+.5f} | {100*kol[a]['headroom']:+.0f} % |")
out["kolf2"] = kol
Pr("")
Pr("`kolf2` oracle rows: " + "; ".join(f"`{a}` {kol[a]['d_pds']:+.5f} (n {kol[a]['n_selected']:.0f})" for a in kol if a.startswith("oracle")) + "\n")

# ================================================================================ price
Pr("## E. Price (READ-FIRST §3: coherent PESSIMISTIC / POINT / OPTIMISTIC; 261's `price()`; k_pds {0.85, 1.22}, k_cos {1.36, 3.76}, DE-three {0.0425, 0.24}, identity ×1.5 on a loss / ×0 on a gain; amplitude leg beyond 0.2 %; class = shape-only within 0.5 % delivered amplitude, else amplitude-changing with its +0.0089 gate)\n")
allcells = {**cells, **sm}
ranked = sorted(allcells, key=lambda a: -allcells[a]["price"]["pessimistic"])
Pr("| cell | Δ`pds` | +/3 | σ | class | PESSIMISTIC | POINT | OPTIMISTIC | × gate | legs pess (`pds` / `mse` / DE-3 / `nmae`) | LT3O | controls < ½ arm | CLEARS? |")
Pr("|---|---:|---:|---:|---|---:|---:|---:|---:|---|---:|---|---|")
verdict = {}
for a in ranked[:12] + [x for x in ranked if x.startswith(("cd4:raw", "kolf2:raw", "union:raw")) and x.endswith(":nm") and x not in ranked[:12]][:6]:
    c = allcells[a]
    pr = c["price"]
    ctrl_ok = (("rev_d_pds" in c) and c["d_pds"] > 0 and c["rev_d_pds"] < 0.5 * c["d_pds"] and c["random_mean"] < 0.5 * c["d_pds"])
    clears = (pr["pessimistic"] >= pr["gate"]) and (c["retention"] >= 0.5) and (c["pds_pos"] == 3) and (c["d_pds"] >= 2 * SIGMA) and ctrl_ok
    verdict[a] = dict(clears=bool(clears), controls_ok=bool(ctrl_ok), **pr)
    lp = pr["legs_pess"]
    Pr(f"| `{a}` | {c['d_pds']:+.5f} | {c['pds_pos']}/3 | {c['sigma']:+.1f} | {pr['cls']} | **{pr['pessimistic']:+.5f}** | {pr['point']:+.5f} | {pr['optimistic']:+.5f} | "
       f"{pr['gate_x']:+.2f}× | {lp['pds']:+.5f} / {lp['mse']:+.5f} / {lp['de3']:+.5f} / {lp['nmae']:+.5f} | {100*c['retention']:+.0f} % | "
       f"{'yes' if ctrl_ok else 'NO'} | {'**CLEARS**' if clears else 'no'} |")
out["verdict"] = verdict
out["ranked_top"] = ranked[:12]
n_clear = sum(v["clears"] for v in verdict.values())
Pr(f"\nCells clearing all conditions: **{n_clear}** of {len(allcells)} priced.\n")

# oracle prices for reference
Pr("### E.1 The oracle rows priced the same way (a ceiling, never an arm)\n")
Pr("| oracle | Δ`pds` | class | PESSIMISTIC | POINT | OPTIMISTIC |")
Pr("|---|---:|---|---:|---:|---:|")
for a, c in orc.items():
    pr = price_cell(c)
    Pr(f"| `{a}` | {c['d_pds']:+.5f} | {pr['cls']} | {pr['pessimistic']:+.5f} | {pr['point']:+.5f} | {pr['optimistic']:+.5f} |")
Pr("")

# ================================================================================ headroom summary
Pr("## F. Oracle-headroom fractions, summary (three cancer folds; the best cell per arbiter, either space, any gate, norm-matched or not)\n")
Pr("| arbiter | best cell | Δ`pds` | fraction of the oracle headroom at the same gate | at gate all | agreement with the oracle's pick | controls (rev / random) |")
Pr("|---|---|---:|---:|---:|---:|---|")
best_per = {}
for arb in ARBS:
    cand = [a for a in cells if a.startswith(arb + ":")]
    if not cand:
        continue
    a = max(cand, key=lambda x: cells[x]["d_pds"])
    c = cells[a]
    best_per[arb] = a
    Pr(f"| `{arb}` | `{a}` | {c['d_pds']:+.5f} | {100*c['headroom_same']:+.1f} % | {100*c['headroom_all']:+.1f} % | {100*c['oracle_agree']:.0f} % | {c['rev_d_pds']:+.5f} / {c['random_mean']:+.5f} |")
out["best_per_arbiter"] = best_per
maxfrac = max(cells[a]["headroom_same"] for a in cells)
maxfrac_all = max(cells[a]["headroom_all"] for a in cells)
Pr(f"\nLargest fraction of the oracle headroom captured by any observable arbiter cell: **{100*maxfrac:+.1f} %** (same gate), **{100*maxfrac_all:+.1f} %** (of the gate-all oracle).\n")
out["max_headroom_fraction"] = dict(same_gate=maxfrac, gate_all=maxfrac_all)

# agreement rates per fold (from the sweeps' own measurement)
Pr("## G. Arbiter agreement with the whitened oracle's pick, per fold (all targets with both defined)\n")
keys = sorted({k for h in FOLDS for k in meta[h]["agree_with_oracle"]})
Pr("| arbiter:space | " + " | ".join(LINE[h] for h in FOLDS) + " |")
Pr("|---|" + "---:|" * len(FOLDS))
for k in keys:
    Pr(f"| `{k}` | " + " | ".join((f"{100*meta[h]['agree_with_oracle'][k]['agree']:.0f} % (n {meta[h]['agree_with_oracle'][k]['n']})" if k in meta[h]["agree_with_oracle"] else "—") for h in FOLDS) + " |")
Pr("")


# ================================================================================ H. where the oracle headroom lives
Pr("## H. Where the oracle's headroom lives (three cancer folds pooled; per-target Δ`pds` of the oracle rows by 262's base·truth group)\n")
Pr("| oracle row | group | n | mean Δ`pds` | share of the oracle's gain | targets selected in the group |")
Pr("|---|---|---:|---:|---:|---:|")
hl = {}
for a in ("oracle_wht:all:nm", "oracle_wht:tau0:nm", "oracle3_wht:all:nm"):
    tot = sum(float((pt[(h, f"{a}|pds")] - pt[(h, f"{BASE}|pds")]).sum()) for h in CANCER)
    for gname, lo, hi in (("(e) base·truth <= 0", -2.0, 0.0), ("marginal band (0, 0.05]", 0.0, 0.05), ("saturated > 0.05", 0.05, 2.0)):
        n = 0; sd = 0.0; ns = 0
        for h in CANCER:
            b = cov[h]["base262_cos_truth"]; m = (b > lo) & (b <= hi)
            dd = pt[(h, f"{a}|pds")] - pt[(h, f"{BASE}|pds")]
            n += int(m.sum()); sd += float(dd[m].sum()); ns += int(((picks[h][a] >= 0) & m).sum())
        hl[f"{a}|{gname}"] = dict(n=n, mean=sd / n if n else float("nan"), share=sd / tot if tot else float("nan"), n_selected=ns)
        Pr(f"| `{a}` | {gname} | {n} | {sd/n if n else float('nan'):+.4f} | {100*sd/tot if tot else float('nan'):.0f} % | {ns} |")
out["headroom_by_group"] = hl
Pr("")

# ================================================================================ I. break-even pick accuracy
Pr("## I. What pick accuracy an arbiter would need (three cancer folds; from the random controls' per-target split into picks that "
   "AGREE with the whitened oracle and picks that do not)\n")
Pr("| gate | norm | n picks (3 seeds) | Δ`pds` per target, pick = oracle's | Δ`pds` per target, pick ≠ oracle's | random mean Δ`pds` (fold) | break-even accuracy p* | accuracy for +0.005 raw `pds` | best observable agreement |")
Pr("|---|---|---:|---:|---:|---:|---:|---:|---:|")
be = {}
for g in GATES:
    for nm in ("", ":nm"):
        right, wrong, nr, nw = 0.0, 0.0, 0, 0
        for h in CANCER:
            o = picks[h][f"oracle_wht:{g}{nm}"]
            for k in range(3):
                a = f"random:{g}:s{k}{nm}"
                p = picks[h][a]
                dd = pt[(h, f"{a}|pds")] - pt[(h, f"{BASE}|pds")]
                m = (p >= 0) & (o >= 0)
                r = m & (p == o); w = m & (p != o)
                right += float(dd[r].sum()); wrong += float(dd[w].sum()); nr += int(r.sum()); nw += int(w.sum())
        dr, dw = right / max(nr, 1), wrong / max(nw, 1)
        n_gate = float(np.mean([meta[h]["gates"][g] for h in CANCER]))
        n_t = float(np.mean([meta[h]["n_targets"] for h in CANCER]))
        # fold-level d pds of a selector with accuracy p on the gated targets: (n_gate / n_t) * (p * dr + (1 - p) * dw)
        pstar = -dw / (dr - dw) if dr != dw else float("nan")
        p5 = (0.005 * n_t / n_gate - dw) / (dr - dw) if dr != dw else float("nan")
        rnd = float(np.mean([compare(f"random:{g}:s{k}{nm}", CANCER)["d_pds"] for k in range(3)]))
        best_ag = max(cells[a]["oracle_agree"] for a in cells if gate_of(a) == g and is_nm(a) == bool(nm) and np.isfinite(cells[a]["oracle_agree"]))
        be[f"{g}{nm}"] = dict(n_right=nr, n_wrong=nw, d_right=dr, d_wrong=dw, random_mean=rnd, p_star=pstar, p_for_0005=p5, best_agreement=best_ag)
        Pr(f"| {g} | {'matched' if nm else 'as-is'} | {nr + nw} | {dr:+.4f} | {dw:+.4f} | {rnd:+.5f} | {100*pstar:.0f} % | {100*p5:.0f} % | {100*best_ag:.0f} % |")
out["break_even"] = be
Pr("")

(WORK / "analysis263.json").write_text(json.dumps(out, indent=1, default=float))
(WORK / "tables263.md").write_text("\n".join(md) + "\n")
print(f"wrote {WORK / 'analysis263.json'} and {WORK / 'tables263.md'}")
