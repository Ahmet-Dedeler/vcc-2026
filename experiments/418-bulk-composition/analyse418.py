"""418 analysis: PLOLO composition arms -> board predictions with 414's instruments, per fold.

References (all PP = 415's reconstruction of the shipped composition):
    PP|none  = 395's bulk (319p2, unkeyed)      PP|K3 = 410b's bulk      PP|K3+sh = 417a's bulk (q 1.25)
pds    board-emulated (414: the shipped artifact's non-scaling beta at lam 0.5) is primary; literal reported.
mse    414's cliff-aware route: 187's per-context identity MSE_c - 1 = [A_c + Kc_c - min(claim_c, B_c) - 2X]/sumD,
       A_c/B_c = 410b's MEASURED values x the PLOLO render-space ratio arm / PP|K3,
       X(arm) = X_410b (1 + k rel_arm)/(1 + k rel_PPK3) sqrt(mean A ratio), rel = cos_agg / cos_agg(PP|none) - 1,
       X_410b solved from 410b's board raw (0.730994) on its measured A_c/B_c. k = 1 and k = 1.407 (key class).
overall = (2.2138 dpds_raw - 1.051 dmse_raw) / 6.
"amp-matched" = the same arm with each context's bulk rescaled to the reference's render-space energy (direction only).
RETRO: historical one-variable composition steps with board pds/mse; X solved from each step's base board raw.

    uv run python analyse418.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
W414 = ROOT / "experiments/414-bulk-reprice/work"
for p in ("experiments/187-joint-amplitude", "experiments/096-mse-floor", "src"):
    sys.path.insert(0, str(ROOT / p))
import joint187  # noqa: E402

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
CTX = ("A", "B", "C")
PDS_SC, MSE_SC = 2.2138, -1.051
K_KEY = 1.407181520135006
KS = (1.0, K_KEY)
BOARD_410B_RAW = 0.730994
# board raw (pds_cosine, expr_mse_unbiased_capped_norm) from data/board-entries/*.json and data/leaderboard snapshots
BOARD = {
    "R:rc_callset(f0.10)": (0.8359866220735785, 0.8738045054832962),   # 31oB9yWy8CSNBDEkY4iW (lb 09-03T1815Z)
    "R:cd4f034": (0.8322779635823114, 0.8365398403169783),              # QSe3JlNUYZ2BzvlcsL14 (lb 09-04T0313Z)
    "R:h1union": (0.8357227796358231, 0.8172630319018631),              # iUziqE1gewfcPKzYUDaB
    "R:cd4svd5": (0.850642883686362, 0.821806190726257),                # Xv0ayvuoHpsKYWv2iu42
    "R:312(plain-1.5)": (0.8561538461538463, 0.8167140977906326),       # Q519I2pLWRCW7PAEx5ul
    "R:312b(plain-3)": (0.8570642883686362, 0.820313865913021),         # 1n5LFk6PSEHq15UHENIE
    "R:319(cap-1.5)": (0.8552471200297287, 0.8158761658635082),         # xuDelEik7zGWMXdNu0R4
    "R:319d3(cap-3)": (0.857699739873653, 0.815582524600044),           # 5LOjb2VLiIBh0JG6aAgg
    "R:254(svd5+p0.34)": (0.8506391675956894, 0.8173155951053511),      # PVHrJHweHD5gvOf43ftL
    "R:319p(p0.34)": (0.8551542177629133, 0.8118335038164238),          # SvzrTOQJAAaNjawzSyrC
    "R:319p2(p0.5)": (0.854738015607581, 0.8098975167495851),           # lqWNXX5YaYArdglu6Stn (= 395's bulk)
}
STEPS = [  # (base, arm, the one variable)
    ("R:rc_callset(f0.10)", "R:cd4f034", "CD4 target f 0.10 -> 0.34"),
    ("R:cd4f034", "R:h1union", "H1 union g 0 -> 0.28"),
    ("R:h1union", "R:cd4svd5", "CD4 deflation k 0 -> 5"),
    ("R:cd4svd5", "R:312(plain-1.5)", "normrel off -> plain -1.5"),
    ("R:312(plain-1.5)", "R:312b(plain-3)", "normrel plain -1.5 -> -3"),
    ("R:312(plain-1.5)", "R:319(cap-1.5)", "normrel plain -> cap (-1.5)"),
    ("R:319(cap-1.5)", "R:319d3(cap-3)", "normrel cap -1.5 -> -3"),
    ("R:cd4svd5", "R:254(svd5+p0.34)", "partial g 0 -> 0.34 (no normrel)"),
    ("R:319(cap-1.5)", "R:319p(p0.34)", "partial g 0 -> 0.34 (cap)"),
    ("R:319p(p0.34)", "R:319p2(p0.5)", "partial g 0.34 -> 0.5"),
]
MAIN_ARMS = ["a:cd4x0.5", "a:cd4x0.75", "a:cd4x1.25", "a:cd4x1.5", "a:cd4x2", "b:h1x0", "b:h1x0.5", "b:h1x1.5",
             "b:h1x2", "c:k0", "c:k3", "c:k8", "d:g0", "d:g0.25", "d:g0.75", "d:g1", "e:exp0", "e:exp-0.75",
             "e:exp-3", "f1:key-base-only", "f2:plain-1.5"]


def mse_c(A, B, Kc, claim, X, sumD):
    return 1.0 + (A + Kc - min(claim, B) - 2.0 * X) / sumD


class Cliff:
    def __init__(self):
        m = joint187.Model(joint187.F_FIT)
        assert abs(m.X - 0.017930409533420316) < 1e-12
        self.claim = {c: m.t[c]["claim"] for c in CTX}
        self.Kc = {c: m.t[c]["Kc"] for c in CTX}
        self.sumD = joint187.sumD0
        self.dec = json.loads((W414 / "decomp414.json").read_text())
        self.X410b = self.solve_X({c: self.dec["410b"][c] for c in CTX}, BOARD_410B_RAW)
        an = json.loads((W414 / "analyse414.json").read_text())
        assert abs(self.X410b - an["Xb"]["410b"]) < 1e-12, "X_410b does not reproduce analyse414"
        self.X395 = self.solve_X({c: self.dec["395"][c] for c in CTX}, 0.809898)
        assert abs(self.X395 - an["Xb"]["395"]) < 1e-12

    def solve_X(self, AB, raw):
        num = np.mean([AB[c]["A"] + self.Kc[c] - min(self.claim[c], AB[c]["B"]) for c in CTX])
        return float((num - (raw - 1.0) * self.sumD) / 2.0)

    def raw(self, AB, X, s=None):
        s = s or {c: 1.0 for c in CTX}
        return float(np.mean([mse_c(s[c] ** 2 * AB[c]["A"], s[c] ** 2 * AB[c]["B"], self.Kc[c], self.claim[c],
                                    s[c] * X, self.sumD) for c in CTX]))


def load(f):
    p = WORK / f"comp418-{f}.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    return d if "peak_rss_gib" in d else None      # only finished folds (the runner writes partial JSON as it goes)


ABC = None     # abc418.json arms, loaded in main(): the render-space energies on the REAL A/B/C panel with the shipped key


def abc_key(name):
    """PLOLO arm tag -> (abc418 arm, render variant)."""
    if name == "PP|none":
        return "PP", "none"
    arm, variant = name.split("|")
    return arm, variant.replace("V1", "K3")      # A/B/C carries the shipped key only; V1's energy ratio is K3's


def price(cl, R, arm, ref, k, matched=False, k_pds=1.0, energy="plolo"):
    """Arm vs ref on one fold.

    The REFERENCE sits at 414's calibration: A_c/B_c = 410b's measured x (PLOLO ref / PP|K3), and
    X_ref = X_410b (1 + k_key rel_ref)/(1 + k_key rel_PPK3) sqrt(mean A ratio) with rel vs PP|none (== X_410b for PP|K3).
    The ARM's own change is carried with k on ITS relative cosine change vs the reference (the retro's k_fit has the
    same definition):  X_arm = X_ref (1 + k r) sqrt(mean A_arm / mean A_ref),  r = cos_agg(arm)/cos_agg(ref) - 1.
    matched: each context's arm bulk rescaled to the reference's render-space energy (direction only).
    k_pds multiplies the board-emulated pds change (1 = literal transfer; the retro supplies class values).
    """
    base = R["PP|none"]["cos_agg"]
    anchor = R["PP|K3"]
    def AB(x):
        if energy == "abc":
            a, v = abc_key(x)
            e, e0 = ABC[a][v], ABC["PP"]["K3"]
            return {c: dict(A=cl.dec["410b"][c]["A"] * e[c]["A"] / e0[c]["A"],
                            B=cl.dec["410b"][c]["B"] * e[c]["B"] / e0[c]["B"]) for c in CTX}
        return {c: dict(A=cl.dec["410b"][c]["A"] * R[x]["AB"][c]["A"] / anchor["AB"][c]["A"],
                        B=cl.dec["410b"][c]["B"] * R[x]["AB"][c]["B"] / anchor["AB"][c]["B"]) for c in CTX}
    relK3 = anchor["cos_agg"] / base - 1
    A410 = np.mean([cl.dec["410b"][c]["A"] for c in CTX])
    ABa, ABr = AB(arm), AB(ref)
    rel_ref = R[ref]["cos_agg"] / base - 1
    Xr = cl.X410b * (1 + K_KEY * rel_ref) / (1 + K_KEY * relK3) * np.sqrt(np.mean([ABr[c]["A"] for c in CTX]) / A410)
    rr = cl.raw(ABr, Xr)
    r = R[arm]["cos_agg"] / R[ref]["cos_agg"] - 1
    amp = np.sqrt(np.mean([ABa[c]["A"] for c in CTX]) / np.mean([ABr[c]["A"] for c in CTX]))
    s = {c: float(np.sqrt(ABr[c]["A"] / ABa[c]["A"])) for c in CTX} if matched else None
    ra = cl.raw(ABa, Xr * (1 + k * r) * amp, s)       # raw() scales A,B by s^2 and X by s per context
    dpb = k_pds * (R[arm]["pds_lam0.5"] - R[ref]["pds_lam0.5"])
    dpl = R[arm]["pds"] - R[ref]["pds"]
    dm = ra - rr
    return dict(d_rel=r, dpds_board=dpb, dpds_lit=dpl, dmse=dm,
                ov=(PDS_SC * dpb + MSE_SC * dm) / 6, ov_lit=(PDS_SC * dpl + MSE_SC * dm) / 6,
                Ar={c: ABa[c]["A"] / ABr[c]["A"] for c in CTX},
                rho={c: ABa[c]["B"] / cl.claim[c] for c in CTX})


def boot(f, arm, ref, n=2000, seed=418, src="arrays418"):
    Z = np.load(WORK / f"{src}-{f}.npz")
    d = Z["d"]
    rng = np.random.default_rng(seed)
    xa, aa, pa = Z[f"{arm}__x"], Z[f"{arm}__a"], Z[f"{arm}__pdsl"]
    xr, ar, pr = Z[f"{ref}__x"], Z[f"{ref}__a"], Z[f"{ref}__pdsl"]
    agg = lambda x, a, i: np.mean([abs(x[c, i].sum()) / np.sqrt(a[c, i].sum() * d[i].sum()) for c in range(3)])
    rel, dp = [], []
    for _ in range(n):
        i = rng.integers(0, d.size, d.size)
        rel.append(agg(xa, aa, i) / agg(xr, ar, i) - 1)
        dp.append(pa[:, i].mean() - pr[:, i].mean())
    return dict(p_rel=float(np.mean(np.asarray(rel) > 0)), p_pds=float(np.mean(np.asarray(dp) > 0)),
                pds_ci=np.quantile(dp, [0.025, 0.975]).tolist())


FAMILY_STEP = {  # arm family -> the historical board step(s) of the same class that calibrate its transfer
    "a": ["R:rc_callset(f0.10) -> R:cd4f034"],
    "b": ["R:cd4f034 -> R:h1union"],
    "c": ["R:h1union -> R:cd4svd5"],
    "d": ["R:319(cap-1.5) -> R:319p(p0.34)", "R:319p(p0.34) -> R:319p2(p0.5)", "R:cd4svd5 -> R:254(svd5+p0.34)"],
}


def retro(cl, RS, folds, P, energy="plolo"):
    """Historical steps: board deltas vs PLOLO. Per fold the step's base A_c/B_c = 395's measured x PLOLO ratio to
    PP|none (= 395's bulk), X_base solved from the base's own board raw; the arm's X = X_base (1 + k r) sqrt(A ratio)."""
    out = {}
    P(f"{'step':34s} {'board dpds':>10s} {'PLOLO dpds':>10s} {'k_pds':>6s} | {'board dmse':>10s} {'pred k=1':>9s} "
      f"{'pred kK':>9s} {'k_fit':>6s} | {'rel cos':>7s} {'A ratio':>7s} | board ov(pds+mse) vs pred k=1 | fold pds signs")
    for b, a, what in STEPS:
        per = []
        for f in folds:
            R = RS[f]
            E = ({x: ABC[x]["none"] for x in (a, b)} | {"PP|none": ABC["PP"]["none"]} if energy == "abc"
                 else {x: R[x]["AB"] for x in (a, b, "PP|none")})
            ABb = {c: dict(A=cl.dec["395"][c]["A"] * E[b][c]["A"] / E["PP|none"][c]["A"],
                           B=cl.dec["395"][c]["B"] * E[b][c]["B"] / E["PP|none"][c]["B"]) for c in CTX}
            ABa = {c: dict(A=ABb[c]["A"] * E[a][c]["A"] / E[b][c]["A"],
                           B=ABb[c]["B"] * E[a][c]["B"] / E[b][c]["B"]) for c in CTX}
            Xb = cl.solve_X(ABb, BOARD[b][1])
            r = R[a]["cos_agg"] / R[b]["cos_agg"] - 1
            amp = np.sqrt(np.mean([ABa[c]["A"] for c in CTX]) / np.mean([ABb[c]["A"] for c in CTX]))
            pred = {k: cl.raw(ABa, Xb * amp * (1 + k * r)) - BOARD[b][1] for k in (0.0, 1.0, K_KEY)}
            per.append(dict(dpds=R[a]["pds"] - R[b]["pds"], pred=pred, rel=r,
                            Ar=float(np.mean([ABa[c]["A"] / ABb[c]["A"] for c in CTX]))))
        bd_p, bd_m = BOARD[a][0] - BOARD[b][0], BOARD[a][1] - BOARD[b][1]
        m = lambda k: float(np.mean([x[k] for x in per]))
        p0 = float(np.mean([x["pred"][0.0] for x in per]))
        pk1 = float(np.mean([x["pred"][1.0] for x in per]))
        pkk = float(np.mean([x["pred"][K_KEY] for x in per]))
        kfit = (bd_m - p0) / (pk1 - p0) if abs(pk1 - p0) > 1e-12 else float("nan")   # fold-mean prediction is linear in k
        kpds = bd_p / m("dpds") if abs(m("dpds")) > 1e-9 else float("nan")
        out[f"{b} -> {a}"] = dict(what=what, board_dpds=bd_p, board_dmse=bd_m, plolo_dpds=m("dpds"), k_pds=kpds,
                                  pred_dmse_k0=p0, pred_dmse_k1=pk1, pred_dmse_kkey=pkk, k_fit=kfit, rel=m("rel"),
                                  Ar=m("Ar"), plolo_dpds_folds=[x["dpds"] for x in per],
                                  pred_dmse_k1_folds=[x["pred"][1.0] for x in per])
        P(f"{what:34s} {bd_p:+10.5f} {m('dpds'):+10.5f} {kpds:6.2f} | {bd_m:+10.5f} {pk1:+9.5f} {pkk:+9.5f} {kfit:6.2f} | "
          f"{m('rel'):+7.4f} {m('Ar'):7.4f} | {(PDS_SC * bd_p + MSE_SC * bd_m) / 6:+.5f} vs {(PDS_SC * m('dpds') + MSE_SC * pk1) / 6:+.5f}"
          f" | {''.join('+' if x['dpds'] > 0 else '-' for x in per)}")
    bo = np.array([(PDS_SC * v["board_dpds"] + MSE_SC * v["board_dmse"]) / 6 for v in out.values()])
    pr = np.array([(PDS_SC * v["plolo_dpds"] + MSE_SC * v["pred_dmse_k1"]) / 6 for v in out.values()])
    bp = np.array([v["board_dpds"] for v in out.values()]); pp = np.array([v["plolo_dpds"] for v in out.values()])
    bm = np.array([v["board_dmse"] for v in out.values()]); pm = np.array([v["pred_dmse_k1"] for v in out.values()])
    sgn = lambda x, y: int(np.sum(np.sign(x) == np.sign(y)))
    summ = dict(sign_overall=sgn(bo, pr), sign_pds=sgn(bp, pp), sign_mse=sgn(bm, pm), n=len(bo),
                corr_overall=float(np.corrcoef(bo, pr)[0, 1]), corr_pds=float(np.corrcoef(bp, pp)[0, 1]),
                corr_mse=float(np.corrcoef(bm, pm)[0, 1]))
    P(f"retro sign agreement: overall(pds+mse) {summ['sign_overall']}/{len(bo)}, pds {summ['sign_pds']}/{len(bp)}, "
      f"mse {summ['sign_mse']}/{len(bm)}; corr overall {summ['corr_overall']:+.3f}, pds {summ['corr_pds']:+.3f}, "
      f"mse {summ['corr_mse']:+.3f}")
    return out, summ


def main():
    global ABC
    cl = Cliff()
    ABC = json.loads((WORK / "abc418.json").read_text())["arms"]
    folds = [f for f in (*CANCER, "kolf2") if load(f) is not None]
    cancer = [f for f in CANCER if f in folds]
    RS = {f: load(f)["arms"] for f in folds}
    for f in folds:          # the dose continuation, where it exists (its PP rows are asserted identical at run time)
        px = WORK / f"comp418x-{f}.json"
        if px.exists() and "peak_rss_gib" in (dx := json.loads(px.read_text())):
            RS[f].update({k: v for k, v in dx["arms"].items() if k not in RS[f]})
    arms_all = MAIN_ARMS + [a for a in ("b:h1x2.5", "b:h1x3", "b:h1x4", "a:cd4x1.75")
                            if all(f"{a}|K3" in RS[f] for f in folds)]
    COMP = {f: load(f)["comp"] for f in folds}
    L = []
    P = lambda *a: (print(*a), L.append(" ".join(str(x) for x in a)))
    out = dict(folds=folds, X410b=cl.X410b, rows={})
    P(f"folds: {folds}; X_410b {cl.X410b:.6f} (analyse414 reproduced); k in {KS}; arms {len(arms_all)}")
    for f in folds:
        R = RS[f]
        sh = price(cl, R, "PP|K3+sh", "PP|K3", K_KEY)
        v1 = price(cl, R, "PP|V1", "PP|K3", K_KEY)
        P(f"[{f}] 410b-analog PP|K3 vs PP|none: rel {R['PP|K3']['cos_agg'] / R['PP|none']['cos_agg'] - 1:+.4f}, "
          f"dpds_bd {R['PP|K3']['pds_lam0.5'] - R['PP|none']['pds_lam0.5']:+.5f} (lit "
          f"{R['PP|K3']['pds'] - R['PP|none']['pds']:+.5f}) | 417a sharpening vs PP|K3: rel {sh['d_rel']:+.4f} dpds_bd "
          f"{sh['dpds_board']:+.5f} dmse {sh['dmse']:+.5f} ov {sh['ov']:+.5f} | 415 V1 vs K3: rel {v1['d_rel']:+.4f} "
          f"ov {v1['ov']:+.5f} | CTRL key shuffled rel {R['CTRL:PP|K3shuf']['cos_agg'] / R['PP|K3']['cos_agg'] - 1:+.4f}")

    # ---- RETRO first: it validates each energy route and calibrates the class constants of the 'hist' columns
    hist = {}
    for energy in ("plolo", "abc"):
        P(f"\n==== RETRO (cancer folds; render energies from {'PLOLO folds' if energy == 'plolo' else 'the REAL A/B/C panel'}): "
          "historical one-variable composition steps, no key. Board raw deltas vs PLOLO")
        out[f"retro_{energy}"], out[f"retro_{energy}_summary"] = retro(cl, RS, cancer, P, energy)
        h = {}
        for fam, steps in FAMILY_STEP.items():
            v = [out[f"retro_{energy}"][st] for st in steps]
            kp = float(np.sum([x["board_dpds"] for x in v]) / np.sum([x["plolo_dpds"] for x in v]))
            km = float(np.sum([x["board_dmse"] - x["pred_dmse_k0"] for x in v]) /
                       np.sum([x["pred_dmse_k1"] - x["pred_dmse_k0"] for x in v]))
            h[fam] = dict(k_pds=1.0 if fam == "d" else kp, k_mse=km, steps=steps)
        hist[energy] = h
        P("class transfer constants (pooled over each family's steps; d's k_pds fixed at 1, its pds moves < 0.0005): " +
          "; ".join(f"{k}: k_pds {v['k_pds']:.2f}, k_mse {v['k_mse']:.2f}" for k, v in h.items()))
    if "kolf2" in folds:
        P("\n==== RETRO (kolf2 only, A/B/C energies; 3 production donors + production operator)")
        out["retro_kolf2"], out["retro_kolf2_summary"] = retro(cl, RS, ["kolf2"], P, "abc")
    out["hist_constants"] = hist

    # ---- main table
    def tag_ref(arm, variant):
        if arm.startswith("f1:"):
            return {"K3": ("f1:key-base-only|K3", "PP|K3"), "K3+sh": ("f1:key-base-only|K3+sh", "PP|K3+sh")}.get(variant)
        return f"{arm}|{variant}", {"K3": "PP|K3", "K3+sh": "PP|K3+sh", "V1": "PP|V1", "V1+sh": "PP|V1+sh"}[variant]

    kk = f"k{K_KEY:.3g}"
    for arm in arms_all:
        row = {}
        fam = arm.split(":")[0]
        for variant in ("K3", "K3+sh", "V1", "V1+sh"):
            tr = tag_ref(arm, variant)
            if tr is None:
                continue
            a, r = tr
            per = {}
            for f in folds:
                e = {}
                for energy in ("plolo", "abc"):
                    pre_ = "" if energy == "plolo" else "abc_"
                    for k in KS:
                        e[f"{pre_}k{k:.3g}"] = price(cl, RS[f], a, r, k, energy=energy)
                    e[f"{pre_}matched_k1"] = price(cl, RS[f], a, r, 1.0, matched=True, energy=energy)
                    if fam in hist[energy]:
                        e[f"{pre_}hist"] = price(cl, RS[f], a, r, hist[energy][fam]["k_mse"],
                                                 k_pds=hist[energy][fam]["k_pds"], energy=energy)
                per[f] = e
            row[variant] = per
        out["rows"][arm] = row

    def summ(arm, variant, key, field="ov", fl=cancer):
        v = [out["rows"][arm][variant][f][key][field] for f in fl if key in out["rows"][arm][variant][f]]
        return (float(np.mean(v)), float(min(v)), float(max(v)), sum(x > 0 for x in v), len(v)) if v else None

    fmt = lambda s_: f"{s_[0]:+.5f} [{s_[1]:+.4f},{s_[2]:+.4f}] {s_[3]}/{s_[4]}" if s_ else f"{'-':>27s}"
    for variant, label in (("K3", "vs 410b (PP|K3)"), ("K3+sh", "vs 417a (PP|K3+sh)"),
                           ("V1", "secondary: the arm's own LODO-V1 key, vs PP|V1"),
                           ("V1+sh", "secondary: own LODO-V1 key + sharpening, vs PP|V1+sh")):
        P(f"\n==== {label}: Δoverall, cancer-fold mean [min, max] +folds/n; kolf2 separate (no vote)")
        P(f"{'arm':18s} {'d rel':>7s} {'dpds_bd':>8s} | {'PLOLO-energy k=1':>27s} | {'A/B/C-energy k=1':>27s} | "
          f"{'A/B/C-energy k=1.407':>27s} | {'hist':>8s} {'match':>8s} | {'A/B/C A ratio':>20s} | kolf2 abc k1/kK")
        rows = []
        for arm in arms_all:
            if variant not in out["rows"][arm]:
                continue
            m = lambda fld, key="abc_k1": float(np.mean([out["rows"][arm][variant][f][key][fld] for f in cancer]))
            ar = [float(np.mean([out["rows"][arm][variant][f]["abc_k1"]["Ar"][c] for f in cancer])) for c in CTX]
            kf = out["rows"][arm][variant].get("kolf2")
            sh_ = summ(arm, variant, "abc_hist")
            rows.append((summ(arm, variant, "abc_" + kk)[0],
                         f"{arm:18s} {m('d_rel'):+7.4f} {m('dpds_board'):+8.5f} | {fmt(summ(arm, variant, 'k1'))} | "
                         f"{fmt(summ(arm, variant, 'abc_k1'))} | {fmt(summ(arm, variant, 'abc_' + kk))} | "
                         f"{(f'{sh_[0]:+.5f}' if sh_ else '-'):>8s} {summ(arm, variant, 'abc_matched_k1')[0]:+.5f} | "
                         f"{ar[0]:.4f}/{ar[1]:.4f}/{ar[2]:.4f} | "
                         + (f"{kf['abc_k1']['ov']:+.5f}/{kf['abc_' + kk]['ov']:+.5f}" if kf else "-")))
        for _, line in sorted(rows, key=lambda t: -t[0]):
            P(line)

    # ---- bootstrap for the top and bottom arms (vs 417a)
    ranked = sorted([a for a in arms_all if "K3+sh" in out["rows"][a]], key=lambda a: -summ(a, "K3+sh", "abc_" + kk)[0])
    out["ranked_vs_417a"] = ranked
    P("\npaired target bootstrap (2000) vs 417a, per fold: P(d rel>0)/P(d pds_bd>0) [d pds_bd 95% CI]")
    out["boot"] = {}
    for a in ranked[:8] + ranked[-2:]:
        an, r = tag_ref(a, "K3+sh")
        bb = {}
        for f in folds:
            src = "arrays418x" if (WORK / f"arrays418x-{f}.npz").exists() and a in ("b:h1x2.5", "b:h1x3", "b:h1x4", "a:cd4x1.75") else "arrays418"
            bb[f] = boot(f, an, r, src=src)
        out["boot"][a] = bb
        P(f"  {a:18s} " + "  ".join(f"{f[:6]} {b['p_rel']:.2f}/{b['p_pds']:.2f} [{b['pds_ci'][0]:+.4f},{b['pds_ci'][1]:+.4f}]"
                                    for f, b in bb.items()))

    # ---- what each knob does to the effect-space row on the REAL panel (direction vs scale)
    P("\nper knob on the REAL A/B/C panel (abc418): solved fractions cd4 / h1 / partial ; post-add s ; rows whose direction "
      "moved ; median cos on those rows ; energy-weighted cos ; Frobenius ratio (pre-key) ; render A ratio vs PP|K3 (A/B/C)")
    for arm in ["PP"] + [a for a in arms_all if not a.startswith("f1")]:
        e = ABC[arm]
        ar = [e["K3"][c]["A"] / ABC["PP"]["K3"][c]["A"] for c in CTX]
        sh = e["shape"]
        P(f"  {arm:18s} {e['fractions']['cd4']:.4f} / {e['fractions']['h1_union']:.4f} / {e['fractions']['cd4_partial']:.4f} ; "
          f"s {e['s']:.4f} ; moved {sh['frac_rows_direction_moved']:.2f} ; cos {sh['median_cos_moved_rows']:.5f} ; "
          f"ew-cos {sh['energy_weighted_cos']:.5f} ; frob {sh['frob_ratio']:.4f} ; A {ar[0]:.4f}/{ar[1]:.4f}/{ar[2]:.4f}")
    (WORK / "analyse418.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "analyse418.log").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
