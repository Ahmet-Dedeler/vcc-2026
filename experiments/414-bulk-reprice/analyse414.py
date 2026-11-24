"""414 analysis: PLOLO arm tables -> board predictions.

pds    two reads: PLOLO literal (a23 read-back) and BOARD-EMULATED (board414/extra414: the shipped artifact's non-scaling
       shared pattern beta added at lam = 0.5, the lam at which PLOLO retrodicts 410b's board pds -0.0036 for the key).
mse    ROUTE 1 (the mandate's): 197's slope, raw = 1.0236 - 2.0745 cos, cos(395) = 0.107, relative PLOLO gains carried
       as 0.107 (1 + k rel); k = 1 literal, k_key = board relative gain of 410b / PLOLO rel of K3.
       ROUTE 2 (cliff-aware): 187's per-context model MSE_c - 1 = [A_c + Kc_c - min(claim_c, B_c) - 2X]/sumD with
       claim_c/Kc_c/sumD from joint187.Model(F_FIT) (KA: X23 0.017930, cliffs 0.9823/1.0129/1.0386), A_c/B_c MEASURED on
       the 395/410/410b artifacts (decomp414.json), X solved per artifact from its board raw; an arm's A_c, B_c are
       410b's times the PLOLO render-space energy ratio vs K3, X(arm) = X_410b (1 + k rel_arm)/(1 + k rel_K3) sqrt(A ratio).
       "lev" = the same arm with each context's bulk amplitude re-levelled onto its cliff s_c = sqrt(claim_c / B_c).
overall = (2.2138 dpds_raw - 1.051 dmse_raw) / 6.

    uv run python analyse414.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/187-joint-amplitude"))
sys.path.insert(0, str(ROOT / "experiments/096-mse-floor"))
sys.path.insert(0, str(ROOT / "src"))
import joint187  # noqa: E402

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
CTX = ("A", "B", "C")
SLOPE197, COS395 = 2.0745, 0.107
PDS_SC, MSE_SC = 2.2138, -1.051
BOARD_RAW = {"395": 0.809898, "410": 0.730213, "410b": 0.730994}
BOARD_PDS = {"395": 0.854738, "410": 0.850914, "410b": 0.851156}
SD187 = {"A": 0.0025, "B": 0.0004, "C": 0.0013}    # 187's retrodiction-region cliff sd in ln(alpha)
rng = np.random.default_rng(414)


def fold_rows(f):
    """One dict per arm for fold f, merging arms414 / board414 / extra414 / ballast414."""
    A = json.loads((WORK / f"arms414-{f}.json").read_text())["arms"]
    Bd = json.loads((WORK / f"board414-{f}.json").read_text())["arms"] if (WORK / f"board414-{f}.json").exists() else {}
    X = json.loads((WORK / f"extra414-{f}.json").read_text())["arms"] if (WORK / f"extra414-{f}.json").exists() else {}
    Bl = json.loads((WORK / f"ballast414-{f}.json").read_text())["arms"] if (WORK / f"ballast414-{f}.json").exists() else {}
    rows = {}
    for name, v in A.items():
        rows[name] = dict(cos_agg=v["cos_agg"], pds=v["pds"], AB=v["AB"], note=v["note"])
        if name in Bd:
            rows[name].update(pds_l05=Bd[name]["lam0.5"]["pds"], pds_l1=Bd[name]["lam1"]["pds"])
    for name, v in X.items():
        if name in rows:
            continue
        AB = {c: dict(A=v["AB"][c]["A"], B=v["AB"][c].get("B", np.nan)) for c in CTX}
        rows[name] = dict(cos_agg=v["cos_agg"], pds=v["pds"], AB=AB, note="extra414",
                          pds_l05=v["pds_lam0.5"], pds_l1=v["pds_lam1"])
    for name, v in Bl.items():
        if name in rows:
            rows[name].update(pds_real=v["pds_realised"], cos_real=v["cos_agg_realised"], r_min_q=v["r_min_q"],
                              frac_K_capped=v["frac_K_capped"], clip_mass_q=v["clip_mass_q"])
    return rows


def per_arm(R, folds):
    out = {}
    names = [n for n in R[folds[0]] if all(n in R[f] for f in folds)]
    for a in names:
        fr = []
        for f in folds:
            b, v = R[f]["base"], R[f][a]
            Br = {c: (v["AB"][c]["B"] / b["AB"][c]["B"]) if np.isfinite(v["AB"][c]["B"]) else
                  v["AB"][c]["A"] / b["AB"][c]["A"] for c in CTX}
            fr.append(dict(rel=v["cos_agg"] / b["cos_agg"] - 1, dpds=v["pds"] - b["pds"],
                           dpds_l05=(v["pds_l05"] - b["pds_l05"]) if "pds_l05" in v else np.nan,
                           dpds_l1=(v["pds_l1"] - b["pds_l1"]) if "pds_l1" in v else np.nan,
                           dpds_real=(v["pds_real"] - b["pds_real"]) if "pds_real" in v else np.nan,
                           rel_real=(v["cos_real"] / b["cos_real"] - 1) if "cos_real" in v else np.nan,
                           Ar={c: v["AB"][c]["A"] / b["AB"][c]["A"] for c in CTX}, Br=Br,
                           kcap=v.get("frac_K_capped", np.nan)))
        m = lambda k: float(np.mean([x[k] for x in fr]))
        out[a] = dict(rel=m("rel"), dpds=m("dpds"), dpds_l05=m("dpds_l05"), dpds_l1=m("dpds_l1"),
                      dpds_real=m("dpds_real"), rel_real=m("rel_real"), kcap=m("kcap"),
                      rel_f=[x["rel"] for x in fr], dpds_f=[x["dpds"] for x in fr], dpds_l05_f=[x["dpds_l05"] for x in fr],
                      Ar={c: float(np.mean([x["Ar"][c] for x in fr])) for c in CTX},
                      Br={c: float(np.mean([x["Br"][c] for x in fr])) for c in CTX}, note=R[folds[0]][a]["note"])
    return out


def boot(f, arm, ref, n=2000):
    Z = np.load(WORK / f"arrays414-{f}.npz")
    if f"{arm}__x" not in Z.files:
        return None
    d = Z["d"]
    xa, aa, pa = Z[f"{arm}__x"], Z[f"{arm}__a"], Z[f"{arm}__pds"]
    xr, ar, pr = Z[f"{ref}__x"], Z[f"{ref}__a"], Z[f"{ref}__pds"]
    agg = lambda x, a, i: np.mean([abs(x[c, i].sum()) / np.sqrt(a[c, i].sum() * d[i].sum()) for c in range(3)])
    rel, dp = [], []
    for _ in range(n):
        i = rng.integers(0, d.size, d.size)
        rel.append(agg(xa, aa, i) / agg(xr, ar, i) - 1)
        dp.append(pa[:, i].mean() - pr[:, i].mean())
    return dict(rel_ci=np.quantile(rel, [0.025, 0.975]).tolist(), rel_p=float(np.mean(np.asarray(rel) > 0)),
                dpds_ci=np.quantile(dp, [0.025, 0.975]).tolist(), dpds_p=float(np.mean(np.asarray(dp) > 0)))


def mse_c(A, B, Kc, claim, X, sumD):
    return 1.0 + (A + Kc - min(claim, B) - 2.0 * X) / sumD


def raw_of(AB, Kc, claim, X, sumD, s=None):
    s = s or {c: 1.0 for c in CTX}
    per = {c: mse_c(s[c] ** 2 * AB[c]["A"], s[c] ** 2 * AB[c]["B"], Kc[c], claim[c], s[c] * X, sumD) for c in CTX}
    return float(np.mean(list(per.values()))), per


def main() -> None:
    R = {f: fold_rows(f) for f in (*CANCER, "kolf2")}
    T = per_arm(R, list(CANCER))
    Kf = per_arm(R, ["kolf2"])
    dec = json.loads((WORK / "decomp414.json").read_text())
    m = joint187.Model(joint187.F_FIT)
    assert abs(m.X - 0.017930409533420316) < 1e-12
    assert all(abs(m.cliff[c] - v) < 1e-9 for c, v in zip(CTX, (0.9823377082708273, 1.0129099639546535, 1.0386484485601954)))
    claim = {c: m.t[c]["claim"] for c in CTX}
    Kc = {c: m.t[c]["Kc"] for c in CTX}
    sumD = joint187.sumD0
    L = []
    P = lambda *a: (print(*a), L.append(" ".join(str(x) for x in a)))
    out = dict(claim=claim, Kc=Kc, sumD=sumD)

    # ------------------------------------------------------------------ known answers
    d393 = json.loads((ROOT / "experiments/393-mse-shared-component/work/decomp393.json").read_text())["per_context"]
    ka393 = max(abs(dec["395"][c][k] - d393[c][k]) for c in CTX for k in ("A", "B"))
    P(f"KA-187: X23 {m.X:.6f}, cliffs {[round(m.cliff[c], 4) for c in CTX]} (187: 0.017930; 0.9823/1.0129/1.0386) PASS")
    P(f"KA-decomp: 395 artifact A_c/B_c == 393's decomposition (the 319p2 bulk), worst |d| {ka393:.1e}")
    ka410 = {f: json.loads((WORK / f"arms414-{f}.json").read_text())["ka"] for f in (*CANCER, "kolf2")}
    P("KA-410 (check410 gains reproduced to 0):", {f: {k: round(v, 3) for k, v in ka410[f]["check410"].items()}
                                                   for f in ka410})
    P("KA-228 (arms228 rows t = 0/0.1/0.25/0.5/1 reproduced to 1e-12): PASS on", list(ka410),
      "| KA-169 (plolo169 q 1.00/1.15):", {f: ka410[f]["plolo169_worst"] for f in ka410})
    out["ka"] = dict(decomp_vs_393=ka393, check410={f: ka410[f]["check410"] for f in ka410})

    # ------------------------------------------------------------------ board anchors
    def solve_X(AB, raw):
        num = np.mean([AB[c]["A"] + Kc[c] - min(claim[c], AB[c]["B"]) for c in CTX])
        return float((num - (raw - 1.0) * sumD) / 2.0)
    Xb = {k: solve_X(dec[k], BOARD_RAW[k]) for k in ("395", "410", "410b")}
    Abar = {k: float(np.mean([dec[k][c]["A"] for c in CTX])) for k in Xb}
    P("\nboard X (187 one-equation solve on measured A_c/B_c):", {k: round(v, 6) for k, v in Xb.items()})
    for k in ("395", "410", "410b"):
        P(f"  {k:5s} A_c {[round(dec[k][c]['A'], 6) for c in CTX]}  rho_c = B_c/claim_c "
          f"{[round(dec[k][c]['B'] / claim[c], 5) for c in CTX]}  cliff s_c* {[round(float(np.sqrt(claim[c] / dec[k][c]['B'])), 5) for c in CTX]}")
    bcr = {k: Xb[k] / Xb["395"] / np.sqrt(Abar[k] / Abar["395"]) - 1 for k in ("410", "410b")}
    P("  board relative aggregate-cosine gain vs 395:", {k: round(v, 4) for k, v in bcr.items()},
      "(197 slope reads", {k: round((BOARD_RAW["395"] - BOARD_RAW[k]) / SLOPE197 / COS395, 4) for k in ("410", "410b")}, ")")
    rel3 = T["K3"]["rel"]
    k_key = bcr["410b"] / rel3
    P(f"  calibration k_key = {bcr['410b']:.4f} / PLOLO K3 {rel3:.4f} = {k_key:.3f}")
    P("\nrender-space energy of the key: PLOLO A ratio (cancer mean) vs artifacts, vs 395 / base")
    for c in CTX:
        P(f"  {c}: PLOLO K3 {T['K3']['Ar'][c]:.4f}  Kinf {T['Kinf']['Ar'][c]:.4f} | artifact 410 "
          f"{dec['410'][c]['A'] / dec['395'][c]['A']:.4f}  410b {dec['410b'][c]['A'] / dec['395'][c]['A']:.4f}")
    # retrodiction of the 410 (welded) board raw with k_key calibrated on 410b
    AB410 = {c: dict(A=dec["395"][c]["A"] * T["K3"]["Ar"][c], B=dec["395"][c]["B"] * T["K3"]["Br"][c]) for c in CTX}
    for kk, kv in (("k=1", 1.0), ("k_key", k_key)):
        X = Xb["395"] * (1 + kv * rel3) * np.sqrt(np.mean([AB410[c]["A"] for c in CTX]) / Abar["395"])
        r, _ = raw_of(AB410, Kc, claim, X, sumD)
        P(f"  retrodict 410/410b board raw from 395 + PLOLO K3 ({kk}): {r:.6f}  (board 410 {BOARD_RAW['410']:.6f}, "
          f"410b {BOARD_RAW['410b']:.6f})")
    P(f"  pds, K3 vs base (cancer mean): PLOLO literal {T['K3']['dpds']:+.5f}; board-emulated lam 0.5 "
      f"{T['K3']['dpds_l05']:+.5f}, lam 1 {T['K3']['dpds_l1']:+.5f}; BOARD 410b - 395 {BOARD_PDS['410b'] - BOARD_PDS['395']:+.6f}")
    out.update(Xb=Xb, board_cos_gain=bcr, k_key=k_key)

    # ------------------------------------------------------------------ arm predictions
    def predict(a, ref, k):
        v, r = T[a], T[ref]
        d197 = -SLOPE197 * COS395 * k * (v["rel"] - r["rel"])
        def AB_of(w):
            return {c: dict(A=dec["410b"][c]["A"] * w["Ar"][c] / T["K3"]["Ar"][c],
                            B=dec["410b"][c]["B"] * w["Br"][c] / T["K3"]["Br"][c]) for c in CTX}
        def X_of(w, ABw):
            return Xb["410b"] * (1 + k * w["rel"]) / (1 + k * T["K3"]["rel"]) * np.sqrt(
                np.mean([ABw[c]["A"] for c in CTX]) / Abar["410b"])
        ABa, ABr = AB_of(v), AB_of(r)
        ra, per = raw_of(ABa, Kc, claim, X_of(v, ABa), sumD)
        rr, _ = raw_of(ABr, Kc, claim, X_of(r, ABr), sumD)
        s_lev = {c: float(np.sqrt(claim[c] / ABa[c]["B"])) for c in CTX}
        rl, _ = raw_of(ABa, Kc, claim, X_of(v, ABa), sumD, s_lev)
        ov = lambda dm, dp: (PDS_SC * dp + MSE_SC * dm) / 6.0
        dp0 = v["dpds"] - r["dpds"]
        dpb = v["dpds_l05"] - r["dpds_l05"]
        return dict(rel=v["rel"], d_rel=v["rel"] - r["rel"], dpds=dp0, dpds_board=dpb,
                    dpds_l1=v["dpds_l1"] - r["dpds_l1"], d197=d197, dcliff=ra - rr, dlev=rl - rr, s_lev=s_lev,
                    rho={c: ABa[c]["B"] / claim[c] for c in CTX},
                    ov_mandate=ov(d197, dp0), ov_primary=ov(ra - rr, dpb if np.isfinite(dpb) else dp0),
                    ov_mse_only=ov(ra - rr, 0.0), ov_lev=ov(rl - rr, dpb if np.isfinite(dpb) else dp0),
                    Ar=v["Ar"], rel_f=v["rel_f"], dpds_f=[x - y for x, y in zip(v["dpds_f"], r["dpds_f"])],
                    dpds_board_f=[x - y for x, y in zip(v["dpds_l05_f"], r["dpds_l05_f"])],
                    kcap=v["kcap"], dpds_real=v["dpds_real"] - r["dpds_real"], rel_real=v["rel_real"])

    table = {a: dict(vs410b=predict(a, "K3", k_key), vs410b_k1=predict(a, "K3", 1.0), vs410c=predict(a, "Kinf", k_key))
             for a in T}
    out["table"], out["plolo_cancer"], out["plolo_kolf2"] = table, T, Kf

    hdr = (f"{'arm':22s} {'rel':>6s} {'d rel':>6s} {'dpds':>8s} {'dpds_bd':>8s} | {'d197':>7s} {'dcliff':>7s} {'dlev':>7s} | "
           f"{'ov_mand':>8s} {'ov_prim':>8s} {'ov_mse':>8s} {'ov_lev':>8s} | A-ratio  | folds rel>0 / pds_bd>0")
    for refname, key in (("410b (K3)", "vs410b"), ("410c (Kinf)", "vs410c")):
        P(f"\n==== vs {refname}; k_key {k_key:.3f}; dpds_bd = board-emulated (lam 0.5); ov_mand = 197 slope k_key + literal pds")
        P(hdr)
        for a, row in sorted(table.items(), key=lambda kv: -kv[1][key]["ov_primary"]):
            r = row[key]
            nrel = sum(x > 0 for x in [p - q for p, q in zip(T[a]["rel_f"], T["Kinf" if key == "vs410c" else "K3"]["rel_f"])])
            npd = sum(x > 0 for x in r["dpds_board_f"] if np.isfinite(x))
            P(f"{a:22s} {r['rel']:+6.3f} {r['d_rel']:+6.3f} {r['dpds']:+8.5f} {r['dpds_board']:+8.5f} | {r['d197']:+7.4f} "
              f"{r['dcliff']:+7.4f} {r['dlev']:+7.4f} | {r['ov_mandate']:+8.5f} {r['ov_primary']:+8.5f} {r['ov_mse_only']:+8.5f} "
              f"{r['ov_lev']:+8.5f} | {np.mean(list(r['Ar'].values())):.3f} | {nrel}/3  {npd}/3")

    P("\nkolf2 (reported apart; production operator, 3 production donors, non-cancer truth), vs Kinf:")
    for a in ("K3", "K0.5_inf", "K0.7_inf", "Kinf+gf@0.25", "gi@0.25+Kinf", "Kinf+sh@1.05", "Kinf+sh@1.1", "Kinf+sh@1.15",
              "Kinf+sh@1.25", "Kinf+sh@1.35", "Kinf+sh@1.5", "Kinf+shL@1.15", "Kinf+sh@0.87", "Kinf+shperm@1.15",
              "K0.5_inf+sh@1.25"):
        if a in Kf:
            v, r = Kf[a], Kf["Kinf"]
            P(f"  {a:20s} d rel {v['rel'] - r['rel']:+.4f} (rel/Kinf {(1 + v['rel']) / (1 + r['rel']) - 1:+.4f})  "
              f"dpds {v['dpds'] - r['dpds']:+.5f}  dpds lam0.5 {v['dpds_l05'] - r['dpds_l05']:+.5f}  "
              f"A ratio {np.mean(list(v['Ar'].values())) / np.mean(list(r['Ar'].values())):.4f}")

    P("\nballast realisability (builder ballast_count emulation, K <= 16 x 995k), cancer mean:")
    for a in ("K3", "Kinf", "K0.5_inf", "Kinf+sh@1.15", "Kinf+sh@1.25", "Kinf+sh@1.5", "K0.5_inf+sh@1.25", "K3+sh@1.15"):
        if a in T and np.isfinite(T[a]["kcap"]):
            rq = [R[f][a]["r_min_q"] for f in CANCER]
            P(f"  {a:18s} K capped {100 * T[a]['kcap']:.1f}% of blocks; r_min median {np.mean([x[2] for x in rq]):.3f}; "
              f"realised-vs-design: d rel {T[a]['rel_real'] - T[a]['rel']:+.4f}, d pds "
              f"{T[a]['dpds_real'] - T[a]['dpds']:+.5f}")

    B = {}
    for a in ("K3", "Kinf_shuf", "gi@0.25+Kinf", "Kinf+gf@0.25", "Kinf+sh@1.05", "Kinf+sh@1.15", "Kinf+sh@1.25",
              "Kinf+shL@1.15"):
        B[a] = {f: boot(f, a, "Kinf") for f in (*CANCER, "kolf2")}
    out["boot_vs_Kinf"] = B
    P("\npaired target bootstrap vs Kinf (2000 draws; rel cos CI, P>0 | dpds CI, P>0):")
    for a, per in B.items():
        P(f"  {a:16s} " + "  ".join(f"{f[:8]} rel [{b['rel_ci'][0]:+.3f},{b['rel_ci'][1]:+.3f}] P{b['rel_p']:.2f} "
                                    f"pds [{b['dpds_ci'][0]:+.4f},{b['dpds_ci'][1]:+.4f}] P{b['dpds_p']:.2f}"
                                    for f, b in per.items() if b))

    # ------------------------------------------------------------------ (c) amplitude on 410c's bulk
    AB = {c: dict(A=dec["410b"][c]["A"] * T["Kinf"]["Ar"][c] / T["K3"]["Ar"][c],
                  B=dec["410b"][c]["B"] * T["Kinf"]["Br"][c] / T["K3"]["Br"][c]) for c in CTX}
    X410c = Xb["410b"] * (1 + k_key * T["Kinf"]["rel"]) / (1 + k_key * rel3) * np.sqrt(
        np.mean([AB[c]["A"] for c in CTX]) / Abar["410b"])
    r0, _ = raw_of(AB, Kc, claim, X410c, sumD)
    P(f"\n(c) AMPLITUDE on the 410c bulk (predicted X {X410c:.6f}, raw {r0:.6f}); rho_c "
      f"{[round(AB[c]['B'] / claim[c], 5) for c in CTX]}")
    grid = np.linspace(0.95, 1.07, 1201)
    amp = dict(X410c=X410c, raw410c=r0, per_ctx={})
    for c in CTX:
        f0 = lambda s, X: mse_c(s * s * AB[c]["A"], s * s * AB[c]["B"], Kc[c], claim[c], s * X, sumD)
        dv = np.array([f0(s, X410c) - f0(1.0, X410c) for s in grid])
        dv395 = np.array([f0(s, Xb["395"]) - f0(1.0, Xb["395"]) for s in grid])
        cl = float(np.sqrt(claim[c] / AB[c]["B"]))
        draws = cl * np.exp(rng.normal(0, SD187[c], 200))
        def ev(s):
            v = 0.0
            for d in draws:
                cm = d * d * AB[c]["B"]
                v += (s * s * AB[c]["A"] - min(cm, s * s * AB[c]["B"]) - 2 * s * X410c) - (AB[c]["A"] - min(cm, AB[c]["B"]) - 2 * X410c)
            return v / len(draws) / sumD
        evs = np.array([ev(s) for s in grid[::5]])
        sl_b = (2 * (AB[c]["A"] - AB[c]["B"]) * cl - 2 * X410c) / sumD
        sl_a = (2 * AB[c]["A"] * cl - 2 * X410c) / sumD
        amp["per_ctx"][c] = dict(cliff_s=cl, alpha_equiv=0.982 * cl, argmin_s=float(grid[dv.argmin()]),
                                 argmin_s_X395=float(grid[dv395.argmin()]), d_raw_opt=float(dv.min()),
                                 ev_argmin_s=float(grid[::5][evs.argmin()]), ev_d_raw=float(evs.min()),
                                 slope_below=float(sl_b), slope_above=float(sl_a))
        P(f"  {c}: cliff at bulk scale {cl:.4f} (alpha-equivalent {0.982 * cl:.4f}); argmin {grid[dv.argmin()]:.4f} "
          f"(same with 395's X: {grid[dv395.argmin()]:.4f}); d raw on context {dv.min():+.5f}; EV-optimal "
          f"{grid[::5][evs.argmin()]:.4f} ({evs.min():+.5f}); d raw per unit scale below/above the cliff {sl_b:+.3f}/{sl_a:+.3f}")
    tot = float(np.mean([amp["per_ctx"][c]["d_raw_opt"] for c in CTX]))
    tev = float(np.mean([amp["per_ctx"][c]["ev_d_raw"] for c in CTX]))
    gl = np.array([raw_of(AB, Kc, claim, X410c, sumD, {c: s for c in CTX})[0] - r0 for s in grid])
    P(f"  per-context to the cliffs: d raw {tot:+.5f} -> overall {MSE_SC * tot / 6:+.5f}; EV-optimal {tev:+.5f} -> "
      f"{MSE_SC * tev / 6:+.5f}; context A alone {amp['per_ctx']['A']['d_raw_opt'] / 3:+.5f} -> "
      f"{MSE_SC * amp['per_ctx']['A']['d_raw_opt'] / 18:+.5f}")
    P(f"  one GLOBAL bulk scale: argmin {grid[gl.argmin()]:.4f}, d raw {gl.min():+.5f} -> overall {MSE_SC * gl.min() / 6:+.5f}")
    amp.update(to_cliffs=tot, ev=tev, global_argmin=float(grid[gl.argmin()]), global_d_raw=float(gl.min()))
    out["amplitude"] = amp

    (WORK / "analyse414.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "analyse414.log").write_text("\n".join(L))


if __name__ == "__main__":
    main()
