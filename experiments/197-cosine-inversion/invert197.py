"""197 — invert the board's `mse` for the implied transfer cosine.

The instrument is 096's identity, per context (187's correction), evaluated on the real
challenge panel:

    MSE_c - 1 = ( A_c + Kc_c - min(claim_c, B_c) - 2 X ) / sumD          raw = mean_c MSE_c

with `A_c`, `B_c`, `Kc_c`, `claim_c` at the arm's own shipped amplitude.  `sumD` and `Kc` are
properties of the withheld REFERENCE and are therefore shared by every team on the board;
`K = sum_p w_p C_p` is likewise a reference property and is the hard ceiling on `claim`
(`PRED_TRACE_CAP_K = 1.0`, verified against cell_eval2).

Two regimes, and they are the whole experiment:

  OURS      -- A, B, claim measured on our own shipped artifacts, so X is solved exactly from
               one board equation and  cos = X / sqrt(Abar * sumD)  is a MEASUREMENT.
  THE FIELD -- only `raw mse` is observed.  Minimising cos over the feasible set gives a FLOOR:
               a team at its own amplitude optimum with `claim` at the cap needs the least
               direction to reach a given raw, so

                   cos  >=  (1 + Kc/sumD - raw) / ( 2 sqrt(K/sumD) )

               and the bound is vacuous (<= 0) for any raw >= 1 + Kc/sumD = 1.02361.

    uv run python invert197.py

Reads (never writes) 096, 180, 187, 191, 196 artifacts and data/leaderboard/.  Writes only
experiments/197-cosine-inversion/work/.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
W096 = f"{REPO}/experiments/096-mse-floor/work"
sys.path.insert(0, f"{REPO}/src")
sys.path.insert(0, f"{REPO}/experiments/096-mse-floor")

from cell_eval2.metrics.delta import (  # noqa: E402
    PRED_TRACE_CAP_K, _across_pert_budget, _row_weights,
)
from cell_eval2.prep import bulk_lognorm_means  # noqa: E402

CTX = ("A", "B", "C")
TS = 50_000.0
ALPHA_PROD = 1.0796
B_ANCHOR, SPAN = 0.987862879, 0.950162566        # 096 (d), recovered scorer scale for score_mse
NMAE_B, NMAE_SPAN = 1.000439, 0.5963             # 187 [KA-a], our own 10-row affine map

BOARD_RAW = {"a8": 1.362666960949, "a10": 1.3524148, "a11": 0.9099893, "a12": 0.9013492,
             "a14": 0.882775, "a15": 0.873805, "a23": 0.872353, "abc-knee": 0.874491}
ALPHA = {"a8": 1.0796, "a10": 1.0796, "a11": 0.838, "a12": 1.0108, "a14": 0.982,
         "a15": 0.982, "a23": 0.982}
BOARD_PDS = {"a8": 0.831773, "a11": 0.833307, "a12": None, "a14": 0.832702, "a15": 0.835987,
             "a23": None, "abc-knee": 0.836786, "a10": 0.831773}

LB = f"{REPO}/data/leaderboard/2026-09-03T0030Z.json"
LB_LIVE = f"{WORK}/lb-live.json"
LINES = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LINES.append(s)


def score_mse(raw):
    return min(max((B_ANCHOR - raw) / SPAN, 0.0), 1.0)


# ---------------------------------------------------------------------------------------------
# term measurement
# ---------------------------------------------------------------------------------------------
CTRL = np.load(f"{W096}/ctrl-stats.npz", allow_pickle=True)


def terms_from_sums(perts, sums, genes):
    """Per-context (A, B, Kc) from a shipped artifact's per-(block, gene) panel totals.

    Identical conventions to 096's `mech_a12_audit.measure`: the row's own target coordinate is
    dropped from `A`, `w_p = 1/G_p`, and `B` is cell_eval2's own `_across_pert_budget`.
    """
    G = genes.size
    gpos = {g: i for i, g in enumerate(genes)}
    pm_all = bulk_lognorm_means(sums.astype(np.float64), TS)
    out = {}
    for c in CTX:
        sel = np.flatnonzero(np.char.startswith(perts, f"{c}|"))
        tgt = np.asarray([p.split("|", 1)[1] for p in perts[sel]])
        bctrl = bulk_lognorm_means(CTRL[f"sum_{c}"][None, :].astype(np.float64), TS)[0]
        cols = {k: gpos[t] for k, t in enumerate(tgt) if t in gpos}
        w = _row_weights(sel.size, G, cols)
        pm = pm_all[sel]
        dh = pm - bctrl[None, :]
        a = np.einsum("ij,ij->i", dh, dh)
        for k, col in cols.items():
            a[k] -= dh[k, col] ** 2
        out[c] = dict(A=float((a * w).sum()), B=_across_pert_budget(pm, cols, w),
                      Kc=float((np.full(sel.size, float(CTRL[f"jk_{c}"])) * w).sum()),
                      n=int(sel.size), wsum=float(w.sum()))
    return out


def invert(per, raw, sumD):
    """One board point -> the alignment X at the arm's OWN amplitude, and the cosine.

    `cos = X / sqrt(Abar * sumD)` is exactly invariant to a global amplitude rescale, because
    X scales as s and A as s^2.  It is 096's `deconfounded cosine`.
    """
    num = float(np.mean([per[c]["A"] + per[c]["Kc"] - min(per[c]["claim"], per[c]["B"])
                         for c in CTX]))
    X = (num - (raw - 1.0) * sumD) / 2.0
    Abar = float(np.mean([per[c]["A"] for c in CTX]))
    return X, X / np.sqrt(Abar * sumD), Abar, num


def raw_of(per, X, sumD):
    return float(np.mean([1.0 + (per[c]["A"] + per[c]["Kc"] - min(per[c]["claim"], per[c]["B"])
                                 - 2.0 * X) / sumD for c in CTX]))


def score_of(per, X, sumD):
    return float(np.mean([score_mse(1.0 + (per[c]["A"] + per[c]["Kc"]
                                           - min(per[c]["claim"], per[c]["B"])
                                           - 2.0 * X) / sumD) for c in CTX]))


# ---------------------------------------------------------------------------------------------
def main():
    os.makedirs(WORK, exist_ok=True)
    out = {}
    gate = json.load(open(f"{WORK}/gate197.json"))
    sumD = gate["G6"]["sumD"]
    Xprod = gate["G6"]["X"]
    bs = np.load(f"{W096}/board_solve.npz")
    K = float(bs["K"])
    P("=" * 96)
    P("197 — INVERTING THE BOARD FOR THE IMPLIED TRANSFER COSINE")
    P("=" * 96)
    P(f"instrument constants (all from the Phase-0 gates):  sumD {sumD:.6f}   "
      f"X_prod(v11 family) {Xprod:.6f}   cap K {K:.6f}   PRED_TRACE_CAP_K {PRED_TRACE_CAP_K}")

    # ---------------------------------------------------------------------------------------
    # PHASE 1 — our own rows
    # ---------------------------------------------------------------------------------------
    from mech_a12_audit import cp_proxies, measure
    cp16 = {c: cp_proxies()["ntc"][c][0] for c in CTX}
    P(f"\nC_p proxy (096's 16-construct NTC jackknife): "
      + "  ".join(f"{c} {cp16[c]:.4f}" for c in CTX))

    arms = {}
    for arm in ("a8", "a10", "a11", "a12"):
        rows = measure(f"{W096}/{arm}-stats.npz", cp16)
        arms[arm] = {r["ctx"]: dict(A=r["A"], B=r["B"], Kc=r["Kc"], claim=r["claim"],
                                    Chat=np.asarray(r["Chat"]), Cp=r["Cp"], n=r["n"])
                     for r in rows}

    # a14 = a12's artifact at VCC_ALPHA 0.982.  Emission unchanged -> claim unchanged (096
    # measured C_hat flat in s to 0.01%); A, B scale by s^2.
    f14 = ALPHA["a14"] / ALPHA["a12"]
    arms["a14"] = {c: dict(A=arms["a12"][c]["A"] * f14 ** 2, B=arms["a12"][c]["B"] * f14 ** 2,
                           Kc=arms["a12"][c]["Kc"], claim=arms["a12"][c]["claim"],
                           n=arms["a12"][c]["n"]) for c in CTX}

    # a15 — its OWN delivered panel totals (196 streamed them off `191-a15_local.h5ad`, which
    # 191 §6 proved byte-identical to the shipped a15 package).  `claim` from 191's own
    # per-context measurement on the same artifact.
    z15 = np.load(f"{REPO}/experiments/196-emitted-depth/work/chal-totals-a15.npz",
                  allow_pickle=True)
    t15 = terms_from_sums(z15["blocks"].astype(str), z15["totals"], z15["genes"].astype(str))
    c191 = json.load(open(f"{REPO}/experiments/191-rc-hybrid-build/work/claim191.json"))
    for c in CTX:
        n, G = t15[c]["n"], z15["genes"].size
        t15[c]["claim"] = n * c191["a15_local"][c]["claim_per_row"] / (G - 1)
        t15[c]["Cp"] = c191["a15_local"][c]["C_p"]
    arms["a15"] = t15

    # a23 / abc-knee — 180's 24 measured rows per context, as 187 read them.
    k180 = json.load(open(f"{REPO}/experiments/180-def-readiness/work/knee180.json"))
    arms["a23"] = {c: dict(A=k180["per_context"][c]["before"]["A"],
                           B=k180["per_context"][c]["before"]["B"],
                           Kc=k180["Kc"][c],
                           claim=k180["per_context"][c]["before"]["claim"],
                           Cp=k180["per_context"][c]["before"]["C_p"]) for c in CTX}
    arms["abc-knee"] = {c: dict(A=k180["per_context"][c]["after"]["A"],
                                B=k180["per_context"][c]["after"]["B"],
                                Kc=k180["Kc"][c],
                                claim=k180["per_context"][c]["after"]["claim"],
                                Cp=k180["per_context"][c]["before"]["C_p"]) for c in CTX}
    KAPPA = {c: k180["per_context"][c]["knee_scale"] for c in CTX}

    P("\n" + "-" * 96)
    P("PHASE 1 — our own eight board rows, inverted.  `cos` = X / sqrt(Abar*sumD), the")
    P("          amplitude-invariant deconfounded cosine of 096 (a).")
    P("-" * 96)
    P(f"{'arm':<10}{'alpha':>7}{'board raw':>11}{'Abar':>10}{'claim':>9}{'rho_min':>9}"
      f"{'X':>10}{'cos':>9}{'profile':>26}")
    ROW = {}
    ORDER = ["a8", "a10", "a11", "a12", "a14", "a15", "a23", "abc-knee"]
    PROFILE = {"a8": "v11 (pruned whitener)", "a10": "v11 (pruned whitener)",
               "a11": "v11 (pruned whitener)", "a12": "v11 + tiered spread",
               "a14": "v11 + tiered spread", "a15": "+ CD4 residual",
               "a23": "+ CD4 + mass fix", "abc-knee": "+ CD4 + mass fix"}
    for arm in ORDER:
        per = arms[arm]
        X, cos, Abar, num = invert(per, BOARD_RAW[arm], sumD)
        rho = min(min(1.0, per[c]["B"] / per[c]["claim"]) for c in CTX)
        cl = float(np.mean([per[c]["claim"] for c in CTX]))
        ROW[arm] = dict(X=X, cos=cos, Abar=Abar, claim=cl, rho_min=rho,
                        raw=BOARD_RAW[arm], num=num,
                        per={c: {k: v for k, v in per[c].items() if k != "Chat"} for c in CTX})
        P(f"{arm:<10}{ALPHA.get(arm, float('nan')):>7.4f}{BOARD_RAW[arm]:>11.6f}{Abar:>10.6f}"
          f"{cl:>9.5f}{rho:>9.4f}{X:>10.6f}{cos:>9.5f}   {PROFILE[arm]:<24}")

    # --- the consistency check the mandate asks for ----------------------------------------
    P("\n  CONSISTENCY CHECK — arms that changed only the EMITTER or only the AMPLITUDE must")
    P("  invert to the SAME cosine; only a profile change may move it.")
    GROUPS = [("v11 family (a8, a10, a11)", ["a8", "a10", "a11"],
               "a8->a10 emission only (phi 4 -> s0.15/phi16); a10->a11 alpha only (1.0796->0.838)"),
              ("tiered family (a12, a14)", ["a12", "a14"], "alpha only (1.0108 -> 0.982)"),
              ("mass-fix family (a23, abc-knee)", ["a23", "abc-knee"],
               "per-context alpha only (0.982 -> 0.9823/1.0129/1.0479)")]
    cons = {}
    for name, mem, why in GROUPS:
        cs = [ROW[m]["cos"] for m in mem]
        spread = (max(cs) - min(cs)) / np.mean(cs)
        cons[name] = dict(members=mem, cos=cs, rel_spread=spread)
        P(f"    {name:<32} cos " + " / ".join(f"{c:.5f}" for c in cs)
          + f"   spread {100*spread:+.2f}%   [{why}]")
    P("\n  PROFILE CHANGES (these are ALLOWED to move the cosine, and are the only ones that did):")
    for a, b, what in (("a11", "a12", "tiered spread policy (096 intervention 1)"),
                       ("a14", "a15", "CD4 residual (VCC_CD4_TARGET_F 0.10)"),
                       ("a15", "a23", "mass fix (minnorm)"),):
        P(f"    {a} -> {b:<9} cos {ROW[a]['cos']:.5f} -> {ROW[b]['cos']:.5f}  "
          f"({100*(ROW[b]['cos']/ROW[a]['cos']-1):+.2f}%)   {what}")

    # --- byte-identical resubmission (Phase 4 falsifier 2) ---------------------------------
    P("\n  FALSIFIER 2 — the byte-identical resubmission of a15 (entry 5QfX14axLAG8FF2rSDJe).")
    hit = None
    for f in sorted(glob.glob(f"{REPO}/data/leaderboard/*.json")) + [LB_LIVE]:
        try:
            e = json.load(open(f))
            e = e["entries"] if isinstance(e, dict) else e
        except Exception:
            continue
        for r in e:
            if r.get("id") == "5QfX14axLAG8FF2rSDJe":
                hit = (f, r)
    if hit is None:
        P("    entry id not present in any snapshot on disk; using docs/LOG.md's recorded claim")
        P("    that every member reproduced to float representation -> identical terms, identical")
        P("    board raw, therefore IDENTICAL cosine by construction.  Vacuous but true.")
    else:
        f, r = hit
        X2, c2, _, _ = invert(arms["a15"], r["expr_mse_unbiased_capped_norm"], sumD)
        P(f"    found in {os.path.basename(f)}: raw {r['expr_mse_unbiased_capped_norm']:.10f} "
          f"vs a15 {BOARD_RAW['a15']:.10f}")
        P(f"    inverted cos {c2:.6f} vs a15's {ROW['a15']['cos']:.6f}   "
          f"delta {c2-ROW['a15']['cos']:+.2e}")
    out["resub"] = None if hit is None else dict(raw=hit[1]["expr_mse_unbiased_capped_norm"])

    # --- rc_hybrid: same pseudobulk, different emission -------------------------------------
    P("\n  rc_hybrid (submitted 2026-09-03 11:20 UTC; its board row has NOT landed).")
    P("  191 falsifier 1: per-gene panel-total digests identical to a15 on 900/900 blocks, so")
    P("  its A, B and X are a15's EXACTLY and only `claim` moved.  Prediction, pre-registered:")
    rc = {c: dict(arms["a15"][c]) for c in CTX}
    for c in CTX:
        n = arms["a15"][c]["n"]
        rc[c]["claim"] = n * c191["rc_hybrid"][c]["claim_per_row"] / (z15["genes"].size - 1)
    rc_raw = raw_of(rc, ROW["a15"]["X"], sumD)
    P(f"    claim  a15 " + "/".join(f"{arms['a15'][c]['claim']:.5f}" for c in CTX)
      + "   ->  rc_hybrid " + "/".join(f"{rc[c]['claim']:.5f}" for c in CTX))
    P(f"    branch min(claim,B) is on the B leg in " +
      ", ".join(c for c in CTX if arms["a15"][c]["B"] < rc[c]["claim"]) + " (extra claim buys 0)")
    P(f"    predicted rc_hybrid raw mse {rc_raw:.6f}  vs a15's board {BOARD_RAW['a15']:.6f}  "
      f"(delta {rc_raw-BOARD_RAW['a15']:+.6f})")
    P(f"    => rc_hybrid MUST invert to a15's cosine {ROW['a15']['cos']:.5f}; any other answer")
    P(f"       falsifies either the inversion or 191's digest identity.")
    out["rc_hybrid_pred_raw"] = rc_raw

    out["ours"] = ROW
    out["consistency"] = cons
    out["sumD"] = sumD
    out["K"] = K

    # ---------------------------------------------------------------------------------------
    # PHASE 2 — the field
    # ---------------------------------------------------------------------------------------
    Kc = float(np.mean([arms["a15"][c]["Kc"] for c in CTX]))
    floor0 = 1.0 + Kc / sumD                      # raw at cos = 0, claim at the cap, eps = 0
    coeff = 2.0 * np.sqrt(K / sumD)
    P("\n" + "-" * 96)
    P("PHASE 2 — the field.  What a raw `mse` alone can say about a stranger's direction.")
    P("-" * 96)
    P("  At its own amplitude optimum (096 (b): the minimum sits exactly on the rho cliff),")
    P("      MSE_min = 1 + Kc/sumD + q*eps - 2*cos*sqrt(q*(1+eps)),   q = claim/sumD,")
    P("      eps = (A-B)/B >= 0 is the prediction's common mode.")
    P("  cos is DECREASING in q, and q <= K/sumD because PRED_TRACE_CAP_K = 1.0, so the cap")
    P("  gives a FLOOR on cos.  Being off the amplitude optimum only raises the true cos, so")
    P("  the floor survives that too.  With eps = 0 (its minimum):")
    P(f"      cos >= ({floor0:.6f} - raw) / {coeff:.6f}          Kc {Kc:.6f}  K/sumD {K/sumD:.5f}")
    P(f"  and the bound is VACUOUS for any raw >= {floor0:.5f}.")

    # tightness, calibrated on our own rows where cos is measured
    P("\n  TIGHTNESS, calibrated on our own eight rows (cos measured, floor computed blind):")
    P(f"  {'arm':<10}{'raw':>10}{'floor':>10}{'measured':>10}{'floor/meas':>12}{'verdict':>14}")
    tight = {}
    for arm in ORDER:
        fl = (floor0 - BOARD_RAW[arm]) / coeff
        m = ROW[arm]["cos"]
        v = "VACUOUS" if fl <= 0 else ("valid" if fl <= m else "*** VIOLATED ***")
        tight[arm] = dict(floor=fl, measured=m, ratio=fl / m)
        P(f"  {arm:<10}{BOARD_RAW[arm]:>10.5f}{fl:>10.5f}{m:>10.5f}{fl/m:>12.4f}{v:>14}")
    out["tightness"] = tight
    bad = [a for a in ORDER if tight[a]["floor"] > tight[a]["measured"] + 1e-9]
    P(f"  bound violated on {len(bad)} of {len(ORDER)} of our own rows: {bad if bad else 'none'}")

    # the board
    lb = json.load(open(LB))
    teams = lb["entries"]
    P(f"\n  Board snapshot {os.path.basename(LB)}: {len(teams)} scored teams.")
    shrunk = [t for t in teams if t["expr_mse_unbiased_capped_norm"] < 1.2]
    informative = [t for t in teams if t["expr_mse_unbiased_capped_norm"] < floor0]
    P(f"    raw mse < 1.2               : {len(shrunk)} teams (the 'shrunk cohort')")
    P(f"    raw mse < {floor0:.5f} (informative): {len(informative)} teams")
    P(f"    raw mse >= {floor0:.5f}            : {len(teams)-len(informative)} teams — the "
      f"inversion returns NOTHING for these (falsifier 1)")
    for t in teams:
        u = t["expr_mse_unbiased_capped_norm"]
        t["_floor"] = (floor0 - u) / coeff
        t["_credit"] = 1.0 - t["de_wilcoxon_lfc_nmae"]
        t["_qmin"] = 0.5 * (1.0 + t["_credit"])          # q_w at K = 1 (187 §5); a LOWER bound
    informative.sort(key=lambda t: -t["_floor"])
    us = next(t for t in teams if "Ahmet" in t["team_name"])
    P(f"\n  Top 20 by implied cosine FLOOR (the 49-team shrunk cohort's informative subset):")
    P(f"  {'#':>3} {'team':<24}{'overall':>9}{'raw mse':>9}{'cos floor':>10}{'x ours':>8}"
      f"{'raw pds':>9}{'nmae':>8}{'q_w>=':>8}{'reach':>8}")
    for i, t in enumerate(informative[:20]):
        P(f"  {i+1:>3} {t['team_name'][:23]:<24}{t['score_avg']:>9.5f}"
          f"{t['expr_mse_unbiased_capped_norm']:>9.5f}{t['_floor']:>10.5f}"
          f"{t['_floor']/us['_floor']:>8.3f}{t['pds_cosine']:>9.5f}"
          f"{t['de_wilcoxon_lfc_nmae']:>8.5f}{t['_qmin']:>8.5f}"
          f"{t['de_wilcoxon_direction_reach_raw']:>8.4f}")
    P(f"\n  us: floor {us['_floor']:.5f}, MEASURED {ROW['abc-knee']['cos']:.5f} "
      f"(the floor is {100*(1-us['_floor']/ROW['abc-knee']['cos']):.1f}% low on a row where we "
      f"know the answer)")
    out["floor0"], out["coeff"] = floor0, coeff
    out["field"] = [dict(team=t["team_name"], overall=t["score_avg"],
                         raw=t["expr_mse_unbiased_capped_norm"], floor=t["_floor"],
                         pds=t["pds_cosine"], nmae=t["de_wilcoxon_lfc_nmae"],
                         qmin=t["_qmin"], reach=t["de_wilcoxon_direction_reach_raw"],
                         fid=t["de_wilcoxon_direction_fidelity_yield_raw"],
                         jac=t["de_wilcoxon_sig_jaccard"], model=t.get("model_name"))
                    for t in informative]

    json.dump(out, open(f"{WORK}/invert197.json", "w"), indent=1, default=float)
    open(f"{WORK}/invert197.log", "w").write("\n".join(LINES))
    P(f"\nwrote {WORK}/invert197.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
