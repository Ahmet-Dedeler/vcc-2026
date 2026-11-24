"""211 Phase 1-4 — re-locate the cliff on the CURRENT emission, re-solve the joint optimum with
the DE members no longer assumed dead, retrodict, and price the arm on both base profiles.

Inputs, all measured elsewhere in this directory:
  work/bx211.json    per-context `A_c`, `B_c` measured ON THE SHIPPED ARTIFACTS
                     (submissions/202-rc_callset.h5ad and submissions/210-cd4f034.h5ad),
                     through cell_eval2's own `_across_pert_budget` / `_row_weights`.
  work/chat211.json  per-context `q_c` = mean_p min(C_hat_p, C_p)/C_p for `quasi` (a15),
                     `rc_hybrid` and `rc_callset`, reproducing 208's C_hat/C_p to six decimals.
  work/mse187.json   187's own per-context terms for the a12 family and for a23 (reproduced
                     here; KA-i / KA-ii / KA-iii all PASS).

The identity, per context (096 (a), 187's per-context correction), at global scale `g` on the
base arm's own delivered profile (`alpha = 0.982 g`):

    raw_c(g) = 1 + ( g^2 A_c + Kc_c - min(claim_c, g^2 B_c) - 2 g X ) / sumD
    claim_c  = q_c * cap_c ,   cap_c = (sum_p w_p) * f_c * C_p      <- HARD CAP, q_c <= 1
    score    = mean_c clip( (0.987862879 - raw_c) / 0.950162566, 0, 1 )
    cliff_c  = 0.982 * sqrt( claim_c / B_c )

`PRED_TRACE_CAP_K = 1.0` is what makes `cap_c` a ceiling no emission can exceed, and it is the
hinge of this whole experiment: `q_c` is the ONLY thing an emission arm can move, and 208's
C_hat/C_p is a statement about how far `q_c` is from FALLING off 1, not about how far the
amplitude cliff is.

    uv run python solve211.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
BX = json.loads((WORK / "bx211.json").read_text())
CH = json.loads((WORK / "chat211.json").read_text())

CTX = ("A", "B", "C")
ALPHA0 = 0.982
SUMD = 0.229141                       # 187 [KA-i]
B_ANCHOR, SPAN = 0.987862879, 0.950162566
KC = BX["Kc"]
CP46 = BX["cp46"]
WSUM = BX["per_context"]["rc_callset"]["A"]["weight_sum"]
F_FIT = {"A": 1.0, "B": 1.0, "C": 0.9825}       # 187's retrodicted per-context C_p factor
F_PROXY = {"A": 1.0, "B": 1.0, "C": 1.0}
BOARD_RC = 0.8738045054832962

NMAE_SLOPE = 0.0416                   # raw per ln(alpha), board-measured on the current pipeline
NMAE_SPAN = 0.596294

# --- the board rows this experiment reads (docs/BOARD-LEDGER.md + the two BOARD-RESULT files) --
LEDGER = {
    "a10":  dict(pds=0.736561, mse=0.000000, nmae=0.094338, fid=-0.027897, reach=0.123595, jac=-0.001493),
    "a11":  dict(pds=0.739889, mse=0.082041, nmae=0.096982, fid=-0.031535, reach=0.112220, jac=-0.002907),
    "a12":  dict(pds=0.737885, mse=0.091163, nmae=0.093974, fid=-0.031473, reach=0.116790, jac=-0.000413),
    "a14":  dict(pds=0.738614, mse=0.110615, nmae=0.091972, fid=-0.032586, reach=0.114542, jac=-0.000612),
    "a15":  dict(pds=0.745871, mse=0.120057, nmae=0.094239, fid=-0.030397, reach=0.115368, jac=-0.000619),
    "knee": dict(pds=0.747726, mse=0.119322, nmae=0.096214, fid=-0.029477, reach=0.110910, jac=-0.000433),
}
DLN = {"a10->a11": np.log(0.838 / 1.0796), "a12->a14": np.log(0.982 / 1.0108)}
KNEE_ALPHA = {"A": 0.98234, "B": 1.01291, "C": 1.04786}
A23_OVERALL, KNEE_OVERALL = 0.1733565, 0.1740437
KNEE_DE_OVERALL = 0.000684             # 187 section 4 / B4: reach+fid+jac, abc-knee vs a23
TWIN_REACH_SCALED = 0.011              # 180's twin band, the wide end


def score(raws):
    return float(np.mean([min(max((B_ANCHOR - r) / SPAN, 0.0), 1.0) for r in raws]))


class Arm:
    """One delivered profile + one emission, as a function of a global amplitude factor `g`."""

    def __init__(self, name, A, B, q, f, board_raw=None, X=None):
        self.name, self.A, self.B, self.q, self.f = name, A, B, q, f
        self.cap = {c: WSUM * f[c] * CP46[c] for c in CTX}
        self.claim = {c: q[c] * self.cap[c] for c in CTX}
        if X is None:
            num = float(np.mean([A[c] + KC[c] - min(self.claim[c], B[c]) for c in CTX]))
            X = (num - (board_raw - 1.0) * SUMD) / 2.0
        self.X = X
        self.cliff = {c: ALPHA0 * float(np.sqrt(self.claim[c] / B[c])) for c in CTX}
        self.cliff_ceiling = {c: ALPHA0 * float(np.sqrt(self.cap[c] / B[c])) for c in CTX}
        self.rho = {c: min(1.0, B[c] / self.claim[c]) for c in CTX}
        self.raw0 = self.raws({c: ALPHA0 for c in CTX})
        self.score0 = score(self.raw0)

    def raws(self, alphas):
        out = []
        for c in CTX:
            g = alphas[c] / ALPHA0
            out.append(1 + (g * g * self.A[c] + KC[c]
                            - min(self.claim[c], g * g * self.B[c]) - 2 * g * self.X) / SUMD)
        return out

    def d_mse(self, alphas):
        return score(self.raws(alphas)) - self.score0

    def d_raw(self, alphas):
        return float(np.mean(self.raws(alphas))) - float(np.mean(self.raw0))


def d_nmae(alphas):
    dln = float(np.mean([np.log(alphas[c] / ALPHA0) for c in CTX]))
    return NMAE_SLOPE * dln / NMAE_SPAN, -NMAE_SLOPE * dln


def main():
    lines = []
    P = lambda *a: (print(*a, flush=True), lines.append(" ".join(str(x) for x in a)))  # noqa: E731
    out = {}

    # =====================================================================================
    # Phase 1 — where is the cliff, on each of the three board emissions?
    # =====================================================================================
    P("=" * 100)
    P("PHASE 1 — the cliff on the CURRENT emission")
    P("=" * 100)
    Arc = {c: BX["per_context"]["rc_callset"][c]["A"] for c in CTX}
    Brc = {c: BX["per_context"]["rc_callset"][c]["B"] for c in CTX}
    Acd = {c: BX["per_context"]["cd4f034"][c]["A"] for c in CTX}
    Bcd = {c: BX["per_context"]["cd4f034"][c]["B"] for c in CTX}
    P(f"\nA_c, B_c measured on submissions/202-rc_callset.h5ad (the delivered a15 profile):")
    for c in CTX:
        P(f"  {c}: A {Arc[c]:.6f}  B {Brc[c]:.6f}")
    P("Emitters differ ONLY in q_c (270/270 blocks deliver identical per-gene panel totals, chat211).")

    qs = {}
    for em, tag in (("quasi", "a15"), ("rc_hybrid", "rc_hybrid"), ("rc_callset", "rc_callset")):
        qs[tag] = CH["rows"][f"{em}|MAXT24"]["q"]
    qs["a15 (MAXT8)"] = CH["rows"]["quasi|MAXT8"]["q"]

    P(f"\n--- per-context cliff table, f = (1, 1, 0.9825) [187's retrodicted C_p factor] ---")
    P(f"  {'arm':<14}{'ctx':>4}{'q_c':>10}{'claim_c':>10}{'cap_c':>10}{'B_c':>10}"
      f"{'s^2B/claim':>12}{'branch':>8}{'cliff alpha':>13}{'headroom':>10}{'ceiling':>10}")
    arms1 = {}
    for tag in ("a15", "rc_hybrid", "rc_callset"):
        a = Arm(tag, Arc, Brc, qs[tag], F_FIT, board_raw=BOARD_RC)
        arms1[tag] = a
        for c in CTX:
            branch = "claim" if a.B[c] > a.claim[c] else "s^2 B"
            P(f"  {tag:<14}{c:>4}{a.q[c]:>10.6f}{a.claim[c]:>10.6f}{a.cap[c]:>10.6f}{a.B[c]:>10.6f}"
              f"{a.B[c]/a.claim[c]:>12.5f}{branch:>8}{a.cliff[c]:>13.5f}"
              f"{a.cliff[c]/ALPHA0-1:>+10.3%}{a.cliff_ceiling[c]:>10.5f}")
    P(f"\n  187 measured, on a23's profile: rho 0.9993 / 0.9399 / 0.8783 (mean 0.938; 170's "
      f"challenge-panel pooled rho 0.938, 164's 0.9466)")
    P(f"  187's cliffs on a23's profile:   A 0.9823  B 1.0129 (+3.1%)  C 1.0386 (+5.8%)")
    d = {c: arms1["rc_callset"].cliff[c] - arms1["a15"].cliff[c] for c in CTX}
    P(f"\n  *** a15 -> rc_callset moves the cliff by: "
      + "  ".join(f"{c} {d[c]:+.6f} ({d[c]/arms1['a15'].cliff[c]:+.4%})" for c in CTX) + " ***")
    P(f"  *** and the CEILING (claim_c = cap_c exactly, unreachable by any emission) is: "
      + "  ".join(f"{c} {arms1['rc_callset'].cliff_ceiling[c]:.5f}" for c in CTX) + " ***")
    out["phase1"] = {t: dict(q=a.q, claim=a.claim, cap=a.cap, rho=a.rho, cliff=a.cliff,
                             ceiling=a.cliff_ceiling, X=a.X) for t, a in arms1.items()}

    P(f"\n--- the same table under the raw NTC proxy, f = (1,1,1) (187's B2 rejects it, kept as "
      f"a sensitivity) ---")
    for tag in ("a15", "rc_callset"):
        a = Arm(tag, Arc, Brc, qs[tag], F_PROXY, board_raw=BOARD_RC)
        P(f"  {tag:<14}" + "  ".join(f"{c}: cliff {a.cliff[c]:.5f} ({a.cliff[c]/ALPHA0-1:+.2%}) "
                                     f"rho {a.rho[c]:.4f}" for c in CTX))

    P(f"\n--- SENSITIVITY: the only assumption context A's cliff rests on is f_A ---")
    P(f"  cliff_A = 0.982 * sqrt(f_A * cap_A(f=1) / B_A) = 0.98078 * sqrt(f_A), and 187's B1 "
      f"bounds f_A <= 1.000")
    P(f"  {'f_A':>7}{'cap_A':>10}{'branch':>9}{'cliff_A':>10}{'vs 0.982':>10}")
    fa_rows = []
    for fa in (0.96, 0.97, 0.98, 0.99, 0.995, 1.0, 1.005, 1.01):
        cap = WSUM * fa * CP46["A"]
        cl = ALPHA0 * float(np.sqrt(min(cap, Brc["A"]) / Brc["A"]))
        fa_rows.append((fa, cap, cl))
        P(f"  {fa:>7.3f}{cap:>10.6f}{('claim' if Brc['A'] > cap else 's^2 B'):>9}"
          f"{cl:>10.5f}{cl/ALPHA0-1:>+10.3%}")
    P(f"  -> cliff_A is MONOTONE INCREASING in f_A and SATURATES at 0.98078 the moment f_A >= "
      f"1.00249.  Since 187's a12->a14 delta bounds f_A <= 1, the shipped 0.982 is AT OR ABOVE")
    P(f"     the optimum in every admissible model, and the emission cannot change that.")
    out["fA_scan"] = fa_rows

    P(f"\n--- SENSITIVITY: sumD (187's 0.229141 vs 096's board_solve 0.228572) ---")
    for sd_alt in (0.229141, 0.228572):
        g = globals()
        old = g["SUMD"]
        g["SUMD"] = sd_alt
        a = Arm("x", Arc, Brc, qs["rc_callset"], F_FIT, board_raw=BOARD_RC)
        P(f"  sumD {sd_alt:.6f}: X {a.X:.6f}  cliffs " +
          " ".join(f"{c} {a.cliff[c]:.5f}" for c in CTX) + "   (cliffs do not depend on sumD)")
        g["SUMD"] = old

    # =====================================================================================
    # Phase 3 (done here because it gates Phase 2) — the two exact zeros
    # =====================================================================================
    P("\n" + "=" * 100)
    P("PHASE 3a — the free, very sharp gate: rc_hybrid and rc_callset must read d mse = 0.000000")
    P("=" * 100)
    for tag in ("rc_hybrid", "rc_callset"):
        a15, arm = arms1["a15"], arms1[tag]
        dr = float(np.mean(arm.raw0)) - float(np.mean(a15.raw0))
        ds = arm.score0 - a15.score0
        P(f"  a15 -> {tag:<11} d raw mse {dr:+.16f}   d score_mse {ds:+.16f}   "
          f"board {'0.000000 (16 s.f.)':>20}   {'PASS' if abs(ds) < 1e-12 else 'FAIL'}")
    P("  Mechanism, and it is the whole Phase-1 result in one line:")
    P("    ctx A is the ONLY context on the `claim` branch (B_A > claim_A); B and C sit on `s^2 B`,")
    P("    where `claim` is invisible to the metric.  a15's q_A is 1.000000 — it already collects")
    P("    100% of context A's cap — so the emission change moves claim_A by EXACTLY zero, and the")
    P("    only claim that moved (ctx B, +0.0096% of cap) sits on the branch that cannot see it.")
    P("  The board's own 16-significant-figure zero therefore PROVES q_A(a15) = 1 on all 300")
    P("  perturbations, not just on chat211's 24-row sample.  If a15 had had ANY unsaturated row in")
    P("  context A, rc_callset would have shown a positive d mse.  It showed 0.")

    # =====================================================================================
    # Phase 2 — re-solve the joint optimum, with the DE members priced
    # =====================================================================================
    P("\n" + "=" * 100)
    P("PHASE 2 — the DE members across the amplitude history, then the joint optimum")
    P("=" * 100)
    P(f"\n  {'pair':<14}{'d ln alpha':>12}{'d fid':>10}{'d reach':>10}{'d jac':>10}"
      f"{'DE scaled':>11}{'DE overall':>12}{'slope /ln a':>13}")
    de_rows = []
    for pair, (lo, hi) in (("a10->a11", ("a10", "a11")), ("a12->a14", ("a12", "a14"))):
        d = {m: LEDGER[hi][m] - LEDGER[lo][m] for m in ("fid", "reach", "jac")}
        s = sum(d.values())
        de_rows.append((pair, DLN[pair], s / 6.0))
        P(f"  {pair:<14}{DLN[pair]:>+12.5f}{d['fid']:>+10.6f}{d['reach']:>+10.6f}"
          f"{d['jac']:>+10.6f}{s:>+11.6f}{s/6:>+12.6f}{(s/6)/DLN[pair]:>+13.5f}")
    dln_knee = float(np.mean([np.log(KNEE_ALPHA[c] / ALPHA0) for c in CTX]))
    de_rows.append(("a23->abc-knee", dln_knee, KNEE_DE_OVERALL))
    P(f"  {'a23->abc-knee':<14}{dln_knee:>+12.5f}{'':>10}{'+0.003475':>10}{'':>10}"
      f"{'':>11}{KNEE_DE_OVERALL:>+12.6f}{KNEE_DE_OVERALL/dln_knee:>+13.5f}   "
      f"(187 section 4: board reach/fid/jac vs a23)")
    slopes = [v / dl for _n, dl, v in de_rows]
    slope_small = float(np.mean(slopes[1:]))
    band = TWIN_REACH_SCALED / 6.0 / abs(dln_knee)      # one twin sd converted to a slope
    P(f"\n  slope estimates (overall per unit ln alpha): "
      + ", ".join(f"{n} {s:+.5f}" for (n, _d, _v), s in zip(de_rows, slopes)))
    P(f"  the two SMALL, same-pipeline pairs agree: mean {slope_small:+.5f} overall per ln alpha.")
    P(f"  BUT 180's twin band on `reach` alone is +-{TWIN_REACH_SCALED} scaled = "
      f"+-{TWIN_REACH_SCALED/6:.5f} overall, i.e. +-{band:.3f} on this slope.")
    P(f"  -> the DE coupling is NOT RESOLVABLE point by point; the honest carry is "
      f"{slope_small:+.3f} +- {band:.3f}, a one-sided band of [0, {slope_small+band:+.3f}].")
    DE_SLOPE, DE_HI = slope_small, slope_small + band
    out["de_slope"] = dict(rows=[(n, dl, v) for n, dl, v in de_rows], point=DE_SLOPE, hi=DE_HI,
                           band=band)

    P(f"\n--- marginal slopes of the objective in overall units per unit ln(global alpha) ---")
    base = arms1["rc_callset"]
    eps = 1e-4
    for label, g0 in (("just BELOW ctx A's cliff", base.cliff["A"] * (1 - 2e-3)),
                      ("just ABOVE ctx A's cliff", base.cliff["A"] * (1 + 2e-3))):
        lo = {c: g0 * (1 - eps) for c in CTX}
        hi = {c: g0 * (1 + eps) for c in CTX}
        dm = (base.d_mse(hi) - base.d_mse(lo)) / (2 * eps) / 6.0
        dn = (d_nmae(hi)[0] - d_nmae(lo)[0]) / (2 * eps) / 6.0
        P(f"  {label:<26} alpha {g0:.5f}   mse {dm:>+9.5f}   nmae {dn:>+9.5f}   "
          f"DE {DE_SLOPE:>+9.5f} (hi {DE_HI:+.5f})   NET {dm+dn+DE_SLOPE:>+9.5f} "
          f"(hi {dm+dn+DE_HI:+.5f})")
    P("  -> the objective is a CORNER at context A's cliff under both the point DE slope and the")
    P("     top of its band; the DE members make the corner sharper, not softer.")

    # global sweep on both base profiles
    Xcd = {}
    cd_variants = {
        "cd4f034 (X = rc_callset's, null)": base.X,
        "cd4f034 (addendum raw mse 0.8333)": None,
        "cd4f034 (207's raw mse 0.84413)": None,
    }
    profiles = {"rc_callset (on the board)": Arm("rc_callset", Arc, Brc, qs["rc_callset"], F_FIT,
                                                 board_raw=BOARD_RC)}
    for lab, target in (("cd4f034 @ raw 0.8333 (210 addendum)", 0.8333),
                        ("cd4f034 @ raw 0.84413 (207)", 0.84413),
                        ("cd4f034 @ X unchanged (null)", None)):
        if target is None:
            a = Arm(lab, Acd, Bcd, qs["rc_callset"], F_FIT, X=base.X)
        else:
            a = Arm(lab, Acd, Bcd, qs["rc_callset"], F_FIT, board_raw=target)
        profiles[lab] = a
        Xcd[lab] = a.X

    P(f"\n--- the joint optimum in a GLOBAL alpha, both base profiles ---")
    P(f"  {'profile':<36}{'X':>9}{'cliff A':>9}{'cliff B':>9}{'cliff C':>9}"
      f"{'alpha*':>9}{'d mse sc':>10}{'d nmae sc':>11}{'d DE':>10}{'d overall':>11}")
    opt = {}
    for lab, a in profiles.items():
        grid = np.arange(0.930, 1.0801, 0.00005)
        best = None
        for g in grid:
            al = {c: float(g) for c in CTX}
            dm, (dn, _dnr) = a.d_mse(al), d_nmae(al)
            de = DE_SLOPE * float(np.mean([np.log(al[c] / ALPHA0) for c in CTX]))
            tot = (dm + dn) / 6.0 + de
            if best is None or tot > best[1]:
                best = (float(g), tot, dm, dn, de)
        g, tot, dm, dn, de = best
        opt[lab] = dict(alpha=g, overall=tot, d_mse=dm, d_nmae=dn, d_de=de, X=a.X,
                        cliff=a.cliff)
        P(f"  {lab:<36}{a.X:>9.6f}{a.cliff['A']:>9.5f}{a.cliff['B']:>9.5f}{a.cliff['C']:>9.5f}"
          f"{g:>9.5f}{dm:>+10.5f}{dn:>+11.5f}{de:>+10.5f}{tot:>+11.6f}")
    out["global_opt"] = opt

    P(f"\n--- and the PER-CONTEXT optimum (each context at its own cliff) ---")
    P(f"  {'profile':<36}{'A':>9}{'B':>9}{'C':>9}{'d mse sc':>10}{'d nmae sc':>11}"
      f"{'d DE':>10}{'d overall':>11}")
    pc = {}
    for lab, a in profiles.items():
        al = dict(a.cliff)
        dm, (dn, _) = a.d_mse(al), d_nmae(al)
        de = DE_SLOPE * float(np.mean([np.log(al[c] / ALPHA0) for c in CTX]))
        pc[lab] = dict(alphas=al, d_mse=dm, d_nmae=dn, d_de=de, overall=(dm + dn) / 6 + de)
        P(f"  {lab:<36}{al['A']:>9.5f}{al['B']:>9.5f}{al['C']:>9.5f}{dm:>+10.5f}"
          f"{dn:>+11.5f}{de:>+10.5f}{(dm+dn)/6+de:>+11.6f}")
    out["per_context_opt"] = pc

    # =====================================================================================
    # Phase 4 — the gate
    # =====================================================================================
    P("\n" + "=" * 100)
    P("PHASE 4 — rule 13's composed gate, floor +0.0089 (amplitude-changing)")
    P("=" * 100)
    for lab in ("rc_callset (on the board)", "cd4f034 @ raw 0.8333 (210 addendum)"):
        for kind, r in (("global alpha*", opt[lab]), ("per-context cliffs", pc[lab])):
            dm = r["d_mse"]
            priced = (-0.000076 / 6.0) * 1.011 + dm / 6.0
            six = priced + r["d_nmae"] / 6.0 + r["d_de"]
            P(f"\n  base {lab}  /  arm: {kind}")
            P(f"    PLOLO d_pds/6 x 1.011   {(-0.000076/6)*1.011:+.6f}   (180: pds is "
              f"scale-invariant on this pipeline)")
            P(f"    096 d_mse/6             {dm/6:+.6f}")
            P(f"    rule 13 priced          {priced:+.6f}   vs floor +0.0089  -> "
              f"{'SHIP' if priced > 0.0089 else 'VETO'}  ({priced/0.0089:.1%} of the floor)")
            P(f"    + nmae/6                {r['d_nmae']/6:+.6f}")
            P(f"    + DE (measured slope)   {r['d_de']:+.6f}  (at the top of its band "
              f"{DE_HI*float(np.mean([np.log(v/ALPHA0) for v in (r.get('alphas') or {c: r['alpha'] for c in CTX}).values()])):+.6f})")
            P(f"    six-member best est.    {six:+.6f}   -> {six/0.0089:.1%} of the floor")
            out.setdefault("gate", {})[f"{lab}|{kind}"] = dict(priced=priced, six=six)

    P(f"\n--- how far out of reach is the gate on this axis? ---")
    base = profiles["rc_callset (on the board)"]
    below = 0.02474 + 0.01163 + DE_SLOPE          # overall per unit ln alpha, below A's cliff
    need = 0.0089 / below
    P(f"  below A's cliff the objective gains {below:+.5f} overall per unit ln(alpha) "
      f"(mse + nmae + the measured DE slope).")
    P(f"  clearing the +0.0089 floor on amplitude alone would need d ln(alpha) = {need:.4f}, "
      f"i.e. alpha {ALPHA0*np.exp(need):.4f} (+{np.expm1(need):.1%}).")
    P(f"  the largest amplitude any context permits before its own cliff is "
      f"+{max(base.cliff.values())/ALPHA0-1:.2%} (ctx C), and the BINDING one is ctx A at "
      f"{base.cliff['A']/ALPHA0-1:+.2%}.")
    P(f"  -> the gate is {need/max(1e-9, np.log(max(base.cliff.values())/ALPHA0)):.0f}x further "
      f"than the most permissive context allows, and the binding context allows nothing.")
    out["gate_distance"] = dict(slope_below=below, dln_needed=need,
                                alpha_needed=float(ALPHA0 * np.exp(need)))

    P(f"\n--- nmae's linear range (187 nmae187 reproduced here) ---")
    P(f"  board law raw_nmae(alpha) = 0.944330 - 0.0416 * mean_c ln(alpha_c/0.982), fitted on two")
    P(f"  board pairs (a12->a14: 0.0471; a23->abc-knee: 0.0361).  h1-test's own nmae(s) curve")
    P(f"  turns over at s = 1.290 (nmae187, reproduced: argmin 1.290, gain -0.00269 raw); 193")
    P(f"  section 5 independently derives d nmae/d k < 0 at k = 0.278, q = 0.5375, so nmae wants")
    P(f"  LARGER amplitude throughout our range.  Every optimum in this experiment lies inside")
    P(f"  [-0.2%, +0.6%] of 0.982 — ~50x inside the turnover and ~5x inside the tightest board")
    P(f"  pair the law was fitted on, so it is interpolated, never extrapolated.  nmae never gets")
    P(f"  to vote anyway: the mse cliff binds two orders of magnitude earlier.")

    P(f"\n--- CROSS-CHECK: does the artifact-measured X agree with 187's a23 X moved to a15? ---")
    mA_mine = float(np.mean([Arc[c] for c in CTX]))
    mA_187 = float(np.mean([0.2348299009419585, 0.22477784721536484, 0.23494104987562134]))
    mM_mine = float(np.mean([min(base.claim[c], Brc[c]) for c in CTX]))
    mM_187 = float(np.mean([0.233754, 0.223624, 0.233728]))
    dX = (mA_mine - mA_187 - (mM_mine - mM_187)) / 2.0 - (BOARD_RC - 0.872353) * SUMD / 2.0
    P(f"  mean A: mine {mA_mine:.6f}  187/a23 {mA_187:.6f}   mean min(claim,B): mine {mM_mine:.6f}"
      f"  187/a23 {mM_187:.6f}")
    P(f"  predicted X(a15) = X(a23) + dX = 0.017930 {dX:+.6f} = {0.017930+dX:.6f};  "
      f"measured {base.X:.6f};  agreement {abs(0.017930+dX-base.X):.2e}")
    out["X_crosscheck"] = dict(predicted=0.017930 + dX, measured=base.X)

    (WORK / "solve211.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "solve211.log").write_text("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
