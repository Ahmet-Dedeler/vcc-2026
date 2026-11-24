"""187 — the rho-cliff amplitude model for `mse`, per context, under four C_p hypotheses.

Identity (096 (a), arithmetic, verified to 8.9e-16 against cell_eval2):

    MSE_c(s) - 1 = [ s^2 A_c + Kc_c - min(claim_c, s^2 B_c) - 2 s X ] / sumD        per context
    raw          = mean_c MSE_c ;  score_mse = mean_c clip((b - MSE_c)/span, 0, 1)

`A_c`, `B_c`, `Kc_c`, `C_hat_p` are ours (measured on the shipped artifact).  `C_p` — the
reference's own per-perturbation jackknife — is withheld, and 096 substitutes the delete-one
jackknife of one 400-cell NTC construct.  A12-ARTIFACT-AUDIT §4 named two corrections to that
proxy and computed a12 under each, BEFORE abc-knee existed:

    LEVEL   : real perturbation groups on h1-test read 1.0195x the 400-cell control draws
    SPREAD  : real C_p varies across perturbations with CV 0.052 (h1-test); min(.,C_p) is
              concave, so a spread C_p always LOWERS claim relative to the mean proxy

Known-answer rows (rule 10), all from the 096 record:
    (i)   the C_p-free LS solve on a8+a10+a11 must return X = 0.019240, sumD = 0.229141
    (ii)  a12 under {proxy, level, spread, both} must return 0.898551 / 0.890602 / 0.904073 /
          0.894727 (A12-ARTIFACT-AUDIT §4 table; spread cells were Monte Carlo, tol 1e-3)
    (iii) a23 -> abc-knee under `proxy` must return 180's +0.005025 scaled (raw 0.867566)

Then the same four hypotheses, with NO free parameter, are applied to
    a14 (a12 at alpha 0.982)              board raw 0.882775
    abc-knee (a23 at 0.98234/1.01291/1.04786)   board raw 0.874491  (delta +0.002138 vs a23)

    uv run python mse187.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import norm

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
W096 = ROOT / "experiments/096-mse-floor/work"
W180 = ROOT / "experiments/180-def-readiness/work"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments/096-mse-floor"))

CTX = ("A", "B", "C")
B_ANCHOR, SPAN = 0.987862879, 0.950162566        # 096 (d): score_mse = clip((b-raw)/span,0,1)
ALPHA_PROD = 1.0796

BOARD_RAW = {"a8": 1.362666960949, "a10": 1.3524148, "a11": 0.9099893, "a12": 0.9013492,
             "a14": 0.882775, "a15": 0.873805, "a23": 0.872353, "knee": 0.874491}
ALPHA = {"a8": 1.0796, "a10": 1.0796, "a11": 0.838, "a12": 1.0108, "a14": 0.982}
H1_LEVEL, H1_CV = 1.0195097955261299, 0.05210618886887627     # 096 cp_bias.npz

HYP = {"proxy": (1.0, 0.0), "level": (H1_LEVEL, 0.0), "spread": (1.0, H1_CV),
       "both": (H1_LEVEL, H1_CV)}


def score(raw_c):
    return float(np.mean([min(max((B_ANCHOR - r) / SPAN, 0.0), 1.0) for r in raw_c]))


def emin(chat, mu, sd):
    """E[min(chat, C)] for C ~ Normal(mu, sd), elementwise in chat.  sd=0 -> min(chat, mu)."""
    chat = np.asarray(chat, float)
    if sd <= 0:
        return np.minimum(chat, mu)
    d = (chat - mu) / sd
    return chat - sd * (norm.pdf(d) + d * norm.cdf(d))


def claim_of(chat_rows, weight_per_row, cp, bias, cv):
    """claim_c = sum_p w_p E[min(C_hat_p, C_p)] with C_p ~ N(bias*cp, cv*bias*cp)."""
    mu = bias * cp
    return float((emin(chat_rows, mu, cv * mu) * weight_per_row).sum())


def mse_ctx(A, B, Kc, claim, X_eff, sumD):
    """One context. A, B are at the arm's own amplitude; X_eff is the alignment at that
    amplitude (X_prod * s, or X_anchor * kappa)."""
    return 1.0 + (A + Kc - min(claim, B) - 2.0 * X_eff) / sumD


# ----------------------------------------------------------------------------------------------
# 1. the a8 / a10 / a11 / a12 artifacts, measured per context (096's own `measure`)
# ----------------------------------------------------------------------------------------------
def load_096_arms():
    from mech_a12_audit import cp_proxies, measure
    cps = cp_proxies()["ntc"]
    cp16 = {c: cps[c][0] for c in CTX}
    arms = {}
    for arm in ("a8", "a10", "a11", "a12"):
        rows = measure(str(W096 / f"{arm}-stats.npz"), cp16)
        per = {}
        for r in rows:
            G = r["G"]
            n = r["n"]
            # 096's row weights: w_p = 1/(G_p) normalised so that sum_p w_p * C_p == K; the
            # `measure` rows expose claim = sum(capped*w); recover w from K = sum(Cp*w).
            w_row = r["K"] / (r["Cp"] * n)          # uniform per row (all rows resolve target)
            per[r["ctx"]] = dict(A=r["A"], B=r["B"], Kc=r["Kc"], Chat=np.asarray(r["Chat"]),
                                 w=w_row, n=n, G=G, Cp=r["Cp"], claim096=r["claim"])
        arms[arm] = per
    return arms, cp16


def arm_terms(per, bias, cv, s=1.0, cp_override=None):
    """Per-context (A, B, Kc, claim) for an arm at amplitude factor s under a C_p hypothesis.
    Emission is unchanged by s (096: C_hat moves 0.01% over s in [0.15, 0.85]), so claim is
    computed at s=1 and A, B scale by s^2."""
    out = {}
    for c in CTX:
        p = per[c]
        cp = p["Cp"] if cp_override is None else cp_override[c]
        w = np.full(p["Chat"].size, p["w"])
        out[c] = dict(A=s * s * p["A"], B=s * s * p["B"], Kc=p["Kc"],
                      claim=claim_of(p["Chat"], w, cp, bias, cv))
    return out


def raw_of(terms, X, sumD, s=1.0, kappa=None):
    """terms carry A, B at the arm's amplitude; X is the anchor alignment and s (scalar) or
    kappa (per context) is the factor that carries it to this arm."""
    per = {}
    for c in CTX:
        k = s if kappa is None else kappa[c]
        per[c] = mse_ctx(terms[c]["A"], terms[c]["B"], terms[c]["Kc"], terms[c]["claim"],
                         k * X, sumD)
    return float(np.mean(list(per.values()))), per


def solve_ls(rows):
    """rows: list of (terms, s, raw).  (raw-1)*sumD + 2 s X = mean_c[A + Kc - min(claim,B)]."""
    M = np.array([[-2.0 * s, -(y - 1.0)] for _, s, y in rows])
    b = np.array([-np.mean([t[c]["A"] + t[c]["Kc"] - min(t[c]["claim"], t[c]["B"]) for c in CTX])
                  for t, _, _ in rows])
    sol, *_ = np.linalg.lstsq(M, b, rcond=None)
    return float(sol[0]), float(sol[1])


def solve_X(terms, raw, sumD, s=1.0):
    """One board point -> X, per-context branch, with s the amplitude factor of X's anchor."""
    num = np.mean([terms[c]["A"] + terms[c]["Kc"] - min(terms[c]["claim"], terms[c]["B"])
                   for c in CTX])
    return float((num - (raw - 1.0) * sumD) / (2.0 * s))


# ----------------------------------------------------------------------------------------------
# 2. a23 / abc-knee from 180's record (24 measured rows per context)
# ----------------------------------------------------------------------------------------------
def load_180():
    k = json.loads((W180 / "knee180.json").read_text())
    per = {}
    for c in CTX:
        b = k["per_context"][c]["before"]
        a = k["per_context"][c]["after"]
        per[c] = dict(
            A=b["A"], B=b["B"], Kc=k["Kc"][c], Chat=np.asarray(b["C_hat_values"]),
            Chat_after=np.asarray(a["C_hat_values"]), weight_sum=b["weight_sum"],
            Cp46=b["C_p"], Cp16=k["cp"][c]["cp_096_cached16"],
            scale=k["per_context"][c]["knee_scale"], kappa=k["geometry"][c]["kappa_mean"],
            A_after=a["A"], B_after=a["B"], claim180=b["claim"], claim180_after=a["claim"],
        )
    return per, k


def terms_180(per, bias, cv, which="before", cp_key="Cp46"):
    out = {}
    for c in CTX:
        p = per[c]
        chat = p["Chat"] if which == "before" else p["Chat_after"]
        w = np.full(chat.size, p["weight_sum"] / chat.size)
        out[c] = dict(A=p["A"] if which == "before" else p["A_after"],
                      B=p["B"] if which == "before" else p["B_after"], Kc=p["Kc"],
                      claim=claim_of(chat, w, p[cp_key], bias, cv))
    return out


def main():
    lines = []
    P = lambda *a: (print(*a, flush=True), lines.append(" ".join(str(x) for x in a)))

    arms, cp16 = load_096_arms()
    P("\n=== per-context terms of the four 096 artifacts (proxy C_p, 096's 16 constructs) ===")
    for arm, per in arms.items():
        for c in CTX:
            t = arm_terms(per, 1.0, 0.0)[c]
            P(f"  {arm:>4} {c}  A {t['A']:.6f}  B {t['B']:.6f}  claim {t['claim']:.6f} "
              f"(096 {per[c]['claim096']:.6f})  Kc {t['Kc']:.6f}  rho {min(1, t['B']/t['claim']):.4f}"
              f"  cap binds {(per[c]['Chat'] > per[c]['Cp']).sum()}/{per[c]['n']}")

    # ---- known answer (i): the C_p-free solve --------------------------------------------
    rows = [(arm_terms(arms[a], 1.0, 0.0), ALPHA[a] / ALPHA_PROD, BOARD_RAW[a])
            for a in ("a8", "a10", "a11")]
    X0, sumD0 = solve_ls(rows)
    P(f"\n[KA-i] C_p-free LS solve a8+a10+a11: X = {X0:.6f}  sumD = {sumD0:.6f}   "
      f"(096: 0.019240 / 0.229141)  {'PASS' if abs(X0-0.019240)<2e-5 and abs(sumD0-0.229141)<2e-4 else 'FAIL'}")
    for t, s, y in rows:
        r, _ = raw_of(t, X0, sumD0, s)
        P(f"       residual {r - y:+.5f}")

    # ---- known answer (ii): a12 under the four hypotheses --------------------------------
    KA_A12 = {"proxy": 0.898551, "level": 0.890602, "spread": 0.904073, "both": 0.894727}
    s12 = ALPHA["a12"] / ALPHA_PROD
    P("\n[KA-ii] a12 under the four C_p hypotheses (X, sumD from the C_p-free solve):")
    P(f"  {'hyp':>7}{'claim':>10}{'rho_min':>9}{'raw':>10}{'096 said':>10}{'diff':>9}{'vs board':>10}")
    res12 = {}
    for h, (bias, cv) in HYP.items():
        t = arm_terms(arms["a12"], bias, cv)
        r, per = raw_of(t, X0, sumD0, s12)
        rho_min = min(min(1.0, t[c]["B"] / t[c]["claim"]) for c in CTX)
        cl = np.mean([t[c]["claim"] for c in CTX])
        ok = abs(r - KA_A12[h]) < (1.5e-3 if cv > 0 else 2e-5)
        res12[h] = r
        P(f"  {h:>7}{cl:>10.6f}{rho_min:>9.4f}{r:>10.6f}{KA_A12[h]:>10.6f}{r-KA_A12[h]:>+9.5f}"
          f"{r-BOARD_RAW['a12']:>+10.5f}   {'PASS' if ok else 'FAIL'}")

    # ---- a14 = a12 at 0.982, same emission, no free parameter ----------------------------
    s14 = ALPHA["a14"] / ALPHA_PROD
    f14 = ALPHA["a14"] / ALPHA["a12"]
    P("\n[retrodict] a14 (a12 at alpha 0.982; A,B x s^2, C_hat unchanged), board 0.882775:")
    P(f"  {'hyp':>7}{'raw pred':>10}{'vs board':>10}{'ctx cliffs (alpha) of the a12 profile':>50}")
    res14 = {}
    f14 = ALPHA["a14"] / ALPHA["a12"]
    for h, (bias, cv) in HYP.items():
        t = arm_terms(arms["a12"], bias, cv, s=f14)
        r, per = raw_of(t, X0, sumD0, s14)
        cliffs = [ALPHA["a14"] * np.sqrt(t[c]["claim"] / t[c]["B"]) for c in CTX]
        res14[h] = r
        P(f"  {h:>7}{r:>10.6f}{r-BOARD_RAW['a14']:>+10.5f}    " + " ".join(f"{c}:{v:.4f}" for c, v in zip(CTX, cliffs)))
    P("  a12 -> a14 delta (board -0.018574):")
    for h in HYP:
        P(f"  {h:>7}  pred delta {res14[h]-res12[h]:+.6f}   err {(res14[h]-res12[h])-(BOARD_RAW['a14']-BOARD_RAW['a12']):+.6f}")

    # ---- a23 -> abc-knee -----------------------------------------------------------------
    per180, k180 = load_180()
    P("\n=== a23 (180's 24 measured rows per context) ===")
    for c in CTX:
        p = per180[c]
        P(f"  {c}: A {p['A']:.6f} B {p['B']:.6f} Kc {p['Kc']:.6f} C_hat mean {p['Chat'].mean():.4f} "
          f"C_p46 {p['Cp46']:.4f} C_p16 {p['Cp16']:.4f} claim180 {p['claim180']:.6f} "
          f"knee scale {p['scale']:.6f} kappa {p['kappa']:.6f}")
    # known answer (iii): 180's own numbers under `proxy` with the all-46 C_p
    t23 = terms_180(per180, 1.0, 0.0)
    for c in CTX:
        assert abs(t23[c]["claim"] - per180[c]["claim180"]) < 1e-9, (c, t23[c]["claim"], per180[c]["claim180"])
    X23 = solve_X(t23, BOARD_RAW["a23"], sumD0)
    scales = {c: per180[c]["scale"] for c in CTX}
    kappas = {c: per180[c]["kappa"] for c in CTX}
    tk = terms_180(per180, 1.0, 0.0, which="after")
    # 180 used A_after/B_after (rebuilt) and claim_after; amplitude enters through those, so s=1
    r_knee, per_knee = raw_of(tk, X23, sumD0, 1.0, kappa=kappas)
    r_a23, per_a23 = raw_of(t23, X23, sumD0, 1.0)
    d_scaled = score(list(per_knee.values())) - score(list(per_a23.values()))
    P(f"\n[KA-iii] proxy hypothesis: X23 {X23:.6f} (180: 0.017932); raw a23 {r_a23:.6f} -> knee "
      f"{r_knee:.6f}; d score_mse {d_scaled:+.6f} (180: +0.005025)  "
      f"{'PASS' if abs(X23-0.017932)<2e-5 and abs(d_scaled-0.005025)<5e-5 else 'FAIL'}")

    P("\n[retrodict] a23 -> abc-knee, board delta raw +0.002138 / scaled -0.002259:")
    P(f"  {'hyp':>7}{'X23':>10}{'a23 rho A/B/C':>26}{'knee alpha A/B/C (true cliff)':>34}"
      f"{'d raw':>10}{'d scaled':>10}{'err(scaled)':>12}")
    resk = {}
    for h, (bias, cv) in HYP.items():
        t23h = terms_180(per180, bias, cv)
        tkh = terms_180(per180, bias, cv, which="after")
        Xh = solve_X(t23h, BOARD_RAW["a23"], sumD0)
        r0, p0 = raw_of(t23h, Xh, sumD0, 1.0)
        r1, p1 = raw_of(tkh, Xh, sumD0, 1.0, kappa=kappas)
        ds = score(list(p1.values())) - score(list(p0.values()))
        rho = [min(1.0, t23h[c]["B"] / t23h[c]["claim"]) for c in CTX]
        cliff = [0.982 * np.sqrt(t23h[c]["claim"] / t23h[c]["B"]) for c in CTX]
        resk[h] = dict(X=Xh, d_raw=r1 - r0, d_scaled=ds, cliffs=cliff, rho=rho,
                       per_ctx_d_raw={c: p1[c] - p0[c] for c in CTX},
                       claim={c: t23h[c]["claim"] for c in CTX})
        P(f"  {h:>7}{Xh:>10.6f}   " + "/".join(f"{v:.4f}" for v in rho) + "      "
          + "/".join(f"{v:.4f}" for v in cliff) + f"{r1-r0:>+10.6f}{ds:>+10.6f}{ds-(-0.002259):>+12.6f}")
        P("           per-context d raw: " + "  ".join(f"{c} {p1[c]-p0[c]:+.6f}" for c in CTX))

    # ---- CV x level sweep: which cell fits a12, a14 and abc-knee together? ----------------
    P("\n[sweep] residuals (pred - board) for a12 raw, a14 raw, a12->a14 delta, abc-knee delta scaled:")
    P(f"  {'bias':>7}{'cv':>7}{'a12':>10}{'a14':>10}{'d(a12->a14)':>13}{'knee d(scaled)':>16}{'rms':>9}")
    sweep = []
    for bias in (0.98, 0.99, 1.0, 1.01, H1_LEVEL, 1.03):
        for cv in (0.0, 0.02, 0.03, 0.04, H1_CV, 0.06, 0.07, 0.08, 0.10):
            t = arm_terms(arms["a12"], bias, cv)
            r12, _ = raw_of(t, X0, sumD0, s12)
            r14, _ = raw_of(arm_terms(arms["a12"], bias, cv, s=f14), X0, sumD0, s14)
            t23h = terms_180(per180, bias, cv)
            tkh = terms_180(per180, bias, cv, which="after")
            Xh = solve_X(t23h, BOARD_RAW["a23"], sumD0)
            r0, p0 = raw_of(t23h, Xh, sumD0, 1.0)
            r1, p1 = raw_of(tkh, Xh, sumD0, 1.0, kappa=kappas)
            ds = score(list(p1.values())) - score(list(p0.values()))
            e12 = r12 - BOARD_RAW["a12"]
            e14 = r14 - BOARD_RAW["a14"]
            ed = (r14 - r12) - (BOARD_RAW["a14"] - BOARD_RAW["a12"])
            ek = ds - (-0.002259)
            rms = float(np.sqrt(np.mean([e12 ** 2, e14 ** 2, ek ** 2])))
            sweep.append(dict(bias=bias, cv=cv, e_a12=e12, e_a14=e14, e_delta=ed, e_knee=ek, rms=rms,
                              X23=Xh, cliffs=[0.982 * np.sqrt(t23h[c]["claim"] / t23h[c]["B"]) for c in CTX]))
            P(f"  {bias:>7.4f}{cv:>7.3f}{e12:>+10.5f}{e14:>+10.5f}{ed:>+13.5f}{ek:>+16.5f}{rms:>9.5f}")

    best = min(sweep, key=lambda d: d["rms"])
    P(f"\n  best cell: bias {best['bias']:.4f} cv {best['cv']:.3f}  rms {best['rms']:.5f}  "
      f"a23 cliffs (alpha) " + "/".join(f"{v:.4f}" for v in best["cliffs"]))

    out = dict(X0=X0, sumD0=sumD0, res12=res12, res14=res14, knee=resk, sweep=sweep, best=best,
               arms={a: {c: dict(A=arms[a][c]["A"], B=arms[a][c]["B"], Kc=arms[a][c]["Kc"],
                                 Cp=arms[a][c]["Cp"], w=arms[a][c]["w"], n=arms[a][c]["n"],
                                 Chat=arms[a][c]["Chat"].tolist()) for c in CTX} for a in arms},
               a23={c: {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in per180[c].items()}
                    for c in CTX})
    (HERE / "work").mkdir(exist_ok=True)
    (HERE / "work/mse187.json").write_text(json.dumps(out, indent=1, default=float))
    (HERE / "work/mse187.log").write_text("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
