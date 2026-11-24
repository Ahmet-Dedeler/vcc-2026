"""187 — the joint (mse + nmae) amplitude optimum on a23's profile, global and per context,
under the retrodicting cliff model, with the composed slot gate.

mse   : per-context rho-cliff identity (mse187 / cliff187), C_p = f_c x NTC proxy with
        f = (1.00, 1.00, 0.9825) — the only cell that reproduces BOTH a12->a14 and a23->abc-knee.
        X23 solved from a23's own board raw (one equation).  raw -> scaled through 096's map,
        per context, clamp never active here (all raws ~0.87 < 0.988).
nmae  : board-measured local law on the current pipeline,
            raw_nmae(alpha_c) = 0.944330 - SLOPE * mean_c ln(alpha_c / 0.982),  SLOPE = 0.0416
        (a12->a14: 0.0471, a23->knee: 0.0361), scaled = -d raw / 0.5963 (our own board map).
        Beyond +-10% the h1-test shape (nmae187: turnover at s~1.29) is used as a saturation
        prior; every optimum below lies inside +-6%, where the linear law is what abc-knee measured.
overall = (d score_mse + d score_nmae) / 6 ; pds is scale-invariant (180: PLOLO -0.000076 raw),
fid / jac / reach are carried at 0 with 180's twin noise as the band.

    uv run python joint187.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import norm

from cliff187 import P23, KAPPA, terms_23, X0, sumD0
from mse187 import CTX, BOARD_RAW, score, solve_X, mse_ctx

HERE = Path(__file__).resolve().parent
ALPHA0 = 0.982
F_FIT = {"A": 1.0, "B": 1.0, "C": 0.9825}
F_PROXY = {"A": 1.0, "B": 1.0, "C": 1.0}
ZERO = {c: 0.0 for c in CTX}
NMAE_SLOPE = 0.0416           # raw per ln(alpha), current pipeline (mean of the two board pairs)
NMAE_SPAN = 0.596294          # our board map: scaled = (b - raw)/span
NMAE_A23 = 0.944330
A23_OVERALL, KNEE_OVERALL = 0.1733597079, 0.1740437
KNEE_ALPHA = {"A": 0.98234, "B": 1.01291, "C": 1.04786}


class Model:
    def __init__(self, f):
        self.t = terms_23(f, ZERO)                      # a23 terms at alpha 0.982
        self.X = solve_X(self.t, BOARD_RAW["a23"], sumD0)
        self.cliff = {c: ALPHA0 * float(np.sqrt(self.t[c]["claim"] / self.t[c]["B"])) for c in CTX}
        self.raw0 = {c: mse_ctx(self.t[c]["A"], self.t[c]["B"], self.t[c]["Kc"], self.t[c]["claim"],
                                self.X, sumD0) for c in CTX}
        self.score0 = score(list(self.raw0.values()))

    def raw_c(self, c, alpha):
        s = alpha / ALPHA0
        t = self.t[c]
        return mse_ctx(s * s * t["A"], s * s * t["B"], t["Kc"], t["claim"], s * self.X, sumD0)

    def d_mse_scaled(self, alphas):
        return score([self.raw_c(c, alphas[c]) for c in CTX]) - self.score0

    def d_raw_mean(self, alphas):
        return float(np.mean([self.raw_c(c, alphas[c]) - self.raw0[c] for c in CTX]))


def d_nmae_scaled(alphas):
    dln = float(np.mean([np.log(alphas[c] / ALPHA0) for c in CTX]))
    d_raw = -NMAE_SLOPE * dln
    return -d_raw / NMAE_SPAN, d_raw


def evaluate(m, alphas):
    dm = m.d_mse_scaled(alphas)
    dn, dn_raw = d_nmae_scaled(alphas)
    return dict(d_mse=dm, d_nmae=dn, d_nmae_raw=dn_raw, d_raw_mse=m.d_raw_mean(alphas),
                overall=(dm + dn) / 6.0)


def main():
    lines = []
    P = lambda *a: (print(*a, flush=True), lines.append(" ".join(str(x) for x in a)))
    out = {}

    models = {"fit f=(1,1,.9825)": Model(F_FIT), "proxy f=(1,1,1)": Model(F_PROXY)}
    for name, m in models.items():
        P(f"\n=== {name}: X23 {m.X:.6f}; a23 cliffs (alpha) " + " ".join(f"{c} {m.cliff[c]:.4f}" for c in CTX)
          + f"; a23 raw per ctx " + " ".join(f"{c} {m.raw0[c]:.5f}" for c in CTX) + f"; score {m.score0:.6f}")

        # known-answer: abc-knee's own alphas
        e = evaluate(m, KNEE_ALPHA)
        P(f"  [KA] abc-knee alphas -> d mse {e['d_mse']:+.6f} scaled (board -0.002259), d nmae {e['d_nmae']:+.6f} "
          f"(board +0.002081), overall {e['overall']:+.6f} (board mse+nmae: -0.000030)")

        # global alpha sweep
        P("  global alpha sweep (deltas vs a23):")
        P(f"  {'alpha':>7}{'d raw mse':>11}{'d mse sc':>10}{'d nmae sc':>10}{'overall':>10}   above-cliff contexts")
        best_g = None
        for a in np.round(np.arange(0.90, 1.1001, 0.005), 4):
            al = {c: float(a) for c in CTX}
            e = evaluate(m, al)
            above = [c for c in CTX if a > m.cliff[c] + 1e-9]
            if best_g is None or e["overall"] > best_g[1]["overall"]:
                best_g = (float(a), e)
            P(f"  {a:>7.3f}{e['d_raw_mse']:>+11.5f}{e['d_mse']:>+10.5f}{e['d_nmae']:>+10.5f}{e['overall']:>+10.5f}   {''.join(above) or '-'}")
        # fine around the top
        fine = np.arange(0.960, 1.0201, 0.0005)
        vals = [(float(a), evaluate(m, {c: float(a) for c in CTX})["overall"]) for a in fine]
        ag, vg = max(vals, key=lambda x: x[1])
        P(f"  global optimum: alpha {ag:.4f}  overall {vg:+.6f} vs a23  ({vg - (KNEE_OVERALL - A23_OVERALL):+.6f} vs abc-knee)")
        out[name] = dict(X=m.X, cliffs=m.cliff, global_opt=ag, global_gain=vg)

        # per-context: each at its own cliff, and with margins
        P("  per-context arms (deltas vs a23; abc-knee realised +0.000684 overall, mse -0.002259, nmae +0.002081):")
        P(f"  {'arm':>28}{'A':>8}{'B':>8}{'C':>8}{'d raw mse':>11}{'d mse sc':>10}{'d nmae sc':>10}{'overall':>10}{'vs knee':>10}")
        arms = {"on the cliffs": {c: m.cliff[c] for c in CTX}}
        for margin in (0.005, 0.01, 0.02, 0.03):
            arms[f"cliffs x (1-{margin:.3f})"] = {c: m.cliff[c] * (1 - margin) for c in CTX}
        arms["A a23, B/C on cliffs"] = {"A": ALPHA0, "B": m.cliff["B"], "C": m.cliff["C"]}
        arms["abc-knee as shipped"] = dict(KNEE_ALPHA)
        arms["a23 (0.982 global)"] = {c: ALPHA0 for c in CTX}
        res = {}
        for nm, al in arms.items():
            e = evaluate(m, al)
            res[nm] = dict(alphas=al, **e)
            P(f"  {nm:>28}{al['A']:>8.4f}{al['B']:>8.4f}{al['C']:>8.4f}{e['d_raw_mse']:>+11.5f}{e['d_mse']:>+10.5f}"
              f"{e['d_nmae']:>+10.5f}{e['overall']:>+10.5f}{e['overall'] - (KNEE_OVERALL - A23_OVERALL):>+10.5f}")
        out[name]["arms"] = res

    # ---- expected value under cliff uncertainty (fit model) ---------------------------------
    m = models["fit f=(1,1,.9825)"]
    P("\n=== expected overall (mse+nmae) under cliff-position uncertainty, fit model ===")
    P("  cliff of context c ~ Normal(cliff_c, sd) in ln(alpha); sd from the retrodiction region "
      "(A: f_A in [0.99,1.00] -> cliff [0.977,0.982]; B: [1.0123,1.0131]; C: [1.0373,1.0400])")
    sd_fit = {"A": 0.0025, "B": 0.0004, "C": 0.0013}
    sd_blind = {c: 0.010 for c in CTX}     # D/E/F: f_c unknown, |f-1| ~ 2% -> cliff sd ~1%
    for label, sd in (("A/B/C, retrodicted", sd_fit), ("blind (D/E/F-like), sd 1% per ctx", sd_blind)):
        P(f"  --- {label} ---")
        P(f"  {'ctx':>4}{'cliff':>8}{'best alpha':>11}{'margin':>9}{'E[d overall]':>14}{'at cliff':>10}{'at -1%':>9}{'at -2%':>9}")
        tot = 0.0
        best_alphas = {}
        for c in CTX:
            cl = m.cliff[c]
            zs = np.linspace(-3.5, 3.5, 141)
            wts = norm.pdf(zs); wts /= wts.sum()
            def ev(alpha):
                # expected per-context mse delta (raw -> scaled/3) + nmae share, marginalising the cliff
                acc = 0.0
                for z, w in zip(zs, wts):
                    cliff = cl * np.exp(z * sd[c])
                    # rebuild the context's claim so that its cliff sits at `cliff`
                    t = m.t[c]
                    claim = (cliff / ALPHA0) ** 2 * t["B"]
                    s = alpha / ALPHA0
                    r = mse_ctx(s * s * t["A"], s * s * t["B"], t["Kc"], claim, s * m.X, sumD0)
                    r0 = mse_ctx(t["A"], t["B"], t["Kc"], claim, m.X, sumD0)
                    acc += w * (-(r - r0) / 0.950162566) / 3.0
                dn = (NMAE_SLOPE * np.log(alpha / ALPHA0) / 3.0) / NMAE_SPAN
                return (acc + dn) / 6.0
            grid = cl * np.exp(np.linspace(-0.06, 0.02, 161))
            e = np.array([ev(a) for a in grid])
            i = int(e.argmax())
            best_alphas[c] = float(grid[i])
            tot += float(e[i])
            P(f"  {c:>4}{cl:>8.4f}{grid[i]:>11.4f}{grid[i]/cl-1:>+9.3%}{e[i]:>+14.6f}{ev(cl):>+10.6f}{ev(cl*0.99):>+9.6f}{ev(cl*0.98):>+9.6f}")
        P(f"  sum of per-context expected gains vs a23: {tot:+.6f} overall   (vs abc-knee: {tot - (KNEE_OVERALL - A23_OVERALL):+.6f})")
        out[f"ev_{label}"] = dict(best=best_alphas, total=tot)

    # ---- gate --------------------------------------------------------------------------------
    P("\n=== composed slot gate (rule 13, audit-181 amendment: profile-class arm, no emission toll) ===")
    arm = models["fit f=(1,1,.9825)"].cliff
    e = evaluate(models["fit f=(1,1,.9825)"], arm)
    pred = (-0.000076 / 6.0) * 1.011 + e["d_mse"] / 6.0
    P(f"  arm: per-context alphas on the retrodicted cliffs " + " ".join(f"{c} {arm[c]:.4f}" for c in CTX))
    P(f"  PLOLO d_pds/6 x 1.011 = {(-0.000076/6)*1.011:+.6f}   096 d_mse/6 = {e['d_mse']/6:+.6f}   "
      f"(nmae, not in the gate formula, would add {e['d_nmae']/6:+.6f})")
    P(f"  predicted = {pred:+.6f}   floor (amplitude class) +0.0089   -> {'SHIP' if pred > 0.0089 else 'VETO'}")
    P(f"  even with nmae counted: {pred + e['d_nmae']/6:+.6f}; even vs a23 with everything optimistic it is "
      f"{(pred + e['d_nmae']/6)/0.0089:.1%} of the floor")
    out["gate"] = dict(arm=arm, predicted=pred, with_nmae=pred + e["d_nmae"] / 6.0)

    (HERE / "work/joint187.json").write_text(json.dumps(out, indent=1, default=float))
    (HERE / "work/joint187.log").write_text("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
