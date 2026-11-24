"""212 Phase 2/3 — the gamma sweep priced with COSINE primary, at the RE-SOLVED alpha.

THE ONE THING 170 DID NOT DO
----------------------------
170 swept gamma with `VCC_ALPHA` pinned at 0.982 while gamma moved the predicted effect's
mean square `A` by 14%. `mse`'s numerator is

    num_c(s) = s^2 A_c + Kc_c - min(claim_c, s^2 B_c) - 2 s X_c

so a gamma move at fixed `s` is a gamma move CONFOUNDED with an amplitude move, and above
the rho cliff the metric refuses to refund the amplitude. That is exactly what 170 measured
and correctly reported as a V with its vertex on 0.65.

But the transfer cosine

    cos = X / ( sqrt(SumD) * mean_c sqrt(A_c) )

is EXACTLY invariant to a global rescale, and at the constrained optimum in `s` the
amplitude terms cancel: the minimiser sits on the FIRST context's rho cliff,
`s* = min_c sqrt(claim_c / B_c)`, all three contexts are then on the `B` leg of the `min`,
and

    raw*(gamma) = 1 + [ s*^2 * mean_c(A_c - B_c) + mean_c Kc_c - 2 s* X ] / SumD

with `s* X = sqrt(claim_A/B_A) * cos * sqrt(SumD) * mean_c sqrt(A_c)`. Since `A_c` and `B_c`
move very nearly proportionally in gamma (measured below), `s* X` is proportional to
`sqrt(claim_A) * cos` and the whole thing collapses to

    raw*(gamma) = const - C * cos(gamma) + (tiny, favourable) ,      C > 0.

**At the re-solved alpha, `mse` is a monotone decreasing function of the transfer cosine and
of nothing else.** So the question 170 answered ("does gamma pay at alpha 0.982?") and the
question this experiment asks ("does gamma pay at its own alpha?") are different questions,
and this script computes both.

WHAT IS TRANSFERRED, AND HOW
----------------------------
`A_c`, `B_c`, `claim_c`, `Kc_c` are MEASURED on the challenge panel at every gamma
(`chal212.py`). `SumD` is 096's board-solved constant. `X` is solved from a15's own board
point at gamma 0.65 -- one equation, one unknown, as 096/170/197 all do. Its gamma
dependence is transferred by the COSINE ratio, not by `X`'s own ratio:

    X(gamma) = cos(0.65) * r(gamma) * sqrt(SumD) * mean_c sqrt(A_c(gamma))

because the cosine is the surface-transferable half and the amplitude belongs to whichever
surface it was measured on. 170 carried `X`'s ratio instead and let the challenge panel's own
`A` move underneath it, which divides h1-test's amplitude growth out by the challenge
panel's -- a different and, for a direction question, wrong contraction. Three independent
`r(gamma)` are carried side by side:

    plolo_a15   PLOLO, 4 folds, a15 render (the incumbent's geometry)   <- PRIMARY
    plolo_a23   PLOLO, 4 folds, a23 render
    h1          170's own h1-test `cos_deconf` = X / sqrt(A * SumD)
    x_ratio_170 170's route, for comparability with its published table

BOARD SCALING (docs/BOARD-LEDGER.md, fitted on 612 live entries):
    scaled = slope * (raw - b),  pds +2.2138, mse -1.051, nmae -1.6532.
PLOLO -> board `pds` uses 162's fitted CALIB 1.011.
`nmae`'s amplitude response is 187's board-measured `d nmae_raw / d ln alpha = -0.0416`.

    uv run python alpha212.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"

CTX = ("A", "B", "C")
SUMD = 0.229141
MSE_MAP = (0.987862879, 0.950162566)
A15_SCORE_MSE = 0.1200566151966586
SLOPE = dict(pds=2.2138, mse=-1.051, nmae=-1.6532)
PLOLO_CALIB = 1.011
DNMAE_DLNA = -0.0416                  # 187, board-measured
COS_197 = 0.079249                    # 197's measured a15/rc_hybrid transfer cosine
GATE_SHAPE, GATE_AMP = 0.0019, 0.0089


def score_mse(raw):
    return float(np.clip((MSE_MAP[0] - raw) / MSE_MAP[1], 0.0, 1.0))


def raw_at(per, X, s):
    """096's identity, per context, with `min(claim, s^2 B)` evaluated AT EVERY s."""
    acc = []
    for c in CTX:
        p = per[c]
        acc.append(s * s * p["A"] + p["Kc"] - min(p["claim"], s * s * p["B"]) - 2 * s * X)
    return 1.0 + float(np.mean(acc)) / SUMD


def solve_s(per, X):
    """Minimise `raw(s)` exactly. The minimum is at a kink or at the smooth interior root."""
    kinks = sorted(np.sqrt(p["claim"] / p["B"]) for p in (per[c] for c in CTX))
    cand = list(kinks) + [X / float(np.mean([per[c]["A"] for c in CTX]))]
    # dense scan to make the kink logic auditable rather than assumed
    grid = np.unique(np.concatenate([np.linspace(0.60, 1.40, 8001), np.asarray(cand)]))
    vals = np.asarray([raw_at(per, X, float(s)) for s in grid])
    i = int(np.argmin(vals))
    return float(grid[i]), float(vals[i]), [float(k) for k in kinks]


def main() -> None:
    chal = json.loads((WORK / "chal212_a15.json").read_text())
    chal = {round(r["gamma"], 3): r for r in chal}
    plolo = json.loads((WORK / "plolo212.json").read_text())
    h1 = {round(r["gamma"], 3): r for r in
          json.loads((REPO / "experiments/170-gamma-mse/work/h1sweep.json").read_text())}
    new = WORK / "h1sweep212.json"                # the five gammas 170 did not run
    if new.exists():
        h1.update({round(r["gamma"], 3): r for r in json.loads(new.read_text())})
    folds = sorted({r["fold"] for r in plolo})
    P = {(r["fold"], round(r["gamma"], 3)): r for r in plolo}

    gammas = sorted(chal)
    base = chal[0.65]

    # ---- X solved from a15's own board point -------------------------------------------
    target_raw = MSE_MAP[0] - A15_SCORE_MSE * MSE_MAP[1]
    num = float(np.mean([base["per"][c]["A"] + base["per"][c]["Kc"]
                         - min(base["per"][c]["claim"], base["per"][c]["B"]) for c in CTX]))
    X0 = (num - SUMD * (target_raw - 1.0)) / 2.0
    cos0 = X0 / (np.sqrt(SUMD) * base["mean_sqrtA"])
    print(f"[alpha] a15 board raw mse {target_raw:.6f} -> X(0.65) {X0:.6f}, "
          f"implied cos {cos0:.6f}  (197 measured 0.079249, rel "
          f"{abs(cos0/COS_197-1):.2%})")

    # ---- the four cosine-ratio routes ---------------------------------------------------
    def macro(g, key):
        return float(np.mean([P[(f, g)][key] for f in folds if (f, g) in P]))

    ratios = {}
    for g in gammas:
        r = {}
        for name, key in (("plolo_a15", "a15_cos"), ("plolo_a23", "a23_cos")):
            if all((f, g) in P for f in folds):
                r[name] = macro(g, key) / macro(0.65, key)
        gg = round(g, 3)
        if gg in h1 and 0.65 in h1:
            r["h1"] = h1[gg]["cos_deconf"] / h1[0.65]["cos_deconf"]
            # 170's own route: carry X's ratio, let the challenge panel's own A move
            r["x_ratio_170"] = ((h1[gg]["X"] / h1[0.65]["X"])
                                * base["mean_sqrtA"] / chal[g]["mean_sqrtA"])
        ratios[g] = r

    routes = ["plolo_a15", "plolo_a23", "h1", "x_ratio_170"]
    rows = []
    for g in gammas:
        per = chal[g]["per"]
        row = dict(gamma=g, A=chal[g]["A"], B=chal[g]["B"], claim=chal[g]["claim"],
                   mean_sqrtA=chal[g]["mean_sqrtA"],
                   rho={c: per[c]["rho"] for c in CTX},
                   s_cliff={c: float(np.sqrt(per[c]["claim"] / per[c]["B"])) for c in CTX},
                   chat_over_cp={c: per[c]["chat_over_cp"] for c in CTX})
        for route in routes:
            if route not in ratios[g]:
                continue
            r = ratios[g][route]
            cos = cos0 * r
            X = cos * np.sqrt(SUMD) * chal[g]["mean_sqrtA"]
            fixed = raw_at(per, X, 1.0)
            s_star, star, kinks = solve_s(per, X)
            row[route] = dict(
                ratio=r, cos=float(cos),
                raw_fixed_alpha=fixed, raw_resolved=star,
                s_star=s_star, alpha_star=0.982 * s_star, kinks=kinks,
                d_overall_mse_fixed=(score_mse(fixed)
                                     - score_mse(raw_at(base["per"], X0, 1.0))) / 6.0,
            )
        rows.append(row)

    # the re-solved baseline is the SAME arm at gamma 0.65 with its own alpha re-solved,
    # so the comparison is like-for-like on both sides.
    b_s, b_star, _ = solve_s(base["per"], X0)
    b_fixed = raw_at(base["per"], X0, 1.0)
    print(f"[alpha] gamma 0.65: raw at alpha 0.982 {b_fixed:.6f} (board 0.873805); "
          f"raw at the re-solved alpha {b_star:.6f} at s* {b_s:.5f} "
          f"(alpha* {0.982*b_s:.5f}); the alpha headroom alone is worth "
          f"{(score_mse(b_star)-score_mse(b_fixed))/6:+.5f} overall")

    for row in rows:
        for route in routes:
            if route in row:
                row[route]["d_overall_mse_resolved"] = (
                    score_mse(row[route]["raw_resolved"]) - score_mse(b_star)) / 6.0

    # ---- pds, nmae ----------------------------------------------------------------------
    for row in rows:
        g = row["gamma"]
        if all((f, g) in P for f in folds):
            dp = macro(g, "a23_pds") - macro(0.65, "a23_pds")
            dp15 = macro(g, "a15_pds") - macro(0.65, "a15_pds")
            row["d_pds_plolo_raw"] = dp
            row["d_overall_pds"] = dp * PLOLO_CALIB * SLOPE["pds"] / 6.0
            row["d_overall_pds_a15render"] = dp15 * PLOLO_CALIB * SLOPE["pds"] / 6.0
            row["q_w"] = macro(g, "a15_sgn_w")
        gg = round(g, 3)
        if gg in h1:
            dn = h1[gg]["nmae"] - h1[0.65]["nmae"]
            row["d_nmae_h1_raw"] = dn
            row["d_overall_nmae_fixed"] = dn * SLOPE["nmae"] / 6.0
            row["d_pds_h1_raw"] = h1[gg]["pds"] - h1[0.65]["pds"]

    with open(WORK / "alpha212.json", "w") as fh:
        json.dump(dict(X0=X0, cos0=float(cos0), base_raw_fixed=b_fixed,
                       base_raw_resolved=b_star, base_s_star=b_s, rows=rows), fh, indent=1)

    # ---- report -------------------------------------------------------------------------
    print("\n=== the two cliffs, per context, on the current profile ===")
    print("gamma |  rho A/B/C           | s_cliff A/B/C        | Chat/C_p A/B/C")
    for row in rows:
        print(f"{row['gamma']:5.3f} | "
              + "/".join(f"{row['rho'][c]:.4f}" for c in CTX) + " | "
              + "/".join(f"{row['s_cliff'][c]:.4f}" for c in CTX) + " | "
              + "/".join(f"{row['chat_over_cp'][c]:.3f}" for c in CTX))

    # ---- nmae at the RE-SOLVED alpha: 187's board slope d nmae_raw / d ln alpha ---------
    for row in rows:
        r = row.get("plolo_a15")
        if r is None or "d_overall_nmae_fixed" not in row:
            continue
        dln = float(np.log(r["s_star"] / b_s))
        row["d_lnalpha"] = dln
        row["d_overall_nmae_alpha"] = DNMAE_DLNA * dln * SLOPE["nmae"] / 6.0
        row["d_overall_nmae_resolved"] = (row["d_overall_nmae_fixed"]
                                          + row["d_overall_nmae_alpha"])

    print("\n=== gamma priced with COSINE primary (route plolo_a15) ===")
    print("gamma |   cos    ratio |  raw@0.982  d_ovr | alpha*  raw*   d_ovr* | "
          "d_pds   d_nmae* | TOTAL* | rule 13")
    for row in rows:
        r = row.get("plolo_a15")
        if r is None:
            continue
        dp = row.get("d_overall_pds", float("nan"))
        dn = row.get("d_overall_nmae_resolved", float("nan"))
        tot = r["d_overall_mse_resolved"] + (dp if dp == dp else 0) + (dn if dn == dn else 0)
        row["total_resolved"] = tot
        verdict = "SHIP" if tot > GATE_AMP else "VETO"
        print(f"{row['gamma']:5.3f} | {r['cos']:.6f} {r['ratio']:6.4f} | "
              f"{r['raw_fixed_alpha']:9.6f} {r['d_overall_mse_fixed']:+.5f} | "
              f"{r['alpha_star']:6.4f} {r['raw_resolved']:.6f} "
              f"{r['d_overall_mse_resolved']:+.5f} | {dp:+.5f} {dn:+.5f} | {tot:+.5f} | "
              f"{verdict} (floor +{GATE_AMP:.4f})")

    print("\n=== transfer-route sensitivity: d_overall on `mse` at the re-solved alpha ===")
    print("gamma | " + " | ".join(f"{r:>14}" for r in routes))
    for row in rows:
        cells = []
        for route in routes:
            r = row.get(route)
            cells.append("            --" if r is None
                         else f"{r['ratio']:.4f} {r['d_overall_mse_resolved']:+.5f}")
        print(f"{row['gamma']:5.3f} | " + " | ".join(cells))

    print("\n=== the marginal exchange, both members per unit of BOARD cosine ===")
    g = [r for r in rows if r["gamma"] in (0.60, 0.70)]
    if len(g) == 2:
        lo, hi = g
        dcos = hi["plolo_a15"]["cos"] - lo["plolo_a15"]["cos"]
        dmse = (hi["plolo_a15"]["d_overall_mse_resolved"]
                - lo["plolo_a15"]["d_overall_mse_resolved"])
        dpds = hi["d_overall_pds"] - lo["d_overall_pds"]
        print(f"  d(overall mse*)/d(cos) = {dmse/dcos:+.4f}    "
              f"d(overall pds)/d(cos) = {dpds/dcos:+.4f}    "
              f"sum = {(dmse+dpds)/dcos:+.4f}   "
              f"(197's published mse exchange: +0.3655 per unit cos)")
        print(f"  they cancel to {abs((dmse+dpds)/dcos)/abs(dmse/dcos):.1%} of the mse leg, "
              f"across gamma 0.60-0.70")

    # ---- the same question on 210's `cd4f034` delivered profile --------------------------
    cd4 = WORK / "cd4f212_summary.json"
    if cd4.exists():
        C = {round(r["gamma"], 3): r for r in json.loads(cd4.read_text())}
        rate = (rows[[r["gamma"] for r in rows].index(0.70)]["plolo_a15"]
                ["d_overall_mse_resolved"]
                - rows[[r["gamma"] for r in rows].index(0.60)]["plolo_a15"]
                ["d_overall_mse_resolved"])
        rate /= (rows[[r["gamma"] for r in rows].index(0.70)]["plolo_a15"]["cos"]
                 - rows[[r["gamma"] for r in rows].index(0.60)]["plolo_a15"]["cos"])
        lift = C[0.65]["cos_f034"] / C[0.65]["cos_f010"]
        base34 = cos0 * lift
        print(f"\n=== 210's `cd4f034` (VCC_CD4_TARGET_F 0.34): the same sweep ===")
        print(f"  the dose raises the absolute PLOLO cosine {lift:.4f}x, so the board "
              f"cosine goes {cos0:.6f} -> {base34:.6f} (197 measured +6.8% for the f=0.10 "
              f"residual alone; 210 measured delivered A at 0.97970x, which moves context "
              f"A's rho cliff from s {b_s:.4f} to {b_s/np.sqrt(0.97970):.4f})")
        print("  `mse` is priced through the exchange rate d(overall mse*)/d(cos) = "
              f"{rate:+.4f}, which is EXACT at the re-solved alpha because the amplitude "
              "cancels there; the challenge-panel A/B/claim at f=0.34 were not re-measured.")
        print("gamma | cos ratio | d_ovr mse* |  d_pds   | d_nmae* |  TOTAL* | "
              "  vs f=0.10")
        for row in rows:
            g = round(row["gamma"], 3)
            if g not in C:
                continue
            dm = rate * base34 * (C[g]["ratio_f034"] - 1.0)
            dp = C[g]["d_pds_f034"] * PLOLO_CALIB * SLOPE["pds"] / 6.0
            dn = row.get("d_overall_nmae_resolved", float("nan"))
            tot = dm + dp + (dn if dn == dn else 0.0)
            row["cd4f034"] = dict(ratio=C[g]["ratio_f034"], d_mse=dm, d_pds=dp, total=tot)
            print(f"{g:5.3f} | {C[g]['ratio_f034']:9.4f} | {dm:+10.5f} | {dp:+.5f} | "
                  f"{dn:+.5f} | {tot:+.5f} | {row.get('total_resolved', float('nan')):+.5f}")
        best = max((r for r in rows if "cd4f034" in r), key=lambda r: r["cd4f034"]["total"])
        print(f"  best on cd4f034: gamma {best['gamma']:.3f} at "
              f"{best['cd4f034']['total']:+.5f} overall, against rule 13's "
              f"+{GATE_SHAPE:.4f} (shape) / +{GATE_AMP:.4f} (amplitude) floors")
        with open(WORK / "alpha212.json", "w") as fh:
            json.dump(dict(X0=X0, cos0=float(cos0), base_raw_fixed=b_fixed,
                           base_raw_resolved=b_star, base_s_star=b_s,
                           exchange_rate_mse_per_cos=float(rate), rows=rows), fh, indent=1)

    print("\n=== what the whole operator costs on cosine, and what the gap needs ===")
    z = [r for r in rows if r["gamma"] == 0.00][0]["plolo_a15"]
    print(f"  gamma 0.65 -> 0.00 (no whitening at all): cos {rows[-4]['plolo_a15']['cos']:.6f}"
          f" -> {z['cos']:.6f}, ratio {z['ratio']:.4f}x")
    print(f"  the gap 197 measured: leaders >= 0.13117 = {0.13117/COS_197:.3f}x ours; "
          f"overall leader 0.12821 = {0.12821/COS_197:.3f}x")
    print(f"  so abolishing the operator reaches {z['ratio']*COS_197:.6f}, still "
          f"{0.13117/(z['ratio']*COS_197):.3f}x short of the board's best direction, "
          f"and it costs {[r for r in rows if r['gamma']==0.0][0]['d_pds_plolo_raw']:+.5f} "
          f"raw pds = "
          f"{[r for r in rows if r['gamma']==0.0][0]['d_overall_pds']:+.5f} overall.")


if __name__ == "__main__":
    main()
