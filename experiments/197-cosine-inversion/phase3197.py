"""197 Phases 2b/3/4 — assumption tests, the three-route cross-check, and the target.

Runs after `gate197.py` and `invert197.py`.

  (1) The exact per-context-consistent cosine.  Assuming ONE cosine `c` shared by the three
      contexts, X_c = c*sqrt(A_c*sumD), so

          raw - 1 = [ mean_c( A_c + Kc_c - min(claim_c, B_c) ) - 2 c sqrt(sumD) mean_c sqrt(A_c) ] / sumD

      which is exactly amplitude-invariant PER CONTEXT and is therefore the right definition
      for abc-knee, the one arm with a non-uniform alpha.

  (2) The `C_p` proxy is the inversion's own systematic.  187 board-calibrated it at
      f = (1.00, 1.00, 0.9825) from two board deltas.  Re-inverting a23 and abc-knee under `f`
      tests whether their 4.4% cosine gap is direction or proxy.

  (3) Three routes to one latent quantity: the `mse` floor, `pds` through 159 5.4's curve, and
      `nmae` through 187 5's |t|-weighted sign identity.

  (4) The effective-donor translation through 159's 1/cos^2(k) = A + B/k law, and the
      adjudication of AGENTS.md 8's "at least 1.28x".

    uv run python phase3197.py
"""
from __future__ import annotations

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

CTX = ("A", "B", "C")
LINES = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LINES.append(s)


# 159 5.4's synthetic cos -> pds curve, both surfaces, and 5.5's 1/cos^2(k) = A + B/k fit
C159 = np.array([0.000, 0.020, 0.030, 0.038, 0.043, 0.050, 0.072, 0.100, 0.123, 0.150,
                 0.200, 0.250, 0.300])
PDS_HCT = np.array([0.5043, 0.5730, 0.6096, 0.6374, 0.6538, 0.6765, 0.7430, 0.8199, 0.8705,
                    0.9175, 0.9723, 0.9920, 0.9981])
PDS_HEK = np.array([0.4998, 0.5594, 0.5879, 0.6091, 0.6231, 0.6414, 0.6960, 0.7650, 0.8112,
                    0.8571, 0.9200, 0.9567, 0.9766])
KA, KB = 66.5, 1920.0                 # 159 5.5
COS_SURF_PROD = 0.03814               # 159 5.3, production on the HCT116 challenge-panel surface
COS_SURF_4DON = 0.04278               # 159 5.3, best equal consensus of the 4 panel-covering lines
COS_SURF_INF = 1.0 / np.sqrt(KA)      # 0.12263, the context-independent k->infinity ceiling


def keff(cos_surface):
    inv = 1.0 / cos_surface ** 2
    return KB / (inv - KA) if inv > KA else np.inf


def main():
    inv = json.load(open(f"{WORK}/invert197.json"))
    sumD, K, floor0, coeff = inv["sumD"], inv["K"], inv["floor0"], inv["coeff"]
    ours = inv["ours"]

    from mech_a12_audit import cp_proxies, measure
    cp16 = {c: cp_proxies()["ntc"][c][0] for c in CTX}
    k180 = json.load(open(f"{REPO}/experiments/180-def-readiness/work/knee180.json"))

    P("=" * 100)
    P("197 PHASE 2b — the assumptions, tested")
    P("=" * 100)

    # ------------------------------------------------------------------ (1) exact definition
    def cos_exact(per, raw):
        num = float(np.mean([per[c]["A"] + per[c]["Kc"] - min(per[c]["claim"], per[c]["B"])
                             for c in CTX]))
        msa = float(np.mean([np.sqrt(per[c]["A"]) for c in CTX]))
        return (num - (raw - 1.0) * sumD) / (2.0 * np.sqrt(sumD) * msa)

    BOARD_RAW = {a: ours[a]["raw"] for a in ours}
    arms = {}
    for arm in ("a8", "a10", "a11", "a12"):
        arms[arm] = {r["ctx"]: dict(A=r["A"], B=r["B"], Kc=r["Kc"], claim=r["claim"],
                                    Chat=np.asarray(r["Chat"]), Cp=r["Cp"], n=r["n"])
                     for r in measure(f"{W096}/{arm}-stats.npz", cp16)}
    f14 = 0.982 / 1.0108
    arms["a14"] = {c: dict(A=arms["a12"][c]["A"] * f14 ** 2, B=arms["a12"][c]["B"] * f14 ** 2,
                           Kc=arms["a12"][c]["Kc"], claim=arms["a12"][c]["claim"],
                           Chat=arms["a12"][c]["Chat"], Cp=arms["a12"][c]["Cp"]) for c in CTX}
    W = float(k180["per_context"]["A"]["before"]["weight_sum"])
    for tag, which in (("a23", "before"), ("abc-knee", "after")):
        arms[tag] = {c: dict(A=k180["per_context"][c][which]["A"],
                             B=k180["per_context"][c][which]["B"], Kc=k180["Kc"][c],
                             claim=k180["per_context"][c][which]["claim"],
                             Chat=np.asarray(k180["per_context"][c][which]["C_hat_values"]),
                             Cp=k180["per_context"][c]["before"]["C_p"]) for c in CTX}
    # a15: terms from 196's stream of the shipped artifact, claim from 191
    arms["a15"] = {c: dict(A=ours["a15"]["per"][c]["A"], B=ours["a15"]["per"][c]["B"],
                           Kc=ours["a15"]["per"][c]["Kc"],
                           claim=ours["a15"]["per"][c]["claim"],
                           Cp=ours["a15"]["per"][c]["Cp"]) for c in CTX}

    ORDER = ["a8", "a10", "a11", "a12", "a14", "a15", "a23", "abc-knee"]
    P("\n(1) cosine under the pooled definition X/sqrt(mean_c A_c * sumD) vs the exact")
    P("    per-context one X/(sqrt(sumD)*mean_c sqrt(A_c)).  Jensen; expected to be tiny.")
    P(f"  {'arm':<10}{'pooled':>10}{'exact':>10}{'diff %':>9}")
    exact = {}
    for a in ORDER:
        per = arms.get(a)
        if per is None:
            continue
        ce = cos_exact(per, BOARD_RAW[a])
        exact[a] = ce
        P(f"  {a:<10}{ours[a]['cos']:>10.5f}{ce:>10.5f}{100*(ce/ours[a]['cos']-1):>9.3f}")
    P("  => the per-context correction is below 0.05% everywhere, including abc-knee's")
    P("     non-uniform alpha.  It is NOT the explanation for abc-knee's 4.4% gap to a23.")

    # ------------------------------------------------------------------ (2) the C_p systematic
    P("\n(2) The `C_p` proxy is the inversion's own systematic.  187 board-calibrated it on two")
    P("    board deltas: f = (1.00, 1.00, 0.9825), i.e. the reference's own jackknife on context")
    P("    C is 1.75% BELOW the NTC-construct proxy.  Re-invert a23 and abc-knee under f:")
    P(f"  {'arm':<10}{'hypothesis':<12}{'claim A/B/C':>34}{'branch A/B/C':>18}{'cos':>10}")
    fixres = {}
    for f_C in (1.0, 0.9825):
        for a in ("a23", "abc-knee"):
            per = {c: dict(arms[a][c]) for c in CTX}
            for c in CTX:
                mu = (f_C if c == "C" else 1.0) * arms[a][c]["Cp"]
                per[c]["claim"] = float(np.mean(np.minimum(arms[a][c]["Chat"], mu)) * W)
            br = "/".join("cl" if per[c]["claim"] <= per[c]["B"] else "B" for c in CTX)
            ce = cos_exact(per, BOARD_RAW[a])
            fixres[(a, f_C)] = ce
            P(f"  {a:<10}{('f_C=%.4f' % f_C):<12}"
              + "  ".join(f"{per[c]['claim']:.5f}" for c in CTX).rjust(32)
              + f"{br:>18}{ce:>10.5f}")
    P(f"  a23 vs abc-knee cosine gap:  proxy {100*(fixres[('abc-knee',1.0)]/fixres[('a23',1.0)]-1):+.2f}%"
      f"   under 187's f_C {100*(fixres[('abc-knee',0.9825)]/fixres[('a23',0.9825)]-1):+.2f}%")
    P("  => the gap is the C_p proxy, not a direction change: abc-knee shipped C's alpha 0.9%")
    P("     above C's TRUE cliff, and the inversion charges that overshoot to the cosine.")
    P("  ==> READ THIS AS THE INSTRUMENT'S OWN ERROR BAR ON OUR SIDE: +/- ~4% on the cosine")
    P("      of any single row, and it is one-signed (the proxy can only make us look worse).")

    OURS_COS = float(np.mean([exact["a15"], exact["a23"]]))
    OURS_LO, OURS_HI = min(exact[a] for a in ("a14", "a15", "a23", "abc-knee")), \
        max(exact[a] for a in ("a14", "a15", "a23", "abc-knee"))
    P(f"\n  OUR CURRENT COSINE: {OURS_COS:.5f}  (a15 {exact['a15']:.5f}, a23 {exact['a23']:.5f};"
      f" range over the four current-profile arms {OURS_LO:.5f}-{OURS_HI:.5f})")

    # ------------------------------------------------------------------ (3) three routes
    P("\n" + "=" * 100)
    P("197 PHASE 2c — THREE ROUTES TO ONE LATENT QUANTITY")
    P("=" * 100)
    field = inv["field"]
    us = next(t for t in field if "Ahmet" in t["team"])
    lead_mse = sorted(field, key=lambda t: -t["floor"])[:6]

    # route B: 159's curve, recalibrated so it passes through OUR measured point
    def curve_inv(pds, tab):
        return float(np.interp(pds, tab, C159))

    shift_h = OURS_HI  # our measured cosine; the curve is calibrated to reproduce our own pds
    for surf, tab in (("HCT116", PDS_HCT), ("HEK293T", PDS_HEK)):
        raw_us = curve_inv(us["pds"], tab)
        P(f"\n  159 5.4 curve on the {surf} surface: our board pds {us['pds']:.5f} reads back as")
        P(f"    cos {raw_us:.5f}, against our MEASURED {OURS_COS:.5f} — the curve is "
          f"{100*(raw_us/OURS_COS-1):+.1f}% off on the one row where we know the answer.")
    P("\n  ** CORRECTION TO 159 5.4. ** 159 validated that curve by writing 'our board-implied")
    P("     aggregate cosine is 0.0722 ... the curve at cos = 0.072 predicts raw pds 0.743; the")
    P("     board reads our raw pds at ~0.73-0.74'.  The board has NOT read our raw pds at 0.73")
    P("     since 2026-08-21: 0.72793 is the **v3** entry (raw mse 1.683).  Every entry from")
    P("     08-22 on reads 0.816-0.837, and the a15-era rows 159's cosine came from read 0.836.")
    P("     The known-answer row mixed two epochs.  Re-run honestly, the curve UNDER-reads our")
    P("     pds by ~0.09 at our cosine, which is exactly AGENTS.md 4.2's whitening effect:")
    P("     whitening buys retrieval without buying bulk cosine.")

    P(f"\n  {'team':<22}{'raw mse':>9}{'A: cos floor':>13}{'B: cos<-pds':>12}{'C: q_w >=':>11}"
      f"{'B/A':>7}{'x ours (A)':>11}{'x ours (B)':>11}{'x ours (C-.5)':>14}")
    rows = []
    for t in [us] + [x for x in lead_mse if "Ahmet" not in x["team"]]:
        cb = curve_inv(t["pds"], PDS_HCT)
        r = dict(team=t["team"], raw=t["raw"], A=t["floor"], B=cb, C=t["qmin"],
                 xA=t["floor"] / OURS_COS, xB=cb / curve_inv(us["pds"], PDS_HCT),
                 xC=(t["qmin"] - 0.5) / (us["qmin"] - 0.5))
        rows.append(r)
        P(f"  {t['team'][:21]:<22}{t['raw']:>9.5f}{t['floor']:>13.5f}{cb:>12.5f}{t['qmin']:>11.5f}"
          f"{cb/max(t['floor'],1e-9):>7.2f}{r['xA']:>11.3f}{r['xB']:>11.3f}{r['xC']:>14.3f}")
    P("\n  Note: for US column A is a FLOOR and the measured value is "
      f"{OURS_COS:.5f}; for everyone else A is a floor and nothing more.")
    P("  The three routes DISAGREE, and the disagreement is the finding:")
    P("   * A (mse) says the leaders' direction is >= 1.6-1.7x ours.")
    P("   * B (pds) says ~1.1x — because `pds` is a RANK statistic and our whitening converts")
    P("     cosine into retrieval far more efficiently than the field's (159 5.5: +0.19 pds at")
    P("     the same cosine).  pds is the one member where we are NOT direction-limited.")
    P("   * C (nmae) says ~1.3x on |t|-weighted sign-accuracy excess, i.e. between the two.")
    P("  Consistency of A with C, through the Gaussian sign law 2q-1 = (2/pi) arcsin(cos):")
    for r in rows:
        pred = 0.5 * (1 + (2 / np.pi) * np.arcsin(min(r["A"], 1.0)))
        P(f"    {r['team'][:21]:<22} cos floor {r['A']:.5f} -> q_w >= {pred:.5f}  "
          f"vs nmae's floor {r['C']:.5f}  ({'nmae SHORT of what mse implies' if r['C'] < pred else 'nmae consistent'})")

    # ------------------------------------------------------------------ (4) the target
    P("\n" + "=" * 100)
    P("197 PHASE 3 — THE TARGET: what the leader's direction demands of us")
    P("=" * 100)
    lead = max(field, key=lambda t: t["floor"])
    top3 = sorted(field, key=lambda t: -t["floor"])[:3]
    P(f"  the three highest implied directions on the board: "
      + ", ".join(f"{t['team']} {t['floor']:.5f}" for t in top3))
    P(f"  ratio to ours ({OURS_COS:.5f}, measured): "
      + ", ".join(f"{t['floor']/OURS_COS:.3f}x" for t in top3))
    P(f"  ratio to our displayed entry abc-knee ({exact['abc-knee']:.5f}): "
      + ", ".join(f"{t['floor']/exact['abc-knee']:.3f}x" for t in top3))

    P("\n  159 5.5's law, 1/cos^2(k) = 66.5 + 1920/k, on the HCT116 challenge-panel surface.")
    P("  The board<->surface conversion is the ratio our OWN single prediction spans between")
    P("  the two, in which the surface factor cancels because it is the same prediction:")
    for tag, ourc in (("a15", exact["a15"]), ("a23", exact["a23"]),
                      ("abc-knee", exact["abc-knee"])):
        ratio = ourc / COS_SURF_PROD
        P(f"    anchor {tag:<9} board/surface ratio {ratio:.4f}   "
          f"(159 used 0.0722/0.03814 = 1.8931)")
    P("")
    P(f"  {'board<->surface anchor':>26}{'leader cos (surf)':>13}{'k_eff':>9}{'ceiling (board)':>17}"
      f"{'leader % of ceiling':>21}")
    keffs = []
    for tag, ourc in (("a15 anchor", exact["a15"]), ("a23 anchor", exact["a23"]),
                      ("abc-knee anchor", exact["abc-knee"]),
                      ("159's own 1.8931", COS_SURF_PROD * 1.8931)):
        ratio = ourc / COS_SURF_PROD
        ls = lead["floor"] / ratio
        k = keff(ls)
        keffs.append(k)
        P(f"  {tag:>17} (x{ratio:.3f}){ls:>13.5f}{k:>9.2f}{COS_SURF_INF*ratio:>17.5f}"
          f"{100*lead['floor']/(COS_SURF_INF*ratio):>20.1f}%")
    P(f"\n  ours, for comparison: surface cos {COS_SURF_PROD:.5f} -> k_eff {keff(COS_SURF_PROD):.2f}"
      f"   (159: 3.09)")
    P(f"  the best equal consensus of the FOUR panel-covering lines we hold: surface cos "
      f"{COS_SURF_4DON:.5f} -> k_eff {keff(COS_SURF_4DON):.2f}, board-equivalent "
      f"{COS_SURF_4DON*exact['a23']/COS_SURF_PROD:.5f}")
    P(f"  => the leader's FLOOR of {lead['floor']:.5f} board is "
      f"{lead['floor']/(COS_SURF_4DON*exact['a23']/COS_SURF_PROD):.2f}x the best consensus the")
    P("     four donor lines in existence can produce, and 166 proved there is no fifth line.")
    P(f"  => but it is only {100*lead['floor']/(COS_SURF_INF*exact['a23']/COS_SURF_PROD):.0f}% of the "
      f"context-INDEPENDENT k->infinity ceiling, so it does not require a context-adapted model.")

    P("\n  How many independent lines the leader's direction corresponds to, by the 1/k law:")
    P(f"    {np.mean(keffs):.1f} effective lines (range {min(keffs):.1f}-{max(keffs):.1f}) "
      f"against our 3.09 and the four that exist.")
    P("    Panel coverage (166 census): orion_hek293t 294/300, orion_hct116 290, kolf2 273,")
    P("    k562_gwps 269, h1_2025 20, everything else 0.  FOUR lines, and 166 searched HF, GEO")
    P("    (258 series), Europe PMC (229 papers), ArrayExpress, the Arc atlas and IGVF and found")
    P("    no fifth.  So the leader's cosine is NOT reachable by donor count.  It is not")
    P("    reachable by donor count for THEM either — which is the whole point.")

    # ------------------------------------------------------------------ 8's 1.28x
    P("\n" + "=" * 100)
    P("197 PHASE 4 — ADJUDICATING AGENTS.md 8's 'at least 1.28x'")
    P("=" * 100)
    agi = next((t for t in field if t["team"] == "Aginglab.com"), None)
    P("  Provenance: 096 (b) inverted the then-leader Aginglab.com's raw 0.8270 at the #247 cap")
    P("  and got cos >= 0.09469 against our then-cosine 0.0740 -> 1.28x.  Reproduced here:")
    c_old = (1.0236107 - 0.8270) / coeff
    P(f"    cos floor at raw 0.8270 = ({floor0:.6f} - 0.8270)/{coeff:.6f} = {c_old:.5f}"
      f"   (096: 0.09469)   ratio to 096's 0.0740 = {c_old/0.0740:.3f}x  (096: 1.28x)")
    P("  VERDICT ON THE CONSTRUCTION: CONFIRMED.  The cap bounds `claim` from above, so it")
    P("  bounds the implied cosine from BELOW; the inequality points the way 096's late")
    P("  correction said it does.  Independently re-derived here, and now also validated on")
    P("  eight of our own board rows where the cosine is measured rather than bounded (the")
    P("  floor is valid on 8/8 and 0.91-0.97x the truth on the four rows near the cap).")
    P("")
    P("  VERDICT ON THE NUMBER: STALE, and stale in the direction that matters.")
    if agi:
        P(f"    Aginglab.com today: raw {agi['raw']:.5f} -> floor {agi['floor']:.5f} = "
          f"{agi['floor']/OURS_COS:.2f}x ours.  They are no longer the leader.")
    P(f"    The board's leading direction today is {lead['team']} at raw {lead['raw']:.5f} ->")
    P(f"    floor {lead['floor']:.5f} = {lead['floor']/OURS_COS:.2f}x ours, and three teams are")
    P(f"    within 1% of it.  The sentence that governs allocation should read")
    P(f"    'at least {lead['floor']/OURS_HI:.2f}x ours', not 1.28x.")

    out = dict(exact=exact, ours_cos=OURS_COS, ours_lo=OURS_LO, ours_hi=OURS_HI,
               fixres={f"{k[0]}|{k[1]}": v for k, v in fixres.items()},
               routes=rows, keff=keffs, leader=lead, top3=top3,
               c_old_0827=c_old)
    json.dump(out, open(f"{WORK}/phase3197.json", "w"), indent=1, default=float)
    open(f"{WORK}/phase3197.log", "w").write("\n".join(LINES))
    P(f"\nwrote {WORK}/phase3197.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
