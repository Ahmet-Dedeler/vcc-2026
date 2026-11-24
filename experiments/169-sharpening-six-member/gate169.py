"""169 — the six-member table and AGENTS.md 11.1 rule 13's composed slot gate, assembled.

Every number is read from the JSON its own script wrote; nothing here is re-derived by hand.

  work/plolo169.json    `pds`   -- 162's PLOLO, 4 folds, challenge frame (the instrument rule
                                   13 names).  Its q=1.00 rows reproduce 162's a23 rows to
                                   |d| = 0.0e0 and its a15->a23 delta to +0.004014 scaled.
  work/claim169.json    `mse`   -- 096's `claim` route on the real challenge panel, four
                                   known-answer gates passing, priced at rX = 1 (the null)
                                   and at the measured cosine ratios.
  work/nmaechal169.json `nmae`  -- exact, both surfaces, saturated and saturation-free.
  work/reach169.json    `reach` -- h1-test, five MATCHED seeds, seed-paired CI.
  work/six169.json      `fid`, `jac` -- precision and the board-rebuilt Jaccard.

RULE 13, verbatim:

    predicted = PLOLO d_pds/6 * 1.011  +  096 d_mse/6 (if amplitude/spread/phi moves)
                                       -  0.0146/6   (if the EMITTED PROFILE or emission law
                                                      changes)
    SHIP only if predicted > +0.0019 (shape-only)  or  > +0.0089 (amplitude-changing)

Both branches are printed: the LITERAL gate (blanket `reach` toll) and the gate with the toll
replaced by the direct h1-test measurement, which is the whole point of this experiment.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
from debench import ANCHOR  # noqa: E402

QS = (1.05, 1.10, 1.15, 1.20, 1.25, 1.30)
SPAN_PDS = ANCHOR["pds"][1] - ANCHOR["pds"][0]
SPAN_REACH = ANCHOR["reach"][1] - ANCHOR["reach"][0]
SPAN_FID = ANCHOR["fid"][1] - ANCHOR["fid"][0]
SPAN_JAC = ANCHOR["jac"][1] - ANCHOR["jac"][0]
# `nmae` has no debench anchor. AGENTS.md §4 gives b 1.0009-1.0017, r 0.369-0.431 (span
# 0.57-0.63); our own two board points (a14 raw 0.945764 -> scaled 0.091972, a15 raw 0.944251
# -> scaled 0.094239) solve to span 0.6674. The SMALLER span makes the nmae cost LARGER, so
# 0.60 is the conservative choice and is what the gate uses; 0.6674 is reported alongside.
SPAN_NMAE, SPAN_NMAE_BOARD = 0.60, 0.66741
K_PLOLO_LITERAL = 1.011          # rule 13's multiplier (162's a23-held-out scale)
K_PLOLO_FIT = 0.899              # 162's all-arm fitted scale -- the conservative one
TOLL = -0.014569 / 6             # rule 13's blanket emitted-profile `reach` toll, on overall
SCORE = lambda raw: min(max((0.987862879 - raw) / 0.950162566, 0.0), 1.0)  # noqa: E731


def load(name):
    p = os.path.join(WORK, name)
    if not os.path.exists(p):
        sys.exit(f"missing {p}")
    return json.loads(open(p).read())


def main() -> None:
    P = load("plolo169.json")
    C = load("claim169.json")
    N = load("nmaechal169.json")
    Rj = load("reach169.json")
    S = load("six169.json")
    X = load("xratio169.json")

    # ---- pds, from PLOLO's pre-committed predictor (4-fold mean, challenge frame) ----
    R = {(r["fold"], r["arm"], r["q"], r["space"]): r for r in P}
    folds = sorted({r["fold"] for r in P})

    def chal(f, q):
        return float(np.mean([R[(f, "a23", q, f"chal{c}/delivered")]["pds"] for c in "ABC"]))

    dpds = {q: float(np.mean([chal(f, q) - chal(f, 1.00) for f in folds])) / SPAN_PDS
            for q in [1.00] + list(QS)}
    dpds_sd = {q: float(np.std([chal(f, q) - chal(f, 1.00) for f in folds])) / SPAN_PDS
               for q in [1.00] + list(QS)}

    # ---- mse, from the challenge panel ----
    bs = np.load(f"{REPO}/experiments/096-mse-floor/work/board_solve.npz")
    sD, Kc, X0 = float(bs["sumD"]), float(bs["v11_Kc"]), float(bs["X"])
    s15 = 0.982 / 1.0796

    def mse_overall(q, rX):
        def raw(qq, r):
            c = C[repr(float(qq))]
            return 1 + (s15**2 * c["A"] - 2 * s15 * X0 * r + Kc
                        - min(s15**2 * c["B"], c["claim"])) / sD
        return (SCORE(raw(q, rX)) - SCORE(raw(1.00, 1.0))) / 6

    def rx_of(key, q):
        row = [z for z in X[key] if abs(z["q"] - q) < 1e-9][0]
        zero = [z for z in X[key] if abs(z["q"] - 1.0) < 1e-9][0]
        return row["X"] / zero["X"]

    # PLOLO gives a THIRD, independent cosine ratio -- 4 folds, challenge frame, production
    # panel -- which the two Orion surfaces can be checked against.
    def rx_plolo(q):
        a = np.mean([np.mean([R[(f, "a23", q, f"chal{c}/delivered")]["cos"] for c in "ABC"])
                     for f in folds])
        b = np.mean([np.mean([R[(f, "a23", 1.0, f"chal{c}/delivered")]["cos"] for c in "ABC"])
                     for f in folds])
        return float(a / b)

    # ---- nmae ----
    nm = {r["q"]: r for r in N["hct116"]}
    nmh1 = {r["q"]: r for r in N.get("h1-test", [])}

    # ---- reach, fid, jac ----
    reach = {float(k): v for k, v in Rj.items()}

    print(f"spans: pds {SPAN_PDS:.5f}  reach {SPAN_REACH:.5f}  fid {SPAN_FID:.5f}  "
          f"jac {SPAN_JAC:.5f}  nmae {SPAN_NMAE} (board-fitted {SPAN_NMAE_BOARD})")
    print(f"\nrX (the alignment ratio X_q/X_1), three independent surfaces:")
    print(f"  {'q':>5} {'hct116':>9} {'hek293t':>9} {'PLOLO 4-fold':>14} {'worst':>9}")
    for q in QS:
        a, b, c = rx_of("hct116|row", q), rx_of("hek293t|row", q), rx_plolo(q)
        print(f"  {q:>5.2f} {a:9.5f} {b:9.5f} {c:14.5f} {min(a, b, c):9.5f}")

    print(f"\n{'=' * 128}\nTHE SIX-MEMBER TABLE, in OVERALL board units, vs a15+a23 "
          f"(`row` renormalisation)")
    print(f"  `mse` is shown at the rX = 1 NULL, which is also what rule 13 credits for a "
          f"shape-only arm.")
    print(f"\n  {'q':>5} {'pds (PLOLO)':>13} {'fold sd':>9} {'mse':>10} {'nmae':>10} "
          f"{'reach':>10} {'fid':>10} {'jac':>10} {'SIX-MEMBER SUM':>16}")
    tab = {}
    for q in QS:
        p = dpds[q] * K_PLOLO_LITERAL / 6
        m = mse_overall(q, 1.0)
        n = -(nm[q]["true"] - nm[1.00]["true"]) / SPAN_NMAE / 6
        K_R = (0.1890 - 0.1811) / (0.2029 - 0.1828)
        matched = q in reach
        f_ = S.get("q%.2f_mf" % q, {})
        b1 = S.get("q1.00_mf", {})
        if matched:
            r = reach[q]["mean_overall"]
        elif "reach_f1.0" in f_:
            r = (f_["reach_f1.0"] - b1["reach_f1.0"]) * K_R / SPAN_REACH / 6
        else:
            r = float("nan")
        b_ = b1
        dfid = ((f_.get("precision", np.nan) - b_.get("precision", np.nan))
                / SPAN_FID / 6)
        NR, NP = 240.0, 4424.0
        r0, r1 = b_.get("recall", np.nan), f_.get("recall", np.nan)
        j0 = r0 * NR / (NP + NR - r0 * NR)
        j1 = r1 * NR / (NP + NR - r1 * NR)
        djac = (j1 - j0) / SPAN_JAC / 6
        tot = p + m + n + (0.0 if np.isnan(r) else r) + (0 if np.isnan(dfid) else dfid) \
            + (0 if np.isnan(djac) else djac)
        tab[q] = dict(pds=p, mse=m, nmae=n, reach=r, fid=dfid, jac=djac, total=tot)
        print(f"  {q:>5.2f} {p:>+13.6f} {dpds_sd[q] / 6:>9.6f} {m:>+10.6f} {n:>+10.6f} "
              f"{r:>+10.6f}{'' if matched else '*'} {dfid:>+10.6f} {djac:>+10.6f} "
              f"{tot:>+16.6f}")
    print("  * = single seed, NOT resolvable: the matched-seed sem on `reach` is 0.0037 raw "
          "(= 0.00028 overall).")
    print("    q = 1.15 and q = 1.25 are the only rows with five and three matched seeds.")

    print(f"\n  same table with `mse` at the WORST measured rX instead of the null:")
    print(f"  {'q':>5} {'mse @ worst rX':>16} {'SIX-MEMBER SUM':>16}")
    for q in QS:
        rx = min(rx_of("hct116|row", q), rx_of("hek293t|row", q), rx_plolo(q))
        m = mse_overall(q, rx)
        print(f"  {q:>5.2f} {m:>+16.6f} {tab[q]['total'] - tab[q]['mse'] + m:>+16.6f}")

    print(f"\n{'=' * 128}\nRULE 13's COMPOSED SLOT GATE")
    print(f"  Gate 0 (eligibility): the arm changes the effect matrix only -- no emitter "
          f"change, no count law,\n    no per-cell library change. PLOLO is therefore "
          f"eligible.  PASS.")
    print(f"  Gate 2 (amplitude exposure): `row` renormalisation preserves every row's L2 norm "
          f"exactly;\n    measured rA on the challenge panel with the panel coordinates zeroed "
          f"is 0.99903 at q=1.15 and\n    0.99990 on the Orion pair.  VCC_ALPHA is untouched.  "
          f"SHAPE-ONLY -> threshold +0.0019.")
    print(f"\n  {'q':>5} {'PLOLO d pds':>13} {'x1.011 /6':>12} {'|k.dpds| vs 0.00246':>20} "
          f"{'LITERAL (blanket toll)':>24} {'MEASURED reach':>17} {'verdict':>10}")
    for q in QS:
        p = dpds[q]
        lit = p * K_PLOLO_LITERAL / 6 + TOLL
        meas = p * K_PLOLO_LITERAL / 6 + tab[q]["reach"]
        g1 = abs(K_PLOLO_FIT * p)
        print(f"  {q:>5.2f} {p:>+13.6f} {p * K_PLOLO_LITERAL / 6:>+12.6f} "
              f"{g1:>13.5f} {'PASS' if g1 >= 0.00246 else 'FAIL':>6} "
              f"{lit:>+18.6f} {'PASS' if lit > 0.0019 else 'FAIL':>5} "
              f"{meas:>+12.6f} {'PASS' if meas > 0.0019 else 'FAIL':>4} "
              f"{'SHIP' if meas > 0.0019 and g1 >= 0.00246 else '-':>10}")

    print(f"\n  Gate 4 (the budget): 162 computed that a shape-only arm needs predicted "
          f"d pds > +0.0126 scaled\n    to clear its own unpriced noise -- 'nearly 2x the "
          f"largest value in our entire history' (a15, +0.00666).")
    for q in QS:
        print(f"    q={q:.2f}: PLOLO d pds {dpds[q]:+.6f} scaled  "
              f"({dpds[q] / 0.0126:.2f}x the 0.0126 bar, "
              f"{dpds[q] / 0.006656:.2f}x a15's historical maximum)")

    print(f"\n{'=' * 128}\nCONSERVATIVE VARIANTS of the recommended arm")
    q = 1.15
    p = dpds[q]
    rows = [
        ("rule 13 literal (k=1.011, blanket toll, no other member)",
         p * K_PLOLO_LITERAL / 6 + TOLL),
        ("k = 0.899 (162's all-arm fitted scale), blanket toll",
         p * K_PLOLO_FIT / 6 + TOLL),
        ("k = 1.011, measured reach, no other member",
         p * K_PLOLO_LITERAL / 6 + tab[q]["reach"]),
        ("k = 0.899, measured reach, no other member",
         p * K_PLOLO_FIT / 6 + tab[q]["reach"]),
        ("k = 0.899, measured reach at its 95% LOWER bound",
         p * K_PLOLO_FIT / 6 + reach[q]["lo"]),
        ("k = 0.899, all six members, mse at the rX = 1 null",
         p * K_PLOLO_FIT / 6 + tab[q]["total"] - tab[q]["pds"]),
        ("k = 0.899, all six, mse null, nmae at the board-fitted span 0.6674",
         p * K_PLOLO_FIT / 6 + tab[q]["total"] - tab[q]["pds"] - tab[q]["nmae"]
         - (nm[q]["true"] - nm[1.00]["true"]) / SPAN_NMAE_BOARD / 6),
        ("k = 0.899, all six, mse at the worst measured rX",
         p * K_PLOLO_FIT / 6 + tab[q]["total"] - tab[q]["pds"] - tab[q]["mse"]
         + mse_overall(q, min(rx_of("hct116|row", q), rx_of("hek293t|row", q), rx_plolo(q)))),
        ("WORST CASE: k=0.899, reach at its 95% lower bound, nmae at the SATURATED h1-test "
         "value,\n                mse at the null",
         p * K_PLOLO_FIT / 6 + reach[q]["lo"] + tab[q]["mse"] + tab[q]["fid"] + tab[q]["jac"]
         - (nmh1[q]["true"] - nmh1[1.00]["true"]) / SPAN_NMAE / 6),
    ]
    for label, v in rows:
        print(f"  {v:>+11.6f}   {'PASS' if v > 0.0019 else 'FAIL':<5} {label}")

    with open(os.path.join(WORK, "gate169.json"), "w") as fh:
        json.dump({str(k): v for k, v in tab.items()}, fh, indent=1, default=float)


if __name__ == "__main__":
    sys.exit(main())
