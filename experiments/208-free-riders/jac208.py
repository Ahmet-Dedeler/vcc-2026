"""208 Phase 1 — ADJUDICATE the two `jac` conventions against the board.

200 could not price its tier-vs-flat-1.00 trade because the two conventions disagree in SIGN
about what a call-set change does to board `jac`, and the only board point available (a21)
tested the *wrong direction* -- a collapse, not an expansion -- while also moving the delivered
profile.  `rc_callset` supplied the missing observation: a LARGE, CLEAN, SINGLE-VARIABLE `jac`
move with a MEASURED `n_pred` change and `pds`/`mse`/`nmae` provably unmoved.

THE MODEL both conventions share (six200.overall's docstring is the source):

    J(R, NP, NR) = R*NR / (NP + NR - R*NR)                                        (M)

`R` = recall (share of the reference's DE genes we call), transferred from h1-test; `NP` = our
call-set size on the challenge; `NR` = the challenge's `n_real`.  The conventions differ ONLY
in `NP`:

    F   fixed-NP    (169's, and 200's incumbent)   NP_arm = NP_base
    S   scaled-NP,  h1 ratio                       NP_arm = NP_base * npred_h1 ratio
    S*  scaled-NP,  DIRECT                         NP_arm = the arm's OWN MEASURED challenge NP

S* is the physical statement: we measured `rc_callset`'s challenge call set at 2,599 and
`rc_hybrid`'s at 4,480 on the built artifacts, against the real control cells, before either
was submitted.  Under F the model is asked to use 4,480 for a set we counted at 2,599.

TWO board points, in OPPOSITE directions:

    a15   -> a21          challenge NP 4,420 -> 1,141   board raw `jac` 0.030080 -> 0.023997
    rc_hy -> rc_callset   challenge NP 4,480 -> 2,599   board raw `jac` 0.03031  -> 0.04043

and a21 ALSO moved the delivered profile (board `nmae` -0.064214 scaled) while the rc pair did
not (`pds`/`mse` bit-identical to sixteen figures).  That difference is the whole story, and
section 3 below separates it by INVERTING the model for each arm's implied BOARD recall and
comparing it to the h1-test recall the convention would have transferred.

    uv run python jac208.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/153-de-bench"))
from debench import ANCHOR, Surface  # noqa: E402

E202 = ROOT / "experiments/202-callset-rider/work"
E200 = ROOT / "experiments/200-production-ablation"
SPAN_J = ANCHOR["jac"][1] - ANCHOR["jac"][0]
SLOPE_J = 2.7263                       # docs/BOARD-LEDGER.md's fitted scaled = slope*(raw-b)

# --- the board, as measured -------------------------------------------------------------
# rc pair: BOARD-RESULT.md (185/202); challenge n_pred from 202 falsifier 2 (probe191, 40
# targets x 3 contexts, against the real control cells).
# a15/a21 pair: docs/LOG.md 2026-09-01 (scaled members) and docs/CURRENT_STATE.md:755
# ("a21: n_pred 4,420 -> 1,141, fid 0.5034 -> 0.4945"), both from the same probe.
BOARD = {
    "a15": dict(jac_raw=0.030080, n_pred=4420, profile_moved=False),
    "a21": dict(jac_raw=None, jac_scaled_d=-0.016583, n_pred=1141, profile_moved=True),
    "rc_hybrid": dict(jac_raw=0.03031, n_pred=4480, profile_moved=False),
    "rc_callset": dict(jac_raw=0.04043, n_pred=2599, profile_moved=False),
}
BOARD["a21"]["jac_raw"] = BOARD["a15"]["jac_raw"] + BOARD["a21"]["jac_scaled_d"] / SLOPE_J
# h1-test recall.  The rc pair is measured below off 202's own tables; a15 -> a21 is quoted
# from AGENTS.md section 11.1 rule 13's amendment / experiment 200 (a21 has no h1-test table
# on disk, it is a `resample_mix` arm).
H1_RECALL_A21 = dict(a15=0.669, a21=0.271)
H1_NPRED_A21 = dict(a15=6034.0, a21=1913.0)


def model(R, NP, NR):
    hit = R * NR
    return hit / (NP + NR - hit)


def invert_recall(J, NP, NR):
    """The board recall the model needs to produce `J` at the arm's own measured `NP`."""
    return J * (NP + NR) / (NR * (1.0 + J))


def solve_nr(R, NP, target, lo=1.0, hi=50000.0):
    f = lambda nr: model(R, NP, nr) - target  # noqa: E731
    if f(lo) * f(hi) > 0:
        return float("nan")
    return float(brentq(f, lo, hi, xtol=1e-9))


def per_pert(s):
    nr, npd, rec = [], [], []
    for p in range(len(s.perts)):
        a, b = s.real_sig[p], s.pred_sig[p] & s.adj[p]
        n_real = int(a.sum())
        nr.append(n_real)
        npd.append(int(b.sum()))
        rec.append(int((a & b).sum()) / n_real if n_real else np.nan)
    return np.asarray(nr, float), np.asarray(npd, float), np.asarray(rec, float)


def main() -> None:
    out = {}
    S = {a: Surface(E202 / f"de202.{a}.npz") for a in ("rc_hybrid", "rc_callset")}
    stats = {}
    print("=== 1. h1-test, 202's own emitted tables (the objects the board scored) ===\n")
    print(f"{'arm':<12}{'n_real med':>11}{'n_real mean':>12}{'n_pred med':>11}{'recall':>9}"
          f"{'jac h1':>9}{'board NP':>10}{'board jac':>11}")
    for a, s in S.items():
        nr, npd, rec = per_pert(s)
        stats[a] = dict(n_real=nr, n_pred=npd, recall=rec,
                        recall_mean=float(np.nanmean(rec)),
                        n_pred_med=float(np.median(npd)),
                        n_real_med=float(np.median(nr)),
                        n_real_mean=float(nr.mean()), jac_h1=float(s.jac()))
        print(f"{a:<12}{np.median(nr):>11.0f}{nr.mean():>12.0f}{np.median(npd):>11.0f}"
              f"{np.nanmean(rec):>9.5f}{s.jac():>9.5f}{BOARD[a]['n_pred']:>10,}"
              f"{BOARD[a]['jac_raw']:>11.5f}")
    Rb, Ra = stats["rc_hybrid"]["recall_mean"], stats["rc_callset"]["recall_mean"]
    NPb, NPa = BOARD["rc_hybrid"]["n_pred"], BOARD["rc_callset"]["n_pred"]
    D_RC = BOARD["rc_callset"]["jac_raw"] - BOARD["rc_hybrid"]["jac_raw"]
    D_A21 = BOARD["a21"]["jac_raw"] - BOARD["a15"]["jac_raw"]
    ratio_h1 = stats["rc_callset"]["n_pred_med"] / stats["rc_hybrid"]["n_pred_med"]
    print(f"\nh1-test n_pred ratio {ratio_h1:.4f}; on the board that would be "
          f"{NPb*ratio_h1:,.0f} against the MEASURED {NPa:,} "
          f"({100*(NPb*ratio_h1/NPa-1):+.1f}%).  The `n_pred` transfer is itself imperfect.")
    print(f"realised board d raw `jac`:  rc pair {D_RC:+.6f}   a15->a21 {D_A21:+.6f}")

    # ---- 2. the three conventions, on BOTH board pairs --------------------------------
    print("\n=== 2. the conventions, retrodicting BOTH board pairs ===")
    rows = []
    for pair, base, arm, NP0, NP1, Rb_, Ra_, ratio, real_d in (
            ("rc_hybrid->rc_callset", "rc_hybrid", "rc_callset", NPb, NPa, Rb, Ra,
             ratio_h1, D_RC),
            ("a15->a21", "a15", "a21", BOARD["a15"]["n_pred"], BOARD["a21"]["n_pred"],
             H1_RECALL_A21["a15"], H1_RECALL_A21["a21"],
             H1_NPRED_A21["a21"] / H1_NPRED_A21["a15"], D_A21)):
        NR = solve_nr(Rb_, NP0, BOARD[base]["jac_raw"])
        print(f"\n{pair}   (NR calibrated on the base point: NR* = {NR:.1f}; "
              f"profile moved: {BOARD[arm]['profile_moved']})")
        print(f"{'convention':<24}{'NP used':>10}{'pred jac':>10}{'pred d':>11}"
              f"{'realised d':>12}{'realised/pred':>15}{'sign':>7}")
        for cname, NPuse in (("F  fixed-NP", NP0), ("S  scaled-NP (h1)", NP0 * ratio),
                             ("S* scaled-NP (direct)", NP1)):
            pred = model(Ra_, NPuse, NR)
            d = pred - BOARD[base]["jac_raw"]
            r = real_d / d if d else float("nan")
            sg = "OK" if np.sign(d) == np.sign(real_d) else "WRONG"
            rows.append(dict(pair=pair, convention=cname, NR=NR, NP=NPuse, pred=pred,
                             d=d, realised=real_d, ratio=r, sign=sg))
            print(f"{cname:<24}{NPuse:>10,.0f}{pred:>10.5f}{d:>+11.5f}{real_d:>+12.5f}"
                  f"{r:>15.2f}{sg:>7}")
        out.setdefault("NR_star", {})[pair] = NR

    # ---- 3. the exponent that reconciles, and the reconciliation itself ----------------
    print("\n=== 3. the NP exponent theta:  NP_arm = NP_base * (NP1/NP0)^theta ===")
    print("     theta = 0 is F, theta = 1 is S*.  Solved so the model hits the realised d.")
    for pair, base, arm, NP0, NP1, Rb_, Ra_, real_d in (
            ("rc_hybrid->rc_callset", "rc_hybrid", "rc_callset", NPb, NPa, Rb, Ra, D_RC),
            ("a15->a21", "a15", "a21", BOARD["a15"]["n_pred"], BOARD["a21"]["n_pred"],
             H1_RECALL_A21["a15"], H1_RECALL_A21["a21"], D_A21)):
        NR = solve_nr(Rb_, NP0, BOARD[base]["jac_raw"])
        tgt = BOARD[base]["jac_raw"] + real_d
        g = lambda th: model(Ra_, NP0 * (NP1 / NP0) ** th, NR) - tgt  # noqa: E731
        try:
            th = float(brentq(g, -3.0, 6.0, xtol=1e-8))
        except ValueError:
            th = float("nan")
        print(f"  {pair:<24} theta = {th:+.3f}")
        out.setdefault("theta", {})[pair] = th

    print("\n=== 4. RECONCILIATION — invert the model for each arm's IMPLIED BOARD recall ===")
    print("     at the arm's OWN measured challenge NP, and compare with h1-test's recall.")
    print(f"{'arm':<12}{'board NP':>10}{'board jac':>11}{'implied R':>11}{'h1 R':>9}"
          f"{'implied/h1':>12}   profile")
    NR_ALL = out["NR_star"]["rc_hybrid->rc_callset"]
    h1r = dict(rc_hybrid=Rb, rc_callset=Ra, a15=H1_RECALL_A21["a15"],
               a21=H1_RECALL_A21["a21"])
    imp = {}
    for a in ("a15", "rc_hybrid", "rc_callset", "a21"):
        R = invert_recall(BOARD[a]["jac_raw"], BOARD[a]["n_pred"], NR_ALL)
        imp[a] = R
        print(f"{a:<12}{BOARD[a]['n_pred']:>10,}{BOARD[a]['jac_raw']:>11.5f}{R:>11.5f}"
              f"{h1r[a]:>9.5f}{R/h1r[a]:>12.3f}   "
              f"{'MOVED' if BOARD[a]['profile_moved'] else 'pinned'}")
    print(f"\n  NR held at {NR_ALL:.1f} for all four (calibrated on rc_hybrid).")
    out["implied_recall"] = imp
    out["h1_recall"] = h1r

    # ---- 5. what this does to 200's trade ---------------------------------------------
    print("\n=== 5. re-pricing 200's tier-vs-flat-1.00 trade under the surviving convention ===")
    tier_vs_flat = dict(reach=-0.000344, fid=+0.001387, jac_fixed=-0.006590,
                        jac_scaled=+0.001576, nmae=-0.000034)
    s_fixed = sum(tier_vs_flat[k] for k in ("reach", "fid", "jac_fixed", "nmae"))
    s_scaled = sum(tier_vs_flat[k] for k in ("reach", "fid", "jac_scaled", "nmae"))
    print("  200's sign convention: POSITIVE = what REMOVING the tier would COST us.")
    print(f"  under F  (fixed-NP) : SUM {s_fixed:+.6f}  -> dropping the tier GAINS "
          f"{-s_fixed:+.6f}")
    print(f"  under S* (scaled-NP): SUM {s_scaled:+.6f}  -> dropping the tier GAINS "
          f"{-s_scaled:+.6f}")
    out["tier_vs_flat"] = dict(sum_fixed=s_fixed, sum_scaled=s_scaled)

    out.update(board=BOARD, rows=rows, ratio_h1=ratio_h1, span_jac=SPAN_J,
               d_rc=D_RC, d_a21=D_A21,
               h1={a: {k: v for k, v in stats[a].items()
                       if not isinstance(v, np.ndarray)} for a in stats})
    json.dump(out, open(WORK / "jac208.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
