"""162 — the deliverable: predicted vs realised board deltas, per member, with LOO error bars.

WHAT IS BEING CLAIMED, AND HOW IT COULD FAIL
--------------------------------------------
The claim is that PLOLO (`plolo.py`) predicts the sign and size of a board delta for an arm it
has never seen.  Every element below exists so that claim can lose:

  * the predictor is the UNWEIGHTED mean over all four folds, fixed in advance.  No fold is
    selected, and the per-fold scatter is reported as the instrument's own error bar.  On a23
    that scatter is larger than the effect, and the table says so.
  * the frame is the CHALLENGE frame (mean of contexts A/B/C), pre-committed by the previous
    runtime on the mechanistic ground that it is the frame production renders into and that its
    render residual (2.88/2.90/2.87 %) reproduces `A23-PREDICTION.md`'s independently measured
    2.83-3.35 %.  `library` and `intended` are reported as a labelled post-hoc sensitivity.
  * the fit is ONE parameter (a scale `k`, through the origin) and it is scored by
    LEAVE-ONE-ARM-OUT, never in sample.
  * the null model is `predict zero`.  A member only counts as retrodicted if its LOO RMS
    residual beats the RMS of the realised deltas themselves.  Three of the four measurable
    members FAIL that test, and the table reports them as failures rather than dropping them.
  * an EXACT permutation control (all 720 relabelings of 6 arms) asks what the LOO error would
    be if the prediction carried no arm-specific information.  This is not invariant by
    construction: shuffling destroys the arm-to-prediction pairing, which is the entire claim.

WHAT IS EXCLUDED, AND WHY
-------------------------
  * `a13`   — donor-agreement shrinkage is not reconstructible in this pipeline.
  * `a2-prune` — realised MEMBER deltas do not exist (v11's decomposition was never recorded).
    It is carried as an overall-only row and excluded from every per-member fit.
  * `a11`'s `mse` — its base `a10` sits at scaled `mse` exactly 0.0, which is the CLAMP floor
    (`mse` is clamped to [0,1]).  The realised `+0.082041` is therefore a censored lower bound
    on the true move, not a measurement.  Excluded from the `mse` fit, shown in the table.
  * `a8`/`a10`/`a21`/`a22` — emission arms.  PLOLO has no emitter, so it predicts exactly 0.0
    for all four.  They are NOT folded into any RMS: predicting zero and reading zero is a test
    of the channel classification, not of the instrument, and averaging them in would deflate
    the error bar with free points.  They are reported in their own block, where they carry the
    single most important negative result in this experiment.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np

import board

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

# Board reference spans |r - b|, `docs/reference/vcc2026-metrics-BRIEF.md` §8.
# A fitted scale `k` absorbs any error in these, so they only set the units of the
# UNCALIBRATED column.  `mse`'s value is not the table midpoint but the exact ratio implied by
# `A23-BOARD-RESULT.md`, which records both the raw and the scaled move (0.001457 / 0.001524).
SPAN = dict(pds=0.452, mse=0.956, nmae=0.600, reach=0.895, fid=0.305, jac=0.365)

FOLDS = ("orion_hct116", "orion_hek293t", "k562_gwps", "kolf2")
GEOMETRY = ("a11", "a12", "a14", "a15", "a19", "a23")     # reconstructible, member-resolved
EMISSION = ("a8", "a10", "a21", "a22")                    # predicted exactly 0 by construction

# member -> (endpoint field, sign).  `sign` maps a rise in the endpoint to a rise in the SCALED
# board member: `pds` up is good; `nmse` and `nmae_k` are errors, so down is good; `reach_k` up
# is good.  Board scaling is then `delta_endpoint * sign / SPAN[member]`.
CHANNEL = dict(pds=("pds", +1.0), mse=("nmse_row", -1.0),
               nmae=("nmae_k240", -1.0), reach=("reach_k240", +1.0))


# ------------------------------------------------------------------ predictor
def load(space_mode: str = "challenge"):
    rows = json.loads((WORK / "arms2.json").read_text())
    d = {(r["fold"], r["arm"], r["space"]): r for r in rows}

    def value(fold, arm, field):
        if space_mode == "challenge":
            return float(np.mean([d[(fold, arm, f"chal{c}/delivered")][field]
                                  for c in "ABC"]))
        return float(d[(fold, arm, f"{space_mode}")][field])
    return d, value


def per_fold_prediction(value, d, arm, member, space_mode, folds=FOLDS):
    """Board-scaled predicted delta of `arm` vs its base, one number per usable fold."""
    field, sign = CHANNEL[member]
    base = board.ARMS[arm]["base"]
    out = {}
    for f in folds:
        probe = (f, arm, "intended")
        if probe not in d or (f, base, "intended") not in d:
            continue
        out[f] = sign * (value(f, arm, field) - value(f, base, field)) / SPAN[member]
    return out


# ------------------------------------------------------------------- fitting
def loo_scale(x, y):
    """Leave-one-out predictions under a one-parameter through-origin fit y = k x."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    yhat = np.empty_like(y)
    ks = np.empty_like(y)
    for i in range(len(x)):
        m = np.ones(len(x), bool)
        m[i] = False
        denom = float(x[m] @ x[m])
        ks[i] = float(x[m] @ y[m]) / denom if denom > 0 else 0.0
        yhat[i] = ks[i] * x[i]
    k_in = float(x @ y) / float(x @ x) if float(x @ x) > 0 else 0.0
    return k_in, ks, yhat


def rms(v):
    return float(np.sqrt(np.mean(np.asarray(v, float) ** 2)))


def permutation_p(x, y):
    """Exact p: over all relabelings of the arms, how often is LOO RMS at least as good?"""
    base = rms(y - loo_scale(x, y)[2])
    hits = tot = 0
    for perm in itertools.permutations(range(len(x))):
        xp = np.asarray(x, float)[list(perm)]
        tot += 1
        if rms(y - loo_scale(xp, y)[2]) <= base:
            hits += 1
    return base, hits / tot, tot


# -------------------------------------------------------------------- report
def member_table(member, space_mode="challenge", exclude=(), folds=FOLDS):
    d, value = load(space_mode)
    arms, x, y, sd = [], [], [], []
    for arm in GEOMETRY:
        if arm in exclude:
            continue
        realised = board.delta(arm).get(member)
        if realised is None:
            continue
        pf = per_fold_prediction(value, d, arm, member, space_mode, folds)
        if len(pf) < len(folds):
            continue
        arms.append(arm)
        x.append(float(np.mean(list(pf.values()))))
        sd.append(float(np.std(list(pf.values()), ddof=1)) if len(pf) > 1 else 0.0)
        y.append(realised)
    x, y, sd = np.array(x), np.array(y), np.array(sd)
    k_in, ks, yhat = loo_scale(x, y)
    resid = y - yhat
    null = rms(y)
    pr, pval, ntot = permutation_p(x, y)
    return dict(member=member, arms=arms, x=x, y=y, sd=sd, k_in=k_in, ks=ks,
                yhat=yhat, resid=resid, loo_rms=rms(resid), null_rms=null,
                sign_hits=int(np.sum(np.sign(x) == np.sign(y))),
                perm_p=pval, perm_n=ntot,
                raw_sign_hits=int(np.sum(np.sign(x) == np.sign(y))))


def print_member(t, exclude_note=""):
    m = t["member"]
    print(f"\n### {m}   (span {SPAN[m]}, {len(t['arms'])} arms{exclude_note})")
    print(f"{'arm':<8} {'predicted':>11} {'fold sd':>10} {'realised':>11} "
          f"{'LOO pred':>10} {'LOO resid':>11} {'k_-a':>7}")
    for i, a in enumerate(t["arms"]):
        print(f"{a:<8} {t['x'][i]:+11.6f} {t['sd'][i]:10.6f} {t['y'][i]:+11.6f} "
              f"{t['yhat'][i]:+10.6f} {t['resid'][i]:+11.6f} {t['ks'][i]:7.3f}")
    verdict = "RETRODICTS" if t["loo_rms"] < t["null_rms"] else "FAILS (no better than zero)"
    print(f"  in-sample scale k = {t['k_in']:.3f}   sign agreement {t['sign_hits']}/"
          f"{len(t['arms'])}")
    print(f"  LOO RMS residual  = {t['loo_rms']:.6f}   vs predict-zero null "
          f"{t['null_rms']:.6f}   -> {verdict}")
    print(f"  exact permutation p = {t['perm_p']:.4f}  ({t['perm_n']} relabelings)")


def main():
    problems = board.check()
    if problems:
        raise SystemExit(f"board transcription check FAILED: {problems}")
    print("board transcription check: PASS (mean of 6 member deltas == overall, 12/12 arms)")

    print("\n" + "=" * 78)
    print("PART 1 — per-member retrodiction, challenge frame, 4-fold mean, LOO-scored")
    print("=" * 78)
    tables = {}
    for m in ("pds", "mse", "nmae", "reach"):
        exclude = ("a11",) if m == "mse" else ()
        note = ", a11 excluded: base at mse clamp floor" if m == "mse" else ""
        tables[m] = member_table(m, exclude=exclude)
        print_member(tables[m], note)

    print("\n" + "=" * 78)
    print("PART 2 — emission arms: PLOLO predicts exactly 0.0 for all six members")
    print("=" * 78)
    print(f"{'arm':<8} " + " ".join(f"{m:>10}" for m in board.MEMBERS) + f"{'overall':>11}")
    for a in EMISSION:
        dd = board.delta(a)
        print(f"{a:<8} " + " ".join(f"{dd[m]:+10.6f}" for m in board.MEMBERS)
              + f"{dd['overall']:+11.7f}")
    em = {m: rms([board.delta(a)[m] for a in EMISSION]) for m in board.MEMBERS}
    print("\nRMS of the realised move PLOLO calls zero:")
    print("  " + "  ".join(f"{m} {em[m]:.6f}" for m in board.MEMBERS))

    print("\n" + "=" * 78)
    print("PART 3 — frame sensitivity (post-hoc; the pre-committed frame is `challenge`)")
    print("=" * 78)
    print(f"{'member':<8} " + " ".join(f"{s:>22}" for s in
                                       ("challenge", "library/delivered", "intended")))
    for m in ("pds", "mse", "nmae", "reach"):
        cells = []
        for s in ("challenge", "library/delivered", "intended"):
            try:
                t = member_table(m, space_mode=s, exclude=("a11",) if m == "mse" else ())
                cells.append(f"LOO {t['loo_rms']:.5f} sgn {t['sign_hits']}/{len(t['arms'])}")
            except Exception as exc:                        # noqa: BLE001
                cells.append(f"n/a ({exc})")
        print(f"{m:<8} " + " ".join(f"{c:>22}" for c in cells))

    print("\n" + "=" * 78)
    print("PART 4 — overall, and the gap the unpriced members leave")
    print("=" * 78)
    d, value = load("challenge")
    print(f"{'arm':<8} {'pred pds+mse':>13} {'real pds+mse':>13} {'real other 4':>13} "
          f"{'real overall':>13} {'pred overall':>13}")
    rows = []
    for a in GEOMETRY:
        dd = board.delta(a)
        pred = 0.0
        for m in ("pds", "mse"):
            pf = per_fold_prediction(value, d, a, m, "challenge")
            pred += tables[m]["k_in"] * float(np.mean(list(pf.values())))
        real_pm = dd["pds"] + dd["mse"]
        real_other = sum(dd[m] for m in ("nmae", "fid", "reach", "jac"))
        rows.append((a, pred / 6, real_pm / 6, real_other / 6, dd["overall"]))
        print(f"{a:<8} {pred/6:+13.6f} {real_pm/6:+13.6f} {real_other/6:+13.6f} "
              f"{dd['overall']:+13.6f} {pred/6:+13.6f}")
    ov = np.array([r[4] for r in rows])
    pm = np.array([r[2] for r in rows])
    pr = np.array([r[1] for r in rows])
    print(f"\n  RMS(realised overall)                      {rms(ov):.6f}")
    print(f"  RMS(realised overall - realised pds+mse/6)  {rms(ov - pm):.6f}"
          "   <- the unpriced four")
    print(f"  RMS(realised overall - predicted pds+mse/6) {rms(ov - pr):.6f}")
    print(f"  sign agreement, predicted vs realised overall: "
          f"{int(np.sum(np.sign(pr) == np.sign(ov)))}/{len(ov)}")

    print("\n" + "=" * 78)
    print("PART 5 — robustness: drop one FOLD from the predictor and refit")
    print("=" * 78)
    print(f"{'member':<8} {'all 4':>18} " + " ".join(f"{'-' + f:>18}" for f in FOLDS))
    for m in ("pds", "mse", "nmae", "reach"):
        exc = ("a11",) if m == "mse" else ()
        cells = []
        for drop in (None,) + FOLDS:
            t = member_table(m, exclude=exc, folds=tuple(f for f in FOLDS if f != drop))
            cells.append(f"{t['loo_rms']:.5f} {t['sign_hits']}/{len(t['arms'])}")
        print(f"{m:<8} " + " ".join(f"{c:>18}" for c in cells))
    print("  (`pds` is the only member whose sign agreement survives every fold drop.)")

    print("\n" + "=" * 78)
    print("PART 6 — K sensitivity of the DE-side pools (K is the one free knob)")
    print("=" * 78)
    for m, base in (("nmae", "nmae_k"), ("reach", "reach_k")):
        cells = []
        for K in (120, 240, 480):
            CHANNEL[m] = (f"{base}{K}", CHANNEL[m][1])
            t = member_table(m)
            cells.append(f"K={K}: LOO {t['loo_rms']:.5f} null {t['null_rms']:.5f} "
                         f"sgn {t['sign_hits']}/6")
        CHANNEL[m] = (f"{base}240", CHANNEL[m][1])
        print(f"  {m:<6} " + " | ".join(cells))
    print("  Both members fail their predict-zero null at every K, so the failure is not a"
          "\n  bad choice of pool size.")

    print("\n" + "=" * 78)
    print("PART 7 — what a slot gate would actually have to clear")
    print("=" * 78)
    tp = tables["pds"]
    sigma = tp["loo_rms"]
    print(f"  priced channel     calibrated d_pds/6, LOO sigma {sigma:.6f} scaled pds "
          f"= {sigma/6:.6f} on overall")
    other5 = np.array([board.delta(a)["overall"] - board.delta(a)["pds"] / 6.0
                       for a in GEOMETRY])
    best = max(abs(tp["k_in"] * tp["x"][i]) for i in range(len(tp["arms"]))) / 6.0
    print(f"  unpriced channels  realised (overall - pds/6) over the 6 geometry arms: "
          f"RMS {rms(other5):.6f}")
    print(f"                     values " + " ".join(f"{v:+.5f}" for v in other5))
    print(f"  largest priced move this validator has ever produced: {best:+.6f} on overall")
    print(f"  -> unpriced/priced ratio {rms(other5)/best:.1f}x")

    amp = ("a11", "a12", "a14", "a19")
    shape = ("a15", "a23")
    for label, grp in (("amplitude-changing", amp), ("shape-only", shape)):
        v = np.array([board.delta(a)["overall"] - board.delta(a)["pds"] / 6.0 for a in grp])
        print(f"     {label:<20} unpriced RMS {rms(v):.6f}  ({', '.join(grp)})")

    print("\n  a23 veto test — could the gate have stopped the arm that lost us a slot?")
    toll = float(np.mean([board.delta(a)["reach"] for a in ("a21", "a22")]))
    x23 = tp["x"][tp["arms"].index("a23")]
    k_no23 = tp["ks"][tp["arms"].index("a23")]
    pred23 = k_no23 * x23 / 6.0
    print(f"    emitted-profile reach toll from a21+a22 ONLY (a23 held out): {toll:+.6f}")
    print(f"    PLOLO d_pds for a23, scale fitted WITHOUT a23 (k={k_no23:.3f}): "
          f"{k_no23*x23:+.6f} -> {pred23:+.6f} on overall")
    print(f"    composite prediction {pred23 + toll/6.0:+.6f}   realised "
          f"{board.delta('a23')['overall']:+.6f}   -> "
          f"{'VETOED (correct)' if pred23 + toll/6.0 < 0 else 'PASSED (wrong)'}")

    out = {m: dict(arms=t["arms"], predicted=list(t["x"]), fold_sd=list(t["sd"]),
                   realised=list(t["y"]), loo_pred=list(t["yhat"]),
                   loo_resid=list(t["resid"]), k_in=t["k_in"], loo_rms=t["loo_rms"],
                   null_rms=t["null_rms"], sign_hits=t["sign_hits"], perm_p=t["perm_p"])
           for m, t in tables.items()}
    (WORK / "retrodict.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {WORK / 'retrodict.json'}")


if __name__ == "__main__":
    main()
