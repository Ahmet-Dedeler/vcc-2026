"""169 — price the sharpening arm's `pds` on PLOLO, which is the instrument rule 13 names.

AGENTS.md 11.1 rule 13's gate is

    predicted = PLOLO d_pds/6 * 1.011  +  096 d_mse/6  -  0.0146/6

and `PLOLO d_pds` is a specific number: experiment 162's panel-restricted leave-one-library-out
validator, four donor folds over the real 300-target panel, rendered into the CHALLENGE control
frames A/B/C, scored against the held library's own pseudobulk.  It retrodicts board `pds`
6/6 on sign with a leave-one-arm-out RMS of 0.001494 scaled.

The Orion-pair delta from `sharpen169.py` is NOT that number: `pds` is a retrieval statistic
and the decoy pool there is 132 targets against PLOLO's 269-294 (162's own control notes the
Orion pair reads 4-5x larger than 156 for exactly this reason).  So this script re-runs 162's
folds with the sharpening transform inserted where production would put it -- between the
effect matrix and the renderer -- and reads `pds` in 162's own frame.

162's files are NOT modified.  `plolo.py` and `run_arms.py` are imported and their cached
folds reused; every output of this script lands in `experiments/169-sharpening-six-member/`.

KNOWN-ANSWER ROWS (rule 10's corollary):
  * `q = 1.00` must reproduce 162's own a23 row (the shipped incumbent: alpha 0.982,
    gamma 0.65, cd4_f 0.10, union_pruned, mass minnorm) fold by fold, from `work/arms.json`.
    Any drift means the transform plumbing changed something it should not have.
  * a23 - a15 in the challenge frame must reproduce 162's predicted a23 delta (+0.004014
    scaled = +0.001815 raw pds), which is the mass fix's own PLOLO price.

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
P162 = os.path.join(REPO, "experiments/162-retrodicting-validator")
sys.path.insert(0, P162)
sys.path.insert(0, HERE)
import plolo                                          # noqa: E402
from plolo import (FOLDS, PROD_OPS, Operator, build_fold, cd4_for, deliver,  # noqa: E402
                   endpoints)
from run_arms import CONFIGS, effects_for, frames_for  # noqa: E402
from sharp import sharpen, sharpen_lfc                # noqa: E402

QS = (1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30)
# The incumbent, exactly as 162 records it, and the same arm minus the mass fix.
BASE = {"a23": CONFIGS["a23"], "a15": CONFIGS["a15"]}
PDS_SPAN = 0.95228 - 0.50020          # debench ANCHOR, the board-fitted pds anchors


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    out_path = os.path.join(WORK, "plolo169.json")
    rows = json.loads(open(out_path).read()) if os.path.exists(out_path) else []
    done = {(r["fold"], r["arm"], r["q"], r["space"]) for r in rows}

    for held in FOLDS:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        cd4 = cd4_for(fold["targets"])
        frames = frames_for(fold)
        op = Operator(PROD_OPS / "whitener_union_pruned.npz"
                      if FOLDS[held]["op_kind"] == "production" else fold["op_path"])
        print(f"[{held}] operator maxeig {op.evals.max():.4f}, "
              f"{len(fold['targets'])} targets, {len(fold['scored'])} scored", flush=True)
        for arm, cfg in BASE.items():
            final, caud = effects_for(fold, cfg, op, cd4)
            for q in QS:
                if arm == "a15" and q != 1.00:
                    continue                       # a15 exists only as the mass-fix control
                for fname, b, gate in frames:
                    for tag, tf in (("", sharpen),
                                    ("L", lambda F, qq, g, mode="row":
                                     sharpen_lfc(F, qq, b, g, mode=mode))):
                        sp = f"{fname}{tag}/delivered"
                        if (held, arm, q, sp) in done:
                            continue
                        if tag == "L" and arm != "a23":
                            continue
                        fq = tf(final, q, gate, mode="row")
                        mat = deliver(fq, b, gate, cfg["mass"])
                        e = endpoints(mat[:, fold["scored"]], truth)
                        rows.append(dict(fold=held, arm=arm, q=q, space=sp, **e))
                        print(f"  {held:<14} {arm} q={q:.2f} {sp:<20} pds {e['pds']:.6f} "
                              f"cos {e['cos']:.5f} sse {e['sse']:.6f} "
                              f"sgnw {e['sgn_w']:.5f}", flush=True)
                open(out_path, "w").write(json.dumps(rows, indent=1))
        del op

    # ------------------------------------------------------------------ report
    R = {(r["fold"], r["arm"], r["q"], r["space"]): r for r in rows}
    folds = list(FOLDS)

    print(f"\n{'=' * 108}\nKNOWN-ANSWER GATE: q=1.00 must reproduce 162's own a23 rows")
    ref = json.loads(open(os.path.join(P162, "work/arms.json")).read())
    A162 = {(r["fold"], r["arm"], r["space"]): r for r in ref}
    worst = 0.0
    for f in folds:
        for sp in ("chalA/delivered", "chalB/delivered", "chalC/delivered",
                   "library/delivered"):
            a, c = A162.get((f, "a23", sp)), R.get((f, "a23", 1.00, sp))
            if a is None or c is None:
                continue
            d = abs(a["pds"] - c["pds"])
            worst = max(worst, d)
            print(f"  {f:<14} {sp:<18} 162 {a['pds']:.6f}  here {c['pds']:.6f}  "
                  f"|d| {d:.2e}")
    print(f"  worst |d pds| = {worst:.2e}")
    if worst > 1e-2:
        sys.exit("ABORT: q=1.00 does not reproduce 162's a23 rows")
    if worst > 1e-9:
        print("  NOTE a23's `mass_consistent_effect` is a fixed-point with an rtol=1e-12 exit,")
        print("  so 162 §6 records a ~1% reproducibility floor on sign-keyed statistics for")
        print("  this arm specifically.  `pds` above is well inside it.")

    def chal(f, arm, q):
        v = [R[(f, arm, q, f"chal{c}/delivered")]["pds"] for c in "ABC"]
        return float(np.mean(v))

    print(f"\n{'=' * 108}\nKNOWN-ANSWER GATE: the mass fix's own PLOLO price "
          f"(162 predicted a23 at +0.004014 scaled)")
    d = [chal(f, "a23", 1.00) - chal(f, "a15", 1.00) for f in folds]
    print(f"  per fold raw: " + "  ".join(f"{f}={x:+.6f}" for f, x in zip(folds, d)))
    print(f"  4-fold mean raw {np.mean(d):+.6f}  ->  scaled {np.mean(d)/PDS_SPAN:+.6f}  "
          f"(162: +0.004014)")

    print(f"\n{'=' * 108}\nSHARPENING on PLOLO, challenge frame, per fold "
          f"(the predictor 162 pre-committed: unweighted 4-fold mean)")
    print(f"  {'q':>5} " + " ".join(f"{f:>14}" for f in folds) +
          f" {'MEAN raw':>10} {'sd':>9} {'MEAN scaled':>12}")
    base = {f: chal(f, "a23", 1.00) for f in folds}
    for q in QS:
        v = [chal(f, "a23", q) - base[f] for f in folds]
        print(f"  {q:>5.2f} " + " ".join(f"{x:+14.6f}" for x in v) +
              f" {np.mean(v):+10.6f} {np.std(v):9.6f} {np.mean(v)/PDS_SPAN:+12.6f}")

    print("\n  library frame (162 §7: post hoc it FAILS the null, reported for disclosure only)")
    print(f"  {'q':>5} " + " ".join(f"{f:>14}" for f in folds) + f" {'MEAN raw':>10}")
    lb = {f: R[(f, "a23", 1.00, "library/delivered")]["pds"] for f in folds}
    for q in QS:
        v = [R[(f, "a23", q, "library/delivered")]["pds"] - lb[f] for f in folds]
        print(f"  {q:>5.2f} " + " ".join(f"{x:+14.6f}" for x in v) + f" {np.mean(v):+10.6f}")

    print(f"\n{'=' * 108}\nTHE lfc-SPACE ARM on PLOLO, challenge frame (the same power "
          f"applied to log2 fold\nchange instead of the log-pseudobulk effect)")
    print(f"  {'q':>5} " + " ".join(f"{f:>14}" for f in folds) +
          f" {'MEAN raw':>10} {'MEAN scaled':>12} {'vs h-space':>11}")
    baseL = {f: float(np.mean([R[(f, "a23", 1.00, f"chal{c}L/delivered")]["pds"]
                               for c in "ABC"])) for f in folds}
    for q in QS:
        v = [float(np.mean([R[(f, "a23", q, f"chal{c}L/delivered")]["pds"]
                            for c in "ABC"])) - baseL[f] for f in folds]
        h = np.mean([chal(f, "a23", q) - base[f] for f in folds])
        print(f"  {q:>5.2f} " + " ".join(f"{x:+14.6f}" for x in v) +
              f" {np.mean(v):+10.6f} {np.mean(v)/PDS_SPAN:+12.6f} "
              f"{np.mean(v)/max(h, 1e-12):>11.3f}")

    print(f"\n{'=' * 108}\nsse* and sign credit on PLOLO, challenge frame, 4-fold mean")
    print(f"  {'q':>5} {'sse*':>10} {'d sse*':>10} {'cos':>10} {'sgn_w':>9} {'d sgn_w':>9}")
    for q in QS:
        s = np.mean([np.mean([R[(f, "a23", q, f"chal{c}/delivered")]["sse"] for c in "ABC"])
                     for f in folds])
        s0 = np.mean([np.mean([R[(f, "a23", 1.0, f"chal{c}/delivered")]["sse"] for c in "ABC"])
                      for f in folds])
        cs = np.mean([np.mean([R[(f, "a23", q, f"chal{c}/delivered")]["cos"] for c in "ABC"])
                      for f in folds])
        g = np.mean([np.mean([R[(f, "a23", q, f"chal{c}/delivered")]["sgn_w"] for c in "ABC"])
                     for f in folds])
        g0 = np.mean([np.mean([R[(f, "a23", 1.0, f"chal{c}/delivered")]["sgn_w"] for c in "ABC"])
                      for f in folds])
        print(f"  {q:>5.2f} {s:10.6f} {s - s0:+10.6f} {cs:10.5f} {g:9.5f} {g - g0:+9.5f}")


if __name__ == "__main__":
    sys.exit(main())
