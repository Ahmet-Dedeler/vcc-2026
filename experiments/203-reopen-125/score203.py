"""203 Phase 1 / 4 — price 125's estimator on the endpoints the competition reads.

Every arm is pushed through the EXACT a15 production chain that 162 validated against
12 board points:

    effect  ->  x alpha  ->  whiten(gamma=0.65, lam=0)  ->  + CD4 residual at f=0.10
            ->  zero the 300 panel columns  ->  render  ->  read back  ->  score

so the only thing that differs between arms is the estimator that produced `effect`.

ARMS
  a15            162's cached `pre` (production `build_baseline`: equal mean of the raw
                 source effects).  KNOWN-ANSWER ROW -- must reproduce `162/work/arms.json`.
  eq122          122's `equal_mean_agreement` on the production axis (row-normalised mean
                 x output scale).  This is 125's own comparator, so `post125 - eq122` is
                 the estimator delta 125 measured, transplanted to the scored axis.
  post125        125's hierarchical posterior, nested on 125's frozen 5-fold target split.
  post125@a*     the same, at its own re-solved amplitude (mandate: "report the arm at its
                 own re-solved alpha").  alpha* matches the delivered median row norm.
  shrink6.78     a15's own effect divided by 6.78 -- the pure-rescale falsifier, in
                 production units.
  eb_*           the three unselected posterior settings, as a sensitivity.

FRAME: 162's pre-committed instrument -- challenge frame (mean of contexts A/B/C,
delivered), unweighted mean over folds.  `intended` and `library/delivered` are recorded
too but are post-hoc.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "src"))

import est203  # noqa: E402
import plolo  # noqa: E402
from plolo import (FOLDS, PROD_OPS, Operator, build_fold, cd4_for, deliver,  # noqa: E402
                   endpoints, solve_fraction)

OUT = HERE / "work"
OUT.mkdir(exist_ok=True)
A15 = dict(alpha=0.982, gamma=0.65, lam=0.0, cd4_f=0.10, mass=None)
SHRINK = 6.7801


def median_row_norm(M, wcols):
    return float(np.median(np.linalg.norm(np.asarray(M)[:, wcols], axis=1))) or 1.0


def run_chain(fold, eff, cd4, op, alpha):
    pre = np.asarray(eff, dtype=np.float64) * alpha
    wh = op.apply(pre, fold["wcols"], A15["gamma"], A15["lam"])
    frac, per_unit = solve_fraction(wh, cd4, A15["cd4_f"])
    final, _ = plolo.E062.add_post_whiten_residual(wh, cd4, frac, fold["wcols"])
    final[:, fold["panel_positions"]] = 0.0
    return final


def score_all(fold, final, truth, frames):
    rows = {"intended": endpoints(final[:, fold["scored"]], truth)}
    for fname, b, gate in frames:
        mat = deliver(final, b, gate, A15["mass"])
        rows[f"{fname}/delivered"] = endpoints(mat[:, fold["scored"]], truth)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="*",
                    default=["k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"])
    ap.add_argument("--out", default=str(OUT / "arms203.json"))
    args = ap.parse_args()
    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    for held in args.folds:
        fold = build_fold(held, limit_gib=7.0)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        cd4 = cd4_for(fold["targets"])
        frames = [("library", fold["b_ctrl"], fold["gate"])]
        for ctx in ("A", "B", "C"):
            b = plolo.challenge_bctrl(ctx)
            e = np.expm1(b)
            cpm = 1e6 * e / max(float(e.sum()), 1e-12)
            frames.append((f"chal{ctx}", b, (cpm > plolo.MIN_CPM) & (b > 0)))

        op_path = (PROD_OPS / "whitener_union_pruned.npz"
                   if FOLDS[held]["op_kind"] == "production" else fold["op_path"])
        op = Operator(op_path)
        print(f"[{held}] operator {op_path.name} maxeig {op.evals.max():.4f}", flush=True)

        z = np.load(OUT / f"sources-{held}.npz", allow_pickle=True)
        panel = est203.prepare(np.asarray(z["effects"], dtype=np.float64),
                               available=z["available"], support=z["support"])
        targets = [str(t) for t in z["targets"]]
        pre_a15 = np.asarray(fold["pre"], dtype=np.float64)
        mrn_a15 = median_row_norm(pre_a15, fold["wcols"])

        eq_eff, eq_conf, _ = est203.equal_mean_agreement(panel)
        arms = {"a15": (pre_a15, A15["alpha"]),
                "eq122": (eq_eff, A15["alpha"]),
                f"shrink{SHRINK:.2f}": (pre_a15 / SHRINK, A15["alpha"])}
        confs = {"eq122": eq_conf}
        for name in est203.SETTINGS:
            e, c, _ = est203.nested_posterior(panel, targets, est203.SETTINGS[name])
            arms[f"post125:{name}"] = (e, A15["alpha"])
            confs[f"post125:{name}"] = c
            a_star = A15["alpha"] * mrn_a15 / median_row_norm(e, fold["wcols"])
            arms[f"post125:{name}@alpha*"] = (e, a_star)

        for arm, (eff, alpha) in arms.items():
            if (held, arm) in done:
                continue
            mrn = median_row_norm(eff, fold["wcols"])
            final = run_chain(fold, eff, cd4, op, alpha)
            res = score_all(fold, final, truth, frames)
            rec = dict(fold=held, arm=arm, alpha=float(alpha),
                       median_row_norm=mrn, shrink_vs_a15=float(mrn_a15 / mrn),
                       n_targets=len(fold["targets"]),
                       n_scored=int(len(fold["scored"])),
                       delivered_median_norm=float(median_row_norm(final, fold["wcols"])),
                       spaces={k: v for k, v in res.items()})
            rows.append(rec)
            ch = {k: float(np.mean([res[f"chal{c}/delivered"][k] for c in "ABC"]))
                  for k in ("pds", "cos", "sgn_w", "sgn", "rr", "sse")}
            print(f"  {held:<14} {arm:<28} alpha {alpha:6.4f} shrink {mrn_a15/mrn:6.3f}x  "
                  f"CHAL pds {ch['pds']:.6f} cos {ch['cos']:.5f} q_w {ch['sgn_w']:.5f}",
                  flush=True)
            out_path.write_text(json.dumps(rows, indent=1))
        del op
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
