"""188 — cache the four PLOLO folds' (prediction, truth) pairs, once.

Everything downstream in this experiment needs the same object: the production a23 effect
matrix built WITHOUT the held library, and that library's own effect matrix as truth, on the
official gene axis, over the real 300-target panel.  162 builds it; 169 re-used it; this
script freezes it so the map fit and the pricing never rebuild a 482 MB operator.

WHAT IS CACHED, AND IN WHICH SPACE
----------------------------------
`pred`  -- 162 `run_arms.effects_for(fold, CONFIGS['a23'], op, cd4)`, i.e.
           consensus x alpha 0.982 -> whitened (gamma 0.65) -> + CD4 at f = 0.10 -> panel
           coordinates zeroed.  This is EXACTLY the object production hands the renderer and
           exactly the object 169's `sharpen` transformed.
`truth` -- `E062.raw_truth(held)` = the held donor's `delta.npy` row minus that donor's mean
           perturbation response, in the donor's own `log1p(TS*p)` space
           (`scripts/build_donor_table.py:9`).  So prediction and truth are the SAME space,
           which is what makes a magnitude map fittable at all.
`b_ctrl`-- the held donor's own control profile on the official axis, needed to convert
           either side to log2 fold change.

No cell-eval2, no build, no submission.  Nothing outside this directory is written.
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
import plolo                                                    # noqa: E402
from plolo import FOLDS, PROD_OPS, Operator, build_fold, cd4_for  # noqa: E402
from run_arms import CONFIGS, effects_for                        # noqa: E402


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    for held in FOLDS:
        out = os.path.join(WORK, f"pair-{held}.npz")
        if os.path.exists(out):
            print(f"[{held}] cached", flush=True)
            continue
        fold = build_fold(held)
        cd4 = cd4_for(fold["targets"])
        op = Operator(PROD_OPS / "whitener_union_pruned.npz"
                      if FOLDS[held]["op_kind"] == "production" else fold["op_path"])
        final, caud = effects_for(fold, CONFIGS["a23"], op, cd4)
        del op
        truth = np.asarray(fold["truth"], dtype=np.float64)
        np.savez_compressed(
            out,
            pred=final.astype(np.float32), truth=truth.astype(np.float32),
            b_ctrl=np.asarray(fold["b_ctrl"], dtype=np.float64),
            gate=np.asarray(fold["gate"], dtype=bool),
            scored=np.asarray(fold["scored"], dtype=np.int64),
            targets=np.asarray(fold["targets"]),
            audit=json.dumps(dict(fold["audit"], **caud)),
        )
        print(f"[{held}] pred {final.shape} truth {truth.shape} "
              f"median |pred| row norm {np.median(np.linalg.norm(final, axis=1)):.4f} "
              f"truth {np.median(np.linalg.norm(truth, axis=1)):.4f}", flush=True)
        del fold, final, truth, cd4

    # the three challenge frames, cached alongside so nothing downstream touches 162's dir
    fr = os.path.join(WORK, "chalframes.npz")
    if not os.path.exists(fr):
        d = {}
        for ctx in "ABC":
            b = plolo.challenge_bctrl(ctx)
            expr = np.expm1(b)
            cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
            d[f"b{ctx}"] = b
            d[f"g{ctx}"] = (cpm > plolo.MIN_CPM) & (b > 0)
        np.savez_compressed(fr, **d)
        print("wrote challenge frames: "
              + ", ".join(f"{c}={int(d['g'+c].sum())} gated" for c in "ABC"))


if __name__ == "__main__":
    sys.exit(main())
