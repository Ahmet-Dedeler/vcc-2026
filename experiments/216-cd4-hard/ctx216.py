"""216 — INBOX-3: the dose curve PER CONTEXT, not just pooled.

`sweep216.py` averaged the three challenge frames inside `score_effects` and only kept the
mean, so the per-context curves INBOX-3 asks for are not recoverable from `work/sweep216.json`.
This re-runs the same grid keeping A, B and C separate.  Nothing else differs -- same folds,
same operator, same alpha, same solve -- so the pooled column here must reproduce
`sweep216.json` exactly, and that is the known-answer row.

Deliverables INBOX-3 names:
  1. each context's own cosine maximiser, and the SPREAD between the three;
  2. each context's `pds` cost curve;
  3. whether A (JURKAT, the T-lineage context) has the highest maximiser;
  4. the loss from shipping the pooled optimum into a context whose optimum differs.

    uv run python ctx216.py
"""
from __future__ import annotations

import gc
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   deliver, endpoints, solve_fraction)
from sweep216 import ALPHA, DOSE_F, GAMMA, LAM  # noqa: E402


def main() -> None:
    t0 = time.time()
    out_path = WORK / "ctx216.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["f"]) for r in rows}

    for held in FOLDS:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        whitened = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)
        del op
        gc.collect()

        frames = {}
        for ctx in ("A", "B", "C"):
            b = plolo.challenge_bctrl(ctx)
            expr = np.expm1(b)
            cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
            frames[ctx] = (b, (cpm > plolo.MIN_CPM) & (b > 0))

        for f_val in DOSE_F:
            if (held, f_val) in done:
                continue
            if f_val > 0:
                frac, _ = solve_fraction(whitened, cd4, f_val)
                final, _ = E062.add_post_whiten_residual(
                    whitened, cd4, frac, fold["wcols"])
            else:
                frac, final = 0.0, whitened.copy()
            final[:, fold["panel_positions"]] = 0.0
            rec = dict(fold=held, f=f_val, fraction=float(frac))
            for ctx, (b, gate) in frames.items():
                for tag, mass in (("a23", "minnorm"), ("a15", None)):
                    e = endpoints(deliver(final, b, gate, mass)[:, fold["scored"]], truth)
                    rec[f"{tag}_{ctx}_cos_agg"] = float(np.sqrt(max(0.0, 1.0 - e["sse"])))
                    rec[f"{tag}_{ctx}_cos"] = float(e["cos"])
                    rec[f"{tag}_{ctx}_pds"] = float(e["pds"])
                    rec[f"{tag}_{ctx}_sgn_w"] = float(e["sgn_w"])
            for tag in ("a23", "a15"):
                for k in ("cos_agg", "cos", "pds", "sgn_w"):
                    rec[f"{tag}_{k}"] = float(np.mean(
                        [rec[f"{tag}_{c}_{k}"] for c in "ABC"]))
            rows.append(rec)
            out_path.write_text(json.dumps(rows))
            print(f"[{held}/f={f_val:g}] A {rec['a23_A_cos_agg']:.6f} "
                  f"B {rec['a23_B_cos_agg']:.6f} C {rec['a23_C_cos_agg']:.6f} "
                  f"pooled {rec['a23_cos_agg']:.6f} ({time.time()-t0:.0f}s)", flush=True)
            del final
        del whitened, pre, cd4, fold, truth
        gc.collect()
    print(f"done in {time.time()-t0:.0f}s, {len(rows)} rows")


if __name__ == "__main__":
    main()
