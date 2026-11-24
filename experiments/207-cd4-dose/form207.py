"""207 Phase 2b — decompose the FORM into its two independent ingredients.

The library form and the residual form differ in TWO ways at once, so neither can be
attributed:

  (a) the library member is WHITENED with the rest; the residual is added un-whitened;
  (b) the library member REPLACES consensus mass (denominator `votes + w`); the residual
      ADDS mass and the sum is then rescaled back to the base's median row norm.

`wres` isolates (a): whiten CD4 with the SAME operator, then add it post-whitening at the
same `f` solve the residual uses.  If `wres` behaves like the library arm and unlike the
residual, whitening the CD4 rows is the harmful step; if it behaves like the residual,
the harmful step is the mass replacement.

    uv run python form207.py
"""
from __future__ import annotations

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

from plolo import E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for, solve_fraction  # noqa: E402
from sweep207 import ALPHA, GAMMA, LAM, sha8, score_effects  # noqa: E402

DOSES = [0.05, 0.10, 0.20, 0.28, 0.40, 0.55, 0.80]


def main() -> None:
    out_path = WORK / "form207.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    for held in FOLDS:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        w = fold["wcols"]
        whitened = op.apply(pre * ALPHA, w, GAMMA, LAM)
        # CD4 pushed through the SAME operator, with the same on-axis preserve.
        cd4w = op.apply(np.asarray(cd4.mean_effect, dtype=np.float64), w, GAMMA, LAM)
        cd4w[np.asarray(cd4.gene_votes) == 0] = 0.0
        cd4w[:, fold["panel_positions"]] = 0.0
        wcd4 = E062.CD4Prepared(mode="cd4_whitened", mean_effect=cd4w,
                                gene_votes=np.asarray(cd4.gene_votes).copy(),
                                target_condition_count=np.asarray(
                                    cd4.target_condition_count).copy(),
                                condition_details=cd4.condition_details)
        for v in DOSES:
            arm = f"wres{v:g}"
            if (held, arm) in done:
                continue
            t0 = time.time()
            frac, per_unit = solve_fraction(whitened, wcd4, v)
            final, preserve = E062.add_post_whiten_residual(whitened, wcd4, frac, w)
            final[:, fold["panel_positions"]] = 0.0
            e = score_effects(final, fold, truth)
            rows.append(dict(fold=held, arm=arm, n_targets=len(fold["targets"]),
                             note=f"WHITENED CD4 added post-whitening at f={v:g}", **e,
                             audit=dict(kind="wres", fraction=float(frac),
                                        f_per_unit=float(per_unit),
                                        effective_f=float(frac * per_unit),
                                        preserve_scale=float(
                                            preserve["post_add_global_preserve_scale"]),
                                        operator=dict(path=str(op_path.relative_to(ROOT)),
                                                      sha8=sha8(op_path))),
                             seconds=round(time.time() - t0, 1)))
            out_path.write_text(json.dumps(rows, indent=2))
            print(f"[{held}/{arm}] frac {frac:.5f} pds {e['pds']:.6f} "
                  f"cos_agg {e['cos_agg']:.6f} q_w {e['a15_sgn_w']:.6f}", flush=True)
        del op


if __name__ == "__main__":
    main()
