"""229 — F1 (target-shuffled) and F3 (row-matched CD4) on the constructions 229 can BUILD.

225 ran both controls on its own frame -- a second residual riding on top of `h1union`.  The
package a registry-driven build emits is a DIFFERENT object (one union, one dose), so its
controls have to be run on it and not inherited.  Same two instruments, verbatim:

  F1  `resid218.shuffled` -- 218's derangement.  Geometry, per-row norms, support pattern and
      carrying count preserved exactly; only WHICH TARGET each row belongs to moves.  Two
      seeds, 20726 and 771103, as in 225.
  F3  `compose218.restrict_rows(cd4, union)` -- CD4's own information on exactly the union's
      rows, sized by the same carrying-row rule at the same dose.  "Does the source beat the
      free alternative on the rows it covers?"  Scored as a DIFFERENCE (225's correction to
      223's ratio convention).

Both are scored against the INCUMBENT `h1union`, which is what the slot is spent against.

    uv run python falsif229.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import price229 as P  # noqa: E402  (imports plolo with 229's own WORK, registers the banks)

from price229 import (CANCER, CD4_F, GRID, H1_G, H1_PAIR, STACK4,  # noqa: E402
                      STACK6, WORK, price)
from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold,  # noqa: E402
                   cd4_for, solve_fraction)
from resid218 import ALPHA, GAMMA, LAM, score_effects, shuffled  # noqa: E402
from compose218 import merged, restrict_rows, solve_fraction_carrying  # noqa: E402

SEEDS = (20726, 771103)
DOSES = (0.28, 0.35)


def main() -> None:
    t0 = time.time()
    MOD, builder_path = P.load_builder()
    rows = []
    for held in CANCER:
        fold = build_fold(held)
        targets = list(fold["targets"])
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(targets)
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        whitened = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)
        del op
        frac_cd4, _ = solve_fraction(whitened, cd4, CD4_F)
        h1, _ = P.production_union(MOD, H1_PAIR, targets)
        frac_h1, _ = solve_fraction_carrying(whitened, h1, H1_G)
        base = merged(cd4, frac_cd4, h1, frac_h1)

        def finish(res):
            f, _ = E062.add_post_whiten_residual(whitened, res, 1.0, fold["wcols"])
            f[:, fold["panel_positions"]] = 0.0
            return f

        f = finish(base)
        b = score_effects(f, fold, truth)
        del f

        for arm_tag, names in (("ARM4", STACK4), ("ARM6", STACK6)):
            R, _ = P.production_union(MOD, names, targets)
            variants = {"": R, "CD4ROWS": restrict_rows(cd4, R)}
            for seed in SEEDS:
                variants[f"SHUF{seed}"] = shuffled(R, seed)
            for g in DOSES:
                for vtag, res in variants.items():
                    fr, _ = solve_fraction_carrying(whitened, res, g)
                    m = merged(cd4, frac_cd4, res, fr)
                    ff = finish(m)
                    e = score_effects(ff, fold, truth)
                    del ff, m
                    rows.append(dict(
                        fold=held, arm=arm_tag, variant=vtag, g=g, fraction=float(fr),
                        d_cos=e["cos"] - b["cos"],
                        rel_cos_pct=100.0 * (e["cos"] - b["cos"]) / max(abs(b["cos"]), 1e-300),
                        d_pds=e["pds"] - b["pds"], d_qw=e["a15_sgn_w"] - b["a15_sgn_w"],
                        carrying=int(np.sum(np.asarray(res.target_condition_count) > 0)),
                    ))
                    print(f"[{held:>14} {arm_tag} {vtag:<10} g={g}] rel {rows[-1]['rel_cos_pct']:+7.3f}%  "
                          f"d_pds {rows[-1]['d_pds']:+.6f}", flush=True)
            del R, variants
        del whitened
    (WORK / "falsif229.json").write_text(json.dumps(rows, indent=1))
    print(f"\n{len(rows)} rows, {time.time() - t0:.0f}s, peak {P.peak_gib():.2f} GiB")


if __name__ == "__main__":
    main()
