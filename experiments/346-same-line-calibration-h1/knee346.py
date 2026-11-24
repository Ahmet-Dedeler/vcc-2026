#!/usr/bin/env python3
"""346 step 10 — build the beta=0.25 arm, the knee of the design-stage trade-off curve.

POST-HOC. The registration priced beta in {0.5, 1.0}. Step 9's curve shows the sign gain is
already 39% of its maximum at beta=0.25 while retrieval — which step 8 measured as the
delivered `pds` member at slope 1.01 — is still within 0.001 of cap15. This emits that one
extra point so the trade-off has a measurement at its best design-stage location. Its interval
is NOT corrected for having been chosen after seeing the curve, and RESULTS.md says so.
Same norm-preserving integration rule, same covered axis, same predictions. Writes to work/.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
sys.path.insert(0, str(HERE))
import effects346 as X   # noqa: E402

BETA = 0.25
NAME = "same_ridge_b025"


def main() -> None:
    with np.load(BASE / "work/effects.npz", allow_pickle=False) as z:
        genes = [str(g) for g in z["genes"]]
        targets = [str(t) for t in z["targets"]]
        cap15 = np.asarray(z["cap15"], dtype=np.float64)
    with np.load(WORK / "predictions.npz", allow_pickle=False) as z:
        covered = [str(g) for g in z["genes"]]
        test = [str(t) for t in z["test_targets"]]
        prediction = np.asarray(z["same_line_full"], dtype=np.float64)
    gi = {g: i for i, g in enumerate(genes)}
    cols = np.array([gi[g] for g in covered], dtype=np.int64)
    ti = {t: i for i, t in enumerate(targets)}
    rows = np.array([ti[t] for t in test], dtype=np.int64)
    new = cap15.copy()
    for j, i in enumerate(rows):
        new[i, cols] = X.blend(cap15[i, cols], prediction[j], BETA)
    outside = np.ones(len(genes), dtype=bool)
    outside[cols] = False
    checks = {
        "beta": BETA,
        "outside_covered_axis_bit_identical": bool(np.array_equal(new[:, outside],
                                                                  cap15[:, outside])),
        "covered_norm_max_abs_error": float(np.abs(
            np.linalg.norm(new[:, cols], axis=1)
            - np.linalg.norm(cap15[:, cols], axis=1)).max()),
        "rows_bit_identical_to_cap15": int(np.sum(np.all(new == cap15, axis=1))),
        "finite": bool(np.isfinite(new).all()),
        "post_hoc": True,
    }
    if not checks["outside_covered_axis_bit_identical"] or \
            checks["covered_norm_max_abs_error"] > 1e-10 or not checks["finite"]:
        raise AssertionError("integration invariants failed for the knee arm")
    np.savez_compressed(WORK / "effects-b025.npz", genes=np.array(genes),
                        targets=np.array(targets), covered_columns=cols,
                        modified_rows=rows, cap15=cap15, **{NAME: new})
    (WORK / "knee.json").write_text(json.dumps(checks, indent=2, sort_keys=True) + "\n")
    print(json.dumps(checks, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
