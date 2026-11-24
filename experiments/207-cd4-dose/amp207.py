"""207 — does the dose move `A`?  The one assumption the `mse` price rests on.

096/197's identity is `MSE_c - 1 = (A_c + Kc_c - min(claim_c, B_c) - 2X) / SigmaD`.  Pricing a
dose change purely through `X` (i.e. through the cosine) is only legitimate if `A` is unchanged.

170 §2 point 4 is the reason this cannot be asserted: `apply_v9_whitening` preserves the MEDIAN
row norm while `A` is a weighted MEAN OF SQUARES, and under a gamma sweep the two came apart by
14%.  `add_post_whiten_residual` preserves the same median, so the same trap applies here.

Measured below, on the DELIVERED object in all three challenge frames:
  * median delivered row norm  (what the preserve scalar pins),
  * mean of SQUARED delivered row norms  (what `A` actually is),
  * and `B` = mean of squared row norms about the across-target mean, whose ratio to `A` is the
    `eps` of the identity.

    uv run python amp207.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from plolo import E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for, deliver, solve_fraction  # noqa: E402
from sweep207 import ALPHA, GAMMA, LAM  # noqa: E402

DOSES = [0.0, 0.05, 0.10, 0.20, 0.28, 0.34, 0.40, 0.55, 0.80]


def main() -> None:
    acc: dict[float, list] = {}
    for held in FOLDS:
        fold = build_fold(held, limit_gib=7.0)
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        op = Operator(PROD_OPS / "whitener_union_pruned.npz"
                      if FOLDS[held]["op_kind"] == "production" else fold["op_path"])
        w = fold["wcols"]
        whitened = op.apply(pre * ALPHA, w, GAMMA, LAM)
        for v in DOSES:
            if v > 0:
                frac, _ = solve_fraction(whitened, cd4, v)
                fin, _ = E062.add_post_whiten_residual(whitened, cd4, frac, w)
            else:
                fin = whitened.copy()
            fin[:, fold["panel_positions"]] = 0.0
            per = []
            for ctx in "ABC":
                b = plolo.challenge_bctrl(ctx)
                e = np.expm1(b)
                g = (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (b > 0)
                D = deliver(fin, b, g, "minnorm")
                n = np.linalg.norm(D, axis=1)
                A = float(np.mean(n ** 2))
                mu = D.mean(0)
                B = float(np.mean(np.linalg.norm(D - mu, axis=1) ** 2))
                per.append((float(np.median(n)), A, B))
            acc.setdefault(v, []).append(np.mean(per, axis=0))
        print(f"[{held}] done", flush=True)
        del op

    base = np.mean(acc[0.10], axis=0)
    out = {}
    print(f"\n{'f':>7}{'median |row|':>14}{'ratio':>9}{'A = mean|row|^2':>18}{'ratio':>9}"
          f"{'B':>14}{'ratio':>9}{'eps = A/B-1':>13}{'implied Δraw mse':>18}")
    for v in DOSES:
        m, A, B = np.mean(acc[v], axis=0)
        # the identity's A-channel contribution, holding X fixed: d raw = dA / SigmaD
        d_raw = (A - base[1]) / base[1] * 0.2331583 / 0.229141
        out[f"{v:g}"] = dict(median=m, A=A, B=B, ratio_A=A / base[1],
                             implied_d_raw_mse=d_raw)
        print(f"{v:>7.2f}{m:>14.6f}{m/base[0]:>9.5f}{A:>18.6f}{A/base[1]:>9.5f}"
              f"{B:>14.6f}{B/base[2]:>9.5f}{A/B-1:>13.6f}{d_raw:>+18.6f}")
    (WORK / "amp207.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {WORK/'amp207.json'}")


if __name__ == "__main__":
    main()
