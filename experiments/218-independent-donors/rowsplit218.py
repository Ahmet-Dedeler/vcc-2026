"""218 — where does a partial-coverage arm's gain actually live?

A residual that touches 61 of 290 rows must improve those 61 and leave the other 229 alone
(up to the global preserve rescale, which is scale-only in effect space but NOT in delivered
space, because the renderer is nonlinear).  If the untouched rows move materially, the arm is
not doing what it says.

Prints, per fold and per dose, the per-target cosine to truth split by carrying / not, and
the implied per-carrying-row relative gain -- the number that has to be biologically
plausible for the aggregate number to be believed.

    uv run python rowsplit218.py
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
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   deliver, solve_fraction)
from resid218 import ALPHA, GAMMA, LAM, source_residual  # noqa: E402
from compose218 import merged, restrict_rows, solve_fraction_carrying  # noqa: E402

CD4_F = 0.34
SRC = "gse295214_h1"
G = (0.17, 0.28, 0.45)


def per_target_cos(P, T):
    a = np.linalg.norm(P, axis=1)
    b = np.linalg.norm(T, axis=1)
    out = np.zeros(P.shape[0])
    ok = (a > 0) & (b > 0)
    out[ok] = (P[ok] * T[ok]).sum(1) / (a[ok] * b[ok])
    return out


def main() -> None:
    rows = []
    for held in FOLDS:
        fold = build_fold(held)
        targets = fold["targets"]
        sc = fold["scored"]
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(targets)
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        whitened = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)
        frac_cd4, _ = solve_fraction(whitened, cd4, CD4_F)
        src = source_residual(SRC, targets, "gene_norm")
        carry = np.asarray(src.target_condition_count) > 0
        rowcd4 = restrict_rows(cd4, src)

        b = plolo.challenge_bctrl("A")
        expr = np.expm1(b)
        gate = (1e6 * expr / expr.sum() > plolo.MIN_CPM) & (b > 0)

        def cos_of(res_second, frac_x):
            final, _ = E062.add_post_whiten_residual(
                whitened, merged(cd4, frac_cd4, res_second, frac_x), 1.0, fold["wcols"])
            final[:, fold["panel_positions"]] = 0.0
            return per_target_cos(deliver(final, b, gate, "minnorm")[:, sc], truth)

        base = cos_of(None, 0.0)
        print(f"[{held}] {len(targets)} targets, {SRC} carries {int(carry.sum())} "
              f"({carry.mean():.1%})")
        print(f"  {'g':>5s} {'arm':>8s} {'carrying cos':>13s} {'rel%':>8s} "
              f"{'NOT carrying':>13s} {'rel%':>8s}  {'macro rel%':>10s}")
        print(f"  {'':5s} {'base':>8s} {base[carry].mean():13.6f} {'':8s} "
              f"{base[~carry].mean():13.6f}")
        for g in G:
            fx, _ = solve_fraction_carrying(whitened, src, g)
            for tag, res in (("real", src), ("rowCD4", rowcd4)):
                c = cos_of(res, fx)
                r_in = c[carry].mean() / base[carry].mean() - 1
                r_out = c[~carry].mean() / base[~carry].mean() - 1
                r_all = c.mean() / base.mean() - 1
                rows.append(dict(fold=held, g=g, arm=tag, n_carry=int(carry.sum()),
                                 cos_in=float(c[carry].mean()),
                                 cos_out=float(c[~carry].mean()),
                                 base_in=float(base[carry].mean()),
                                 base_out=float(base[~carry].mean()),
                                 rel_in=float(r_in), rel_out=float(r_out),
                                 rel_macro=float(r_all)))
                print(f"  {g:5.2f} {tag:>8s} {c[carry].mean():13.6f} {100*r_in:+8.2f} "
                      f"{c[~carry].mean():13.6f} {100*r_out:+8.2f}  {100*r_all:+10.2f}")
        del op, whitened
    (WORK / "rowsplit218.json").write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {WORK/'rowsplit218.json'}")


if __name__ == "__main__":
    main()
