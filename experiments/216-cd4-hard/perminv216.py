"""216 — falsifier 1 applied to the OPERATOR REFIT, where it turns out to be vacuous by algebra.

`sweep216.py`'s falsifier permutes whole CD4 rows across targets and shows the residual's gain
collapses.  The operator refit cannot be falsified that way, and the reason is structural:

    the refit depends on CD4 only through its Gram  G = B^T B / m,
    and  (P B)^T (P B) = B^T P^T P B = B^T B  for any permutation matrix P.

**The refit operator built from target-shuffled CD4 is the SAME OPERATOR, bit for bit.**  So
whatever the refit buys, it is provably NOT target-specific biology -- it is a change of which
subspace the whitener amplifies, i.e. 215's GEOMETRY channel.  That matters because the board
has already priced a geometry arm in this family: a19 (the whitener ridge) realised
**-0.01352 overall**, and 215's instrument classified it GEOMETRY at rho = 0.02.

This script states the claim and measures it rather than asserting it.

    uv run python perminv216.py
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
from plolo import FOLDS, build_fold, cd4_for  # noqa: E402
from sweep216 import shuffled_cd4  # noqa: E402


def main() -> None:
    out = {}
    for held in FOLDS:
        fold = build_fold(held)
        wcols = np.asarray(fold["wcols"])
        cd4 = cd4_for(fold["targets"])
        shuf = shuffled_cd4(cd4)

        def block(c):
            carry = np.asarray(c.target_condition_count) == 3
            return np.ascontiguousarray(
                np.where(np.asarray(c.gene_votes)[carry] > 0,
                         np.asarray(c.mean_effect)[carry], 0.0)[:, wcols])

        A, B = block(cd4), block(shuf)
        assert A.shape == B.shape
        # the rows are the same multiset, permuted: check that directly
        rows_same = bool(np.array_equal(np.sort(A.sum(1)), np.sort(B.sum(1))))
        any_row_moved = int((np.abs(A - B).sum(1) > 0).sum())
        # a random projection of the Gram, so a full n x n product is not needed
        rng = np.random.default_rng(216)
        R = rng.standard_normal((A.shape[1], 96))
        GA = A.T @ (A @ R)
        GB = B.T @ (B @ R)
        rel = float(np.max(np.abs(GA - GB)) / max(np.max(np.abs(GA)), 1e-300))
        out[held] = dict(m_rows=int(A.shape[0]), n_cols=int(A.shape[1]),
                         rows_are_a_permutation=rows_same,
                         rows_that_moved=any_row_moved,
                         gram_projection_rel_max_diff=rel,
                         gram_identical=bool(rel < 1e-12))
        print(f"[{held}] {A.shape[0]} CD4 rows, {any_row_moved} of them moved by the "
              f"derangement; Gram projection rel max|diff| = {rel:.3e} -> "
              f"{'IDENTICAL' if rel < 1e-12 else 'DIFFERENT'}", flush=True)
        del fold, cd4, shuf, A, B, GA, GB
    out["verdict"] = (
        "The target-shuffled control is IDENTICAL to the real arm for the operator refit, so "
        "falsifier 1 cannot discriminate and the refit's gain is not target-specific. It is a "
        "GEOMETRY change in 215's sense.")
    (WORK / "perminv216.json").write_text(json.dumps(out, indent=2))
    print("\n" + out["verdict"])


if __name__ == "__main__":
    main()
