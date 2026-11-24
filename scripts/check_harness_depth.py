"""Acceptance test for a scoring harness: are the two pools at the same sequencing depth?

`build_orion_harness.py` shipped a greedy first-come-first-served control fill until
2026-08-27, which drew all 18,400 controls from the first ~20 (deep) batches while
perturbation cells accumulated across all 109. The result was control cells 38% deeper
than perturbation cells in 132/132 perturbations. Because the Wilcoxon is a rank test and
CPM normalisation does not fix the per-cell zero fraction, the reference's DE calls came
out 90.8% negative against ~0.50 on every clean dataset -- and every DE number measured on
that file is suspect. See docs/HARNESS-DEPTH-CONFOUND.md.

Run this on any harness before trusting a DE result from it.

    uv run python scripts/check_harness_depth.py data/harness/orion-hct116.h5ad
"""
from __future__ import annotations

import sys

import anndata as ad
import numpy as np

CTRL_LABELS = {"non-targeting", "non_targeting", "ntc", "control", "nt"}
# Judgement calls, set from the reference datasets in docs/HARNESS-DEPTH-CONFOUND.md:
# clean files land at 0.985 (h1-test) to 1.020 (Replogle RPE1); the broken one is 0.725,
# and Norman at 0.885 already shows a graded 0.797 negative fraction.
PASS_LO, PASS_HI = 0.95, 1.05


def main(path: str, pert_col: str = "target_gene") -> int:
    A = ad.read_h5ad(path, backed="r")
    lab = np.asarray(A.obs[pert_col].astype(str).values, dtype="U64")
    X = A.X
    tot = np.zeros(A.shape[0])
    for i in range(0, A.shape[0], 20_000):
        m = X[i : i + 20_000]
        tot[i : i + 20_000] = np.asarray(
            (m.todense() if hasattr(m, "todense") else m).sum(1)
        ).ravel()

    is_ctrl = np.isin(np.char.lower(lab), list(CTRL_LABELS))
    if not is_ctrl.any():
        print(f"no control cells found in {pert_col!r}; labels seen: {set(lab[:50])}")
        return 2
    c_med, p_med = np.median(tot[is_ctrl]), np.median(tot[~is_ctrl])
    ratio = p_med / c_med

    targets = np.unique(lab[~is_ctrl])
    per = np.array([np.median(tot[lab == t]) / c_med for t in targets])

    print(f"{path}")
    print(f"  cells                {A.shape[0]:,} x {A.shape[1]:,}")
    print(f"  control cells        {int(is_ctrl.sum()):,}  median UMI {c_med:,.0f}")
    print(f"  perturbation cells   {int((~is_ctrl).sum()):,}  median UMI {p_med:,.0f}")
    print(f"  pert/control depth   {ratio:.4f}")
    print(
        f"  per-perturbation     median {np.median(per):.4f}  "
        f"IQR [{np.percentile(per, 25):.4f}, {np.percentile(per, 75):.4f}]"
    )
    print(
        f"  shallower than ctrl  {int((per < 1).sum())}/{per.size} "
        f"= {(per < 1).mean() * 100:.1f}%"
    )
    ok = PASS_LO <= ratio <= PASS_HI
    print(f"\n  {'PASS' if ok else 'FAIL'} — depth ratio {ratio:.4f} "
          f"{'within' if ok else 'OUTSIDE'} [{PASS_LO}, {PASS_HI}]")
    if not ok:
        print("  Do NOT read absolute fid / reach / jac / n_real or DE direction off this file.")
        print("  Rebuild with `--sampling reservoir`; see docs/HARNESS-DEPTH-CONFOUND.md.")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    raise SystemExit(main(*sys.argv[1:]))
