"""Build real+predicted DE tables for ANY (harness, prediction) pair, scorer-exact.

experiments/092 hardcoded the Orion HCT116 paths, which is how six experiments ended up
tuned on a depth-confounded file without anyone being able to re-point them. This takes
both paths as arguments.

    uv run python experiments/095-clean-harness/de_tables.py \
        --real data/harness/orion-hct116-thinned.h5ad \
        --pred <pred.h5ad> --out experiments/095-clean-harness/de_thinned.npz
"""
from __future__ import annotations

import argparse
import os
import sys

import anndata as ad
import numpy as np
from scipy.stats import false_discovery_control, norm

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "092-call-set-truncation")
)
from de_kernel import ctrl_tiesums, wilcoxon_z  # noqa: E402

CTRL, PERT_COL, ALPHA, MIN_CPM = "non-targeting", "target_gene", 0.05, 5.0


def norm_log(x: np.ndarray) -> np.ndarray:
    """The scorer's transform: CPM to 1e6 per cell, then log1p."""
    s = x.sum(1)
    s[s == 0] = 1.0
    return np.log1p(x * (1e6 / s[:, None]))


def dense(A, rows):
    m = A.X[rows]
    return np.asarray(m.todense() if hasattr(m, "todense") else m, dtype=np.float64)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", required=True)
    ap.add_argument("--pred", default=None, help="omit to emit reference tables only")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    R = ad.read_h5ad(a.real)
    rlab = np.asarray(R.obs[PERT_COL].astype(str).values, dtype="U64")
    ci = np.flatnonzero(rlab == CTRL)
    C = norm_log(dense(R, ci))
    cpm = np.expm1(C).mean(0)
    keep = cpm > MIN_CPM
    genes = np.asarray([str(g) for g in R.var_names])[keep]
    cs = np.ascontiguousarray(np.sort(C[:, keep].T, axis=1))
    tie = ctrl_tiesums(cs)
    ref_lin = np.expm1(C[:, keep]).mean(0)
    del C
    perts = np.array(sorted(set(rlab) - {CTRL}))
    print(f"real {a.real}: {len(perts)} perturbations, {ci.size} controls, "
          f"{keep.sum()}/{keep.size} genes past the {MIN_CPM}-CPM gate")

    def table(A, lab, tag):
        G = keep.sum()
        P, Q, L = (np.full((perts.size, G), np.nan) for _ in range(3))
        for i, t in enumerate(perts):
            rows = np.flatnonzero(lab == t)
            if rows.size == 0:
                continue
            B = norm_log(dense(A, rows))
            z = wilcoxon_z(cs, tie, np.ascontiguousarray(B[:, keep].T))
            p = 2 * norm.sf(np.abs(z))
            P[i], Q[i] = p, false_discovery_control(p, method="bh")
            L[i] = np.log2((np.expm1(B[:, keep]).mean(0) + 1e-9) / (ref_lin + 1e-9))
            del B
            if (i + 1) % 25 == 0:
                print(f"  {tag} {i+1}/{perts.size}", flush=True)
        return P, Q, L

    out = {"genes": genes, "perts": perts}
    out["real_p"], out["real_padj"], out["real_lfc"] = table(R, rlab, "real")
    del R
    if a.pred:
        A = ad.read_h5ad(a.pred)
        plab = np.asarray(A.obs[PERT_COL].astype(str).values, dtype="U64")
        out["pred_p"], out["pred_padj"], out["pred_lfc"] = table(A, plab, "pred")
    np.savez_compressed(a.out, **out)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
