"""164 §1 — what the SHIPPED tier policy actually puts in each block, on the real panel.

157 §5b measured the zero-effect call artefact with `gated = mean-per-cell-CPM > 5.0` and
`KEEP_TOP = 8000` taken over that set.  Production does not gate there.
`emission_spread_vector(..., min_cpm=10.0)` gates on the **pseudobulk** CPM
`1e6 * expm1(b_ctrl)/sum(expm1(b_ctrl))` at **10.0**, and every gene outside that gate keeps
`tail_spread`.  The scorer's DE universe is a *different* set again
(`filter_gene_min_cpm_cell = 5.0`, mean per-cell CPM).

So the tail block, as the scorer sees it, is

    {production-gated, |effect| rank > 8000}   UNION   {5 < mean-cell-CPM, pseudobulk CPM <= 10}

and the second term is invisible to 157's arithmetic.  This script measures both terms on
`data/challenge/context_{A,B,C}.h5ad`.  If the production gate holds <= 8000 genes the tail
block is *entirely* the second term, which makes it independent of the predicted effect and
therefore exactly computable without a build.

Reads only control cells.  No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import scipy.sparse as sp

TS = 5e4
ROOT = "~/Code/virtual-cell-challenge-2026"


def main() -> None:
    print(f"{'ctx':>4}{'scorer DE (cell-CPM>5)':>24}{'prod gate (pb-CPM>10)':>23}"
          f"{'both':>8}{'DE & not prod-gated':>21}{'prod-gated & not DE':>21}")
    rows = {}
    for ctx in "ABC":
        A = ad.read_h5ad(f"{ROOT}/data/challenge/context_{ctx}.h5ad")
        X = A.X.tocsc().astype(np.float64)
        lib = np.asarray(X.sum(1)).ravel()
        lib[lib == 0] = 1.0
        # scorer's gate: mean over cells of per-cell CPM  (de_compute filter_gene_min_cpm_cell)
        D = (sp.diags(1e6 / lib) @ X).tocsc()
        cell_cpm = np.asarray(D.mean(0)).ravel()
        de = cell_cpm > 5.0
        # production's gate: pseudobulk CPM, threshold 10
        totals = np.asarray(X.sum(0)).ravel()
        b_ctrl = np.log1p(TS * totals / totals.sum())
        expr = np.expm1(b_ctrl)
        pb_cpm = 1e6 * expr / expr.sum()
        prod = (pb_cpm > 10.0) & (b_ctrl > 0)
        print(f"{ctx:>4}{int(de.sum()):>24d}{int(prod.sum()):>23d}"
              f"{int((de & prod).sum()):>8d}{int((de & ~prod).sum()):>21d}"
              f"{int((prod & ~de).sum()):>21d}")
        rows[ctx] = (de, prod, pb_cpm, cell_cpm, totals, lib, D)
        del A, X

    print()
    print("Interpretation: keep_top = 8000 is applied to the PRODUCTION gate only.")
    for ctx in "ABC":
        de, prod, pb_cpm, *_ = rows[ctx]
        n_prod = int(prod.sum())
        n_tail_prod = max(n_prod - 8000, 0)
        n_band = int((de & ~prod).sum())
        print(f"  ctx {ctx}: production gate {n_prod}; effect-ranked tail inside it "
              f"{n_tail_prod}; scorer-visible band below the production gate {n_band}; "
              f"total scorer-visible tail {n_tail_prod + n_band} of {int(de.sum())} "
              f"DE-universe genes ({(n_tail_prod + n_band) / int(de.sum()):.1%})")


if __name__ == "__main__":
    main()
