"""Step 6: prove `pds` and `mse` cannot move, using cell_eval2's own functions.

cell_eval2 0.16.0 builds the expression comparator from RAW per-group column SUMS
(`prep._grouped_sums` -> `prep.bulk_lognorm_means`), not from per-cell normalized values,
so if the group sums are equal as integers the bulk is bit-identical and every metric
computed off it -- `pds_cosine` (a rank of cosine distances between bulk rows) and
`expr_mse_unbiased`'s leading term -- is bit-identical too.

The one cell-level term in that family is `moments.jackknife_correction`, which is
SUBTRACTED from the raw mse. More spread raises it, so it can only lower `expr_mse`; and
the official mse is clamped, already sitting at 0 from far above its own baseline
(2.175 normalized against 0.989), so it cannot move in either direction.
"""
import os
import sys

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib093 as L  # noqa: E402

from cell_eval2.prep import (_group_row_index, _grouped_sums,  # noqa: E402
                             bulk_lognorm_means)
from cell_eval2.moments import jackknife_correction  # noqa: E402

TS_BULK = 5e4
gate = L.Gate()
tot = np.load(os.path.join(L.SCRATCH, "totals093.npy"))
sel = [0, 37, 88]
names = [gate.perts[i] for i in sel]
print(f"perturbations {names}")

policies = {"production s=0.25": np.full(gate.G_full, 0.25),
            "global s=1.0": np.full(gate.G_full, 1.0),
            "keep850 |eff| tail 1.0": None}
key_i = {}
blocks = {}
for lab in policies:
    rowsets = []
    for i in sel:
        if lab.startswith("keep850"):
            vec = np.full(gate.G_full, 1.0)
            score = L.rank_key(tot[i], gate)
            vec[gate.keep_idx[np.argpartition(-score, 849)[:850]]] = 0.25
        else:
            vec = policies[lab]
        rng = np.random.default_rng(4242 + i)
        rowsets.append(L.emit_block(tot[i], vec, rng))
    blocks[lab] = rowsets

labels = np.repeat(np.asarray(names), L.N_EMIT)
uniq, order, bounds = _group_row_index(labels)
ref_sums = None
for lab, rowsets in blocks.items():
    X = sp.csr_matrix(np.vstack(rowsets).astype(np.float32))
    sums = _grouped_sums(X, order, bounds, uniq.size)
    bulk = bulk_lognorm_means(sums, TS_BULK)
    codes = np.searchsorted(uniq, labels)
    jk = jackknife_correction(X, codes, uniq.size, TS_BULK)
    if ref_sums is None:
        ref_sums, ref_bulk = sums, bulk
        same_sum = same_bulk = True
    else:
        same_sum = bool(np.array_equal(sums, ref_sums))
        same_bulk = bool(np.array_equal(bulk, ref_bulk))
    print(f"{lab:<24} group sums identical: {str(same_sum):<5} "
          f"bulk_lognorm identical: {str(same_bulk):<5} "
          f"max|dbulk| {np.abs(bulk - ref_bulk).max():.3e}   "
          f"jackknife C_p mean {jk.mean():.6g}")
    del X

# The per-cell library sizes DO move (the scattered part reallocates counts across cells),
# which is the only route by which the DE side's per-cell CPM normalization can shift lfc.
for lab, rowsets in blocks.items():
    libs = np.concatenate([r.sum(axis=1) for r in rowsets])
    print(f"{lab:<24} per-cell library mean {libs.mean():.1f}  sd {libs.std():.2f}  "
          f"cv {libs.std() / libs.mean():.5f}")
