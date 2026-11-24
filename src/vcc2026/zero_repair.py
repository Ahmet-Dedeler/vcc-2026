"""Zero-gene repair (experiment 428; copied verbatim from experiments/428-zero-gene-audit/repair428.py). A gated gene (control mean CPM > 5) that a perturbation's 400 cells never express has a
delivered arithmetic-mean lfc of log2(1e-9 / ctrl), about -35. If the reference calls that gene, the nmae numerator gains
~30 for it. The block total is 0 because the carrier design's knockdown is so deep that pass 1 emits nothing, and
ratio-mode ballast copies the carriers' composition, so it gives 0 too.

Repair: for every gated non-target gene whose block total is 0, give the CARRIERS
n = max(1, round(ctrl_cpm * 2**floor_log2 * sum(carrier libraries) / 1e6)) single counts. Each count goes into a
distinct random carrier cell, and each cell pays for it with one count from its own largest sub-gate gene (control
CPM <= 5, never tested). So every library, every tested gene other than the repaired one, and the ballast are unchanged.
The repaired gene's delivered lfc lands near `floor_log2`. It stays far below the controls, so it remains called down,
with the same sign.
"""
from __future__ import annotations

import numpy as np


def zero_repair(block: np.ndarray, n_car: int, gated: np.ndarray, ctrl_cpm_gated: np.ndarray, target_col: int,
                rng: np.random.Generator, floor_log2: float = -5.0) -> tuple[np.ndarray, dict]:
    """block: (cells, genes) int counts, carriers first (rows < n_car). gated: full-axis bool (scorer gate).
    ctrl_cpm_gated: control mean CPM on the gated genes (same order as np.flatnonzero(gated)). target_col: the
    perturbation's own column (-1 if absent), never repaired. Returns (block, diag); a no-op returns the input object."""
    gidx = np.flatnonzero(gated)
    tot = block[:, gidx].sum(0)
    zero = tot == 0
    if target_col >= 0:
        zero &= gidx != target_col
    if not zero.any():
        return block, {"zero_repaired": 0, "zero_counts_moved": 0}
    out = np.array(block, copy=True)
    sub = np.flatnonzero(~gated)
    if target_col >= 0:
        sub = sub[sub != target_col]
    car_libs = out[:n_car].sum(1)
    S = float(car_libs.sum())
    moved = 0
    names = []
    for j in np.flatnonzero(zero):
        g = int(gidx[j])
        n = int(max(1, round(float(ctrl_cpm_gated[j]) * 2.0 ** floor_log2 * S / 1e6)))
        n = min(n, n_car)
        cells = rng.choice(n_car, size=n, replace=False)
        for c in cells:
            row = out[c, sub]
            k = int(np.argmax(row))
            if row[k] <= 0:
                continue                              # no sub-gate count to pay with: skip this cell
            out[c, sub[k]] -= 1
            out[c, g] += 1
            moved += 1
        names.append(g)
    if not np.array_equal(out.sum(1), block.sum(1)):
        raise AssertionError("zero repair moved a library")
    if out.min() < 0:
        raise AssertionError("zero repair made a negative count")
    keep_cols = np.setdiff1d(gidx, np.asarray(names, dtype=np.int64))
    if not np.array_equal(out[:, keep_cols], block[:, keep_cols]):
        raise AssertionError("zero repair touched another tested gene")
    return out, {"zero_repaired": len(names), "zero_counts_moved": moved}
