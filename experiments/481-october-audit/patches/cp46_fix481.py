"""481 F2 — proposed replacement for `experiments/417-bulk-sharpen/ntclevel417.py::cp46` (NOT applied; analysis copy).

The re-level's claim is Cp46 x sum(w): the mean over the context's NTC constructs of cell_eval2's delete-1 jackknife on
each construct's cells, standing in for the reference's per-perturbation noise C_p on `cells_per_pert` (400) cells. The
jackknife scales as 1/n (481 `cpsize481.py`, context D: x2.065 at 200 cells, x0.5002 at 800), so a construct of n cells
is put on the reference's footing by the factor n / cells_per_pert. At n == cells_per_pert the factor is exactly 1.0,
so every A/B/C (46 x 400) value -- and every re-level string ever printed -- is unchanged bit for bit.

Guards: refuses a construct with fewer than `min_cells` cells (default 100: the jackknife of a tiny group is noise), and
returns the per-construct sizes so the caller can print them. `cells_per_pert` comes from the bundle manifest.

The diff against ntclevel417.py is in `ntclevel417_cp46.diff` next to this file.
"""
from __future__ import annotations

import numpy as np


def cp46_rows(matrix, ids: np.ndarray, cells_per_pert: int, jackknife, *, min_cells: int = 100, ts: float = 50000.0):
    """Mean/sd of the size-normalised jackknife over NTC constructs, plus the sizes. `jackknife(csr_block) -> float`."""
    vals, sizes = [], []
    for name in dict.fromkeys(ids):
        rows = np.flatnonzero(ids == name)
        n = int(rows.size)
        if n < min_cells:
            raise SystemExit(f"NTC construct {name!r} has {n} cells (< {min_cells}): its jackknife cannot proxy C_p")
        vals.append(jackknife(matrix[rows]) * (n / float(cells_per_pert)))
        sizes.append(n)
    v = np.asarray(vals, dtype=np.float64)
    return float(v.mean()), float(v.std()), int(v.size), sizes
