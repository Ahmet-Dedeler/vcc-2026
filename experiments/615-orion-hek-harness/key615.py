#!/usr/bin/env python3
"""615 stage 2: the PROXY design key for the Orion raw HEK293T harness (coordinator 18:2xZ: proxy first, as 614).
Both Orion lines are excluded. Of the K562 + Nadig donors, only k562_gwps covers panel targets (269 / 300); nadig_jurkat,
nadig_hepg2 and k562_essential cover 0 / 300, so the proxy is the k562_gwps delta row (log1p space at TS 50k, the
production key's own units) on the HEK truth's gene axis, times ONE scalar s fitted ON H1 ONLY:
  s = median |cap15| / median |k562_gwps delta| over H1's (target, gene) cells where both are nonzero
  (cap15 = 323's production key, effects.npz). No HEK293T data is read.
Targets k562_gwps does not cover get a zero row (no key), as uncovered targets do in production.
Output: work/key615.npz (eff, genes, targets) in the reference's target order (targets615.txt sorted), key615.json."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/614-k562-real-emitter"))
from key614 import donor_rows  # noqa: E402  (read-only: the same row mapper 614 used)

DONORS = ("k562_gwps",)


def main():
    z = np.load(ROOT / "experiments/323-current-family-validation/work/effects.npz", allow_pickle=False)
    cap, g1, t1 = np.asarray(z["cap15"], float), z["genes"].astype(str), z["targets"].astype(str)
    c1 = np.nanmean(np.stack([donor_rows(d, t1, g1) for d in DONORS]), 0)
    both = np.isfinite(c1) & (c1 != 0) & (cap != 0)
    s = float(np.median(np.abs(cap[both])) / np.median(np.abs(c1[both])))
    rho = float(np.corrcoef(cap[both], c1[both])[0, 1])
    sgn = float(np.mean(np.sign(cap[both]) == np.sign(c1[both])))
    print(f"H1 scale fit: s = {s:.4f} over {int(both.sum())} cells; Pearson(cap15, k562_gwps) {rho:.3f}; "
          f"sign agreement {sgn:.3f}")
    with h5py.File(HERE / "work/hek-truth.h5", "r") as h:
        genes = h["var/_index"][:].astype(str)
    targets = np.array(sorted(l.strip() for l in open(HERE / "work/targets615.txt") if l.strip()))
    rows = np.stack([donor_rows(d, targets, genes) for d in DONORS])
    has = np.isfinite(rows).any(2)
    cons = np.nanmean(rows, 0) if rows.shape[0] > 1 else rows[0]
    eff = s * np.nan_to_num(cons, nan=0.0)
    np.savez(HERE / "work/key615.npz", eff=eff, genes=genes, targets=targets)
    info = {"donors": list(DONORS), "scale_s": s, "h1_cells": int(both.sum()), "h1_pearson": rho,
            "h1_sign_agreement": sgn, "targets": int(len(targets)), "targets_with_key": int(has.any(0).sum()),
            "key_nonzero_frac": float((eff != 0).mean()),
            "key_median_abs_nz": float(np.median(np.abs(eff[eff != 0]))) if (eff != 0).any() else None,
            "h1_cap15_median_abs_nz": float(np.median(np.abs(cap[cap != 0]))),
            "excluded": "orion_hct116, orion_hek293t (same platform as the truth); nadig_* / k562_essential: 0 panel cover"}
    (HERE / "work/key615.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info))


if __name__ == "__main__":
    main()
