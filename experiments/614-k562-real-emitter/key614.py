#!/usr/bin/env python3
"""614 stage 4: the PROXY design key for K562 (TRUSTGATE.md).
For each K562 target it averages the ORION donors' delta rows (log1p space at TS 50k, the production key's own units)
on the K562 gene axis, times ONE scalar s fitted ON H1 ONLY:
  s = median |cap15| / median |3-donor mean delta| over H1's (target, gene) cells where both are nonzero;
  cap15 = 323's production key (effects.npz); the 3 donors are k562_gwps, orion_hct116, orion_hek293t (323's DONORS).
No K562 truth is read, and there is no K562 donor and no PerturbMatch K562 replacement.
Output: work/key614.npz (eff, genes, targets), in the K562 reference's target order (work/targets614.txt sorted)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DON = ROOT / "data/donors"


def donor_rows(name, targets, genes):
    P = [l.strip() for l in open(DON / name / "perts.txt")]
    G = [l.strip() for l in open(DON / name / "genes.txt")]
    pi = {p: i for i, p in enumerate(P)}
    gc = {}
    for i, g in enumerate(G):
        if g and g not in gc:
            gc[g] = i
    M = np.load(DON / name / "delta.npy", mmap_mode="r")
    cols = np.array([gc.get(g, -1) for g in genes])
    ok = cols >= 0
    out = np.full((len(targets), len(genes)), np.nan)
    for k, t in enumerate(targets):
        if t in pi:
            out[k, ok] = np.asarray(M[pi[t]], dtype=np.float64)[cols[ok]]
    return out


def main():
    z = np.load(ROOT / "experiments/323-current-family-validation/work/effects.npz", allow_pickle=False)
    cap, g1, t1 = np.asarray(z["cap15"], float), z["genes"].astype(str), z["targets"].astype(str)
    c3 = np.nanmean(np.stack([donor_rows(d, t1, g1) for d in ("k562_gwps", "orion_hct116", "orion_hek293t")]), 0)
    both = np.isfinite(c3) & (c3 != 0) & (cap != 0)
    s = float(np.median(np.abs(cap[both])) / np.median(np.abs(c3[both])))
    rho = float(np.corrcoef(cap[both], c3[both])[0, 1])
    sgn = float(np.mean(np.sign(cap[both]) == np.sign(c3[both])))
    print(f"H1 scale fit: s = {s:.4f} over {int(both.sum())} cells; Pearson(cap15, 3-donor mean) {rho:.3f}; sign agreement {sgn:.3f}")
    import h5py
    with h5py.File(HERE / "work/k562-gwps-c250.h5", "r") as h:
        genes = h["var/_index"][:].astype(str)
    targets = np.array(sorted(l.strip() for l in open(HERE / "work/targets614.txt") if l.strip()))
    rows = np.stack([donor_rows(d, targets, genes) for d in ("orion_hct116", "orion_hek293t")])
    has = np.isfinite(rows).any(2)
    cons = np.nanmean(rows, 0)
    eff = s * np.nan_to_num(cons, nan=0.0)
    np.savez(HERE / "work/key614.npz", eff=eff, genes=genes, targets=targets)
    info = {"scale_s": s, "h1_cells": int(both.sum()), "h1_pearson": rho, "h1_sign_agreement": sgn,
            "targets": int(len(targets)), "targets_with_any_orion_row": int(has.any(0).sum()),
            "targets_with_both_orion_rows": int(has.all(0).sum()),
            "key_nonzero_frac": float((eff != 0).mean()), "key_median_abs_nz": float(np.median(np.abs(eff[eff != 0]))),
            "h1_cap15_median_abs_nz": float(np.median(np.abs(cap[cap != 0])))}
    (HERE / "work/key614.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info))


if __name__ == "__main__":
    main()
