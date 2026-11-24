#!/usr/bin/env python3
"""346 step 4 — design-stage |lfc|-weighted sign accuracy, computed exactly as 343's signs343.py.

Same gate (reference p_adj < .05, finite non-zero lfc, gene != target, at least 10 gate genes),
same weights (|reference lfc|), same "absent from the effects axis counts as a design miss"
convention, same per-target mean over the 58 eligible targets. Known answer: cap15 must return
343's recorded design_w_all = 0.600185097698462.
Reads 323's reference and this experiment's effect archive; writes only into work/.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
KNOWN_343 = 0.600185097698462
ORDER = ("cap15", "source_mean_b1", "diag_cross_b1", "same_diag_b1", "crossline_full_b1",
         "same_ridge_b05", "same_ridge_b1", "shuffled_b1")


def main() -> None:
    with np.load(BASE / "work/reference.npz", allow_pickle=False, mmap_mode="r") as z:
        genes_all = np.asarray(z["genes"]).astype(str)
        keep = np.asarray(z["keep"], dtype=bool)
        targets = np.asarray(z["targets"]).astype(str)
    genes = genes_all[keep]
    gi = {g: i for i, g in enumerate(genes)}
    with np.load(WORK / "effects346.npz", allow_pickle=False) as z:
        if not np.array_equal(np.asarray(z["targets"]).astype(str), targets):
            raise ValueError("effect target order drift")
        egenes = np.asarray(z["genes"]).astype(str)
        arms = {a: np.asarray(z[a], dtype=np.float64) for a in ORDER}
        covered = set(egenes[np.asarray(z["covered_columns"], dtype=np.int64)])
    eg = {g: i for i, g in enumerate(egenes)}
    col = np.array([eg.get(g, -1) for g in genes])
    miss = col < 0
    on_covered = np.array([g in covered for g in genes])
    mapped = {}
    for a, e in arms.items():
        m = e[:, np.where(miss, 0, col)]
        m[:, miss] = 0.        # sign 0 never matches: counted as a design miss, as in 343
        mapped[a] = m
    # cap15's delivered call set, for 343's called/uncalled split (known-answer check only)
    with np.load(BASE / "work/stats-cap15.npz", allow_pickle=False) as z:
        qhat = np.asarray(z["q"], dtype=np.float64)

    ref = pl.read_parquet(BASE / "work/reference_de.parquet")
    acc = {a: {"design_w_all": [], "design_w_covered": [], "design_u_all": []} for a in ORDER}
    cap_called, cap_uncalled, mass_covered, n_gate = [], [], [], []
    for ti, t in enumerate(targets):
        sub = ref.filter(pl.col("target") == t)
        g = sub["feature"].to_numpy().astype(str)
        q = sub["p_adj"].to_numpy()
        l = sub["log2_fold_change"].to_numpy().astype(float)
        gate = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
        if gate.sum() < 10:
            continue
        idx = np.array([gi[x] for x in g[gate]])
        l = l[gate]
        w = np.abs(l)
        D = w.sum()
        cov = on_covered[idx]
        n_gate.append(int(gate.sum()))
        mass_covered.append(float(w[cov].sum() / D))
        called = qhat[ti, idx] < .05
        for a in ORDER:
            se = np.sign(mapped[a][ti, idx]) == np.sign(l)
            acc[a]["design_w_all"].append(float((w * se).sum() / D))
            acc[a]["design_u_all"].append(float(se.mean()))
            acc[a]["design_w_covered"].append(
                float((w[cov] * se[cov]).sum() / max(w[cov].sum(), 1e-12))
                if cov.any() else np.nan)
            if a == "cap15":
                cap_called.append(float((w[called] * se[called]).sum()
                                        / max(w[called].sum(), 1e-12))
                                  if called.any() else np.nan)
                cap_uncalled.append(float((w[~called] * se[~called]).sum()
                                          / max(w[~called].sum(), 1e-12))
                                    if (~called).any() else np.nan)
    table = {a: {k: float(np.nanmean(v)) for k, v in d.items()} for a, d in acc.items()}
    known = {
        "cap15_design_w_all": table["cap15"]["design_w_all"],
        "expected_343_design_w_all": KNOWN_343,
        "absolute_difference": abs(table["cap15"]["design_w_all"] - KNOWN_343),
        "cap15_design_w_called": float(np.nanmean(cap_called)),
        "expected_343_design_w_called": 0.6492991728970661,
        "cap15_design_w_uncalled": float(np.nanmean(cap_uncalled)),
        "expected_343_design_w_uncalled": 0.5580069530066262,
        "cap15_design_u_all": table["cap15"]["design_u_all"],
        "expected_343_design_u_all": 0.569406545928883,
        "n_targets": len(n_gate), "expected_343_n_targets": 58,
        "mean_gate_genes": float(np.mean(n_gate)), "expected_343_n_gate": 1119.5,
    }
    known["pass"] = bool(
        known["absolute_difference"] < 1e-12
        and abs(known["cap15_design_w_called"] - known["expected_343_design_w_called"]) < 1e-12
        and abs(known["cap15_design_w_uncalled"] - known["expected_343_design_w_uncalled"]) < 1e-12
        and abs(known["cap15_design_u_all"] - known["expected_343_design_u_all"]) < 1e-12
        and known["n_targets"] == 58 and abs(known["mean_gate_genes"] - 1119.5) < 1e-9)
    out = {
        "known_answer": known,
        "arms": table,
        "gate_mass_on_covered_axis_mean": float(np.nanmean(mass_covered)),
        "deltas_vs_cap15": {a: {k: table[a][k] - table["cap15"][k] for k in table[a]}
                            for a in ORDER if a != "cap15"},
        "definition": "343 signs343.py design_w_all: |reference lfc|-weighted fraction of the "
                      "reference's significant gate genes whose DESIGNED effect sign matches, "
                      "averaged over the 58 targets with >=10 gate genes",
    }
    (WORK / "signs.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out, indent=2, sort_keys=True))
    if not known["pass"]:
        raise AssertionError("343 sign known-answer did not reproduce; see work/signs.json")


if __name__ == "__main__":
    main()
