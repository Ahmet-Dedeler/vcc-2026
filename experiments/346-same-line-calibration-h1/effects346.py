#!/usr/bin/env python3
"""346 step 3 — norm-preserving integration of the calibrated directions into 323's cap15.

For evaluation target i and covered gene index set C (5,275 indices on effects.npz['genes']):

    d      = cap15[i, C];  u_d = d/||d||_2
    p      = arm prediction for target i on C;  u_p = p/||p||_2
    m      = (1-beta)*u_d + beta*u_p
    new[C] = m * (||d||_2 / ||m||_2)
    new[j] = cap15[i, j]  for every gene j outside C

so the modified row's L2 norm on C equals cap15's exactly, every coordinate outside C is
bit-identical, and beta=0 returns cap15 bit-identically. Amplitude is not what changes.
Evaluation targets without a full three-donor consensus keep cap15 unchanged.
Writes only into this experiment's work/.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"

# arm name -> (prediction key in predictions.npz, beta)
BUILD = {
    "diag_cross_b1": ("crossline_diag", 1.0),
    "same_ridge_b1": ("same_line_full", 1.0),
    "same_ridge_b05": ("same_line_full", 0.5),
    "shuffled_b1": ("labels_shuffled_full", 1.0),
    # design-stage-only comparators (not emitted)
    "source_mean_b1": ("source_mean", 1.0),
    "crossline_full_b1": ("crossline_full", 1.0),
    "same_diag_b1": ("same_line_diag", 1.0),
    "identity_b0": ("same_line_full", 0.0),
}
EMITTED = ("cap15", "diag_cross_b1", "same_ridge_b1", "same_ridge_b05", "shuffled_b1")


def blend(design_row, prediction, beta):
    d = np.asarray(design_row, dtype=np.float64)
    nd = np.linalg.norm(d)
    if nd <= 0:
        raise ValueError("zero design row on the covered axis")
    if beta == 0.0:
        return d.copy()
    p = np.asarray(prediction, dtype=np.float64)
    npn = np.linalg.norm(p)
    if npn <= 0:
        raise ValueError("zero prediction on the covered axis")
    m = (1.0 - beta) * (d / nd) + beta * (p / npn)
    nm = np.linalg.norm(m)
    if nm <= 1e-9:
        raise ValueError("degenerate blend direction")
    return m * (nd / nm)


def main() -> None:
    with np.load(BASE / "work/effects.npz", allow_pickle=False) as z:
        genes = [str(g) for g in z["genes"]]
        targets = [str(t) for t in z["targets"]]
        cap15 = np.asarray(z["cap15"], dtype=np.float64)
    with np.load(WORK / "predictions.npz", allow_pickle=False) as z:
        covered = [str(g) for g in z["genes"]]
        test = [str(t) for t in z["test_targets"]]
        predictions = {k: np.asarray(z[k], dtype=np.float64)
                       for k in z.files if k not in
                       ("genes", "test_targets", "training_targets",
                        "pan_test_targets", "cross_gain")}
    gi = {g: i for i, g in enumerate(genes)}
    if len(set(covered)) != len(covered):
        raise ValueError("duplicate covered gene")
    cols = np.array([gi[g] for g in covered], dtype=np.int64)
    ti = {t: i for i, t in enumerate(targets)}
    rows = np.array([ti[t] for t in test], dtype=np.int64)
    outside = np.ones(len(genes), dtype=bool)
    outside[cols] = False

    arms, checks = {"cap15": cap15}, {}
    for name, (key, beta) in BUILD.items():
        new = cap15.copy()
        for j, i in enumerate(rows):
            new[i, cols] = blend(cap15[i, cols], predictions[key][j], beta)
        norm_error = float(np.abs(np.linalg.norm(new[:, cols], axis=1)
                                  - np.linalg.norm(cap15[:, cols], axis=1)).max())
        outside_identical = bool(np.array_equal(new[:, outside], cap15[:, outside]))
        unchanged_rows = int(np.sum(np.all(new == cap15, axis=1)))
        cos = np.sum(new[rows][:, cols] * cap15[rows][:, cols], axis=1) / (
            np.linalg.norm(new[rows][:, cols], axis=1)
            * np.linalg.norm(cap15[rows][:, cols], axis=1))
        checks[name] = {
            "beta": beta, "prediction": key,
            "covered_norm_max_abs_error": norm_error,
            "outside_covered_axis_bit_identical": outside_identical,
            "rows_bit_identical_to_cap15": unchanged_rows,
            "median_cosine_with_cap15_on_covered_axis": float(np.median(cos)),
            "finite": bool(np.isfinite(new).all()),
        }
        if not outside_identical or norm_error > 1e-10 or not checks[name]["finite"]:
            raise AssertionError(f"{name}: integration invariants failed")
        arms[name] = new
    if not np.array_equal(arms["identity_b0"], cap15):
        raise AssertionError("beta=0 must return cap15 bit-identically")
    checks["identity_b0"]["beta_zero_bit_identical_to_cap15"] = True

    np.savez_compressed(WORK / "effects346.npz",
                        genes=np.array(genes), targets=np.array(targets),
                        covered_columns=cols, modified_rows=rows,
                        **{k: v for k, v in arms.items() if k != "identity_b0"})
    summary = {
        "n_covered_genes": len(covered), "n_modified_targets": len(rows),
        "n_evaluation_targets": len(targets),
        "unmodified_targets": sorted(set(targets) - set(test)),
        "covered_axis_norm_share_of_cap15_row": {
            "median": float(np.median(np.linalg.norm(cap15[:, cols], axis=1)
                                      / np.linalg.norm(cap15, axis=1))),
            "min": float(np.min(np.linalg.norm(cap15[:, cols], axis=1)
                                / np.linalg.norm(cap15, axis=1))),
            "max": float(np.max(np.linalg.norm(cap15[:, cols], axis=1)
                                / np.linalg.norm(cap15, axis=1)))},
        "integration_checks": checks,
        "emitted_arms": list(EMITTED),
        "rule": "new[C] = ((1-beta)*u_d + beta*u_p) * ||d||/||(1-beta)*u_d + beta*u_p||; "
                "off-C coordinates untouched",
    }
    (WORK / "effects.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
