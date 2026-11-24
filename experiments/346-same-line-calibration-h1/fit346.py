#!/usr/bin/env python3
"""346 step 2 — fit the same-line (Pan H1) calibration maps and predict the 85 evaluation targets.

Estimator primitives are imported verbatim from 341/334; nothing is copied or retuned.
Recipient responses for the evaluation targets are never read during fitting; the poison
replay proves it by refitting with those rows replaced by noise and requiring byte identity.
Writes only into this experiment's work/.
"""
from __future__ import annotations

import csv
import json
import os
import resource
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
for _p in (str(ROOT / "experiments/334-text-conditioned-transfer"),
           str(ROOT / "experiments/341-known-line-calibration")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import core334 as P    # noqa: E402
import core341 as C    # noqa: E402

SEED = 3460910
PAN = ROOT / "experiments/115-gse295214-h1/work/compact-donor"
DONORS = ("data/donors/k562_gwps", "data/donors/orion_hct116",
          "data/donors/orion_hek293t")
LINES = ("k562_gwps", "orion_hct116", "orion_hek293t", "gse295214_h1")
ARMS = ("source_mean", "crossline_diag", "crossline_full", "same_line_diag",
        "same_line_full", "labels_shuffled_full")


def names(path: Path) -> list[str]:
    return path.read_text().splitlines()


def bank(path: Path, targets, genes):
    """334's loader shape: exact-name row/column selection with a scalar-index fixture."""
    perts, gene_names = names(path / "perts.txt"), names(path / "genes.txt")
    raw = np.load(path / "delta.npy", mmap_mode="r")
    if raw.shape != (len(perts), len(gene_names)):
        raise ValueError(f"{path}: delta shape {raw.shape} does not match its name lists")
    pi = {p: i for i, p in enumerate(perts)}
    gi = {g: i for i, g in enumerate(gene_names)}
    rows = [pi[t] for t in targets]
    cols = [gi[g] for g in genes]
    out = np.array(raw[np.ix_(rows, cols)], dtype=np.float64)
    basal = np.array(np.load(path / "b_ctrl.npy")[cols], dtype=np.float64)
    if not (np.isfinite(out).all() and np.isfinite(basal).all()):
        raise ValueError(f"{path}: non-finite values in the selected block")
    np.testing.assert_array_equal(
        out[:3, :3], [[raw[pi[t], gi[g]] for g in genes[:3]] for t in targets[:3]])
    return out, basal


def build_axes():
    panel = sorted({r["target_gene"] for r in
                    csv.DictReader((ROOT / "data/challenge/pert_counts.csv").open())})
    canonical = [r["gene_name"] for r in
                 csv.DictReader((ROOT / "data/challenge/gene_names.csv").open())]
    with np.load(BASE / "work/reference.npz", allow_pickle=False, mmap_mode="r") as z:
        evaluation = [str(t) for t in z["targets"]]
    pan_perts = names(PAN / "perts.txt")
    if names(PAN / "genes.txt") != canonical:
        raise ValueError("Pan bank is not on the canonical gene axis")
    mask = np.load(PAN / "measured_mask.npy", mmap_mode="r")
    row0 = np.asarray(mask[0], dtype=bool)
    step = max(1, mask.shape[0] // 8)
    for start in range(0, mask.shape[0], step):     # row-constant mask, checked in slabs
        if not bool((np.asarray(mask[start:start + step], dtype=bool) == row0).all()):
            raise ValueError("Pan measured mask is not row-constant")
    measured = {g for g, k in zip(canonical, row0) if k}

    donor_perts = [set(names(ROOT / d / "perts.txt")) for d in DONORS]
    donor_genes = [set(names(ROOT / d / "genes.txt")) for d in DONORS]
    excluded = set(panel) | set(evaluation)
    training = sorted((set.intersection(*donor_perts) & set(pan_perts)) - excluded)
    covered = sorted((set.intersection(*donor_genes) & set(canonical) & measured)
                     - set(training) - excluded)
    test = [t for t in evaluation if all(t in s for s in donor_perts)]
    uncovered = [t for t in evaluation if t not in set(test)]
    pan_test = [t for t in test if t in set(pan_perts)]
    return dict(panel=panel, canonical=canonical, evaluation=evaluation, training=training,
                covered=covered, test=test, uncovered=uncovered, pan_test=pan_test,
                pan_perts=pan_perts)


def assemble(axes, poison: bool):
    """(4, n_targets, n_genes) directions. Pan rows for the test targets are never fitted."""
    order = axes["training"] + axes["test"]
    genes = axes["covered"]
    raws, basals = [], []
    for d in DONORS:
        a, b = bank(ROOT / d, order, genes)
        raws.append(a)
        basals.append(b)
    pan_rows, pan_basal = bank(PAN, axes["training"] + axes["pan_test"], genes)
    pan = np.zeros((len(order), len(genes)), dtype=np.float64)
    pan[:len(axes["training"])] = pan_rows[:len(axes["training"])]
    have = {t: i for i, t in enumerate(axes["pan_test"], start=len(axes["training"]))}
    for j, t in enumerate(axes["test"], start=len(axes["training"])):
        if t in have:
            pan[j] = pan_rows[have[t]]
    if poison:
        rng = np.random.default_rng(SEED + 77)
        pan[len(axes["training"]):] = 1e6 * rng.standard_normal(
            (len(axes["test"]), len(genes)))
    raws.append(pan)
    basals.append(pan_basal)
    raw = np.asarray(raws)
    del raws
    direction, rms = P.rms_direction(raw)
    if poison:
        rms[3, len(axes["training"]):] = 1e6
    return direction, np.asarray(basals), rms, raw


def fit_all(axes, direction, basal, rms):
    n_train = len(axes["training"])
    train = np.arange(n_train)
    test = np.arange(n_train, direction.shape[1])
    outer = 3
    donors = [0, 1, 2]
    baseline = direction[donors].mean(0)
    inner = P.source_blocks(direction, basal, rms, outer, train)
    cross_gain = P.diagonal(inner)
    x = np.concatenate([f[..., 0] for _, f, _ in inner])
    cross_y = np.concatenate([r for _, _, r in inner])
    same_y = np.tile(direction[outer, train], (len(donors), 1))
    permutation = P.derangement(n_train, SEED)
    shuffled_y = np.tile(direction[outer, train[permutation]], (len(donors), 1))
    models = {"crossline_full": C.fit(x, cross_y, gain=cross_gain),
              "same_line_full": C.fit(x, same_y),
              "labels_shuffled_full": C.fit(x, shuffled_y)}
    del cross_y, same_y, shuffled_y
    query = baseline[test]
    predictions = {"source_mean": query.copy(),
                   "crossline_diag": query * cross_gain,
                   "same_line_diag": query * models["same_line_full"]["gain"]}
    predictions.update({a: C.predict(m, query) for a, m in models.items()})
    # 341's subset-linearity identity: the mean of the leave-one-out subset predictions
    # equals the prediction on the full donor mean, for every linear map.
    subset = np.mean([C.predict(models["same_line_full"],
                                direction[[d for d in donors if d != omit]][:, test].mean(0))
                      for omit in donors], axis=0)
    subset_error = float(np.abs(subset - predictions["same_line_full"]).max())
    residual = float(np.abs(
        (x / models["same_line_full"]["scale"]).T
        @ ((x / models["same_line_full"]["scale"]) @ models["same_line_full"]["beta"])
        + len(x) * models["same_line_full"]["beta"]
        - (x / models["same_line_full"]["scale"]).T
        @ (np.tile(direction[outer, train], (len(donors), 1)) - x * models["same_line_full"]["gain"])
    ).max())
    return models, predictions, dict(cross_gain=cross_gain, permutation=permutation,
                                     subset_linear_max_difference=subset_error,
                                     normal_equation_max_residual=residual,
                                     n_rows=int(len(x)), n_features=int(x.shape[1]))


def main() -> None:
    started, clock = datetime.now(timezone.utc).isoformat(), time.monotonic()
    WORK.mkdir(parents=True, exist_ok=True)
    axes = build_axes()
    print(json.dumps({k: len(v) for k, v in axes.items() if isinstance(v, list)}), flush=True)

    direction, basal, rms, raw = assemble(axes, poison=False)
    models, predictions, receipt = fit_all(axes, direction, basal, rms)
    print(f"fit done ({time.monotonic()-clock:.1f}s) "
          f"subset={receipt['subset_linear_max_difference']:.3e} "
          f"residual={receipt['normal_equation_max_residual']:.3e}", flush=True)

    # Held-target poison replay: refit with the recipient's evaluation-target rows replaced
    # by 1e6-scale noise. Every fitted parameter and every prediction must be byte-identical.
    pdirection, pbasal, prms, _ = assemble(axes, poison=True)
    np.testing.assert_array_equal(direction[:3], pdirection[:3])
    np.testing.assert_array_equal(rms[:3], prms[:3])
    np.testing.assert_array_equal(basal, pbasal)
    pmodels, ppredictions, preceipt = fit_all(axes, pdirection, pbasal, prms)
    for arm in models:
        for key in models[arm]:
            np.testing.assert_array_equal(models[arm][key], pmodels[arm][key])
    for arm in predictions:
        np.testing.assert_array_equal(predictions[arm], ppredictions[arm])
    np.testing.assert_array_equal(receipt["cross_gain"], preceipt["cross_gain"])
    del pdirection, pbasal, prms, pmodels, ppredictions

    # Secondary, same-study diagnostic: how the arms do against Pan's OWN held-out responses
    # for the 30 evaluation targets that have a direct Pan row. Never used in fitting.
    n_train = len(axes["training"])
    pan_index = {t: j for j, t in enumerate(axes["test"])}
    held = np.array([pan_index[t] for t in axes["pan_test"]])
    truth = direction[3, n_train:][held]
    same_study = {}
    if len(held) >= 3:
        for arm, value in predictions.items():
            e = P.endpoints(value[held], truth)
            s = C.sign_endpoints(value[held], truth)
            same_study[arm] = {"own_cosine": float(e["cosine"].mean()),
                               "retrieval": float(e["retrieval"].mean()),
                               "weighted_sign": float(s["weighted_sign"].mean()),
                               "top64_weighted_sign": float(s["top64_weighted_sign"].mean())}

    for name, values in (("fitted-training-targets.txt", axes["training"]),
                         ("fitted-covered-genes.txt", axes["covered"]),
                         ("fitted-test-targets.txt", axes["test"])):
        (WORK / name).write_text("\n".join(values) + "\n")
    np.savez_compressed(WORK / "predictions.npz",
                        genes=np.array(axes["covered"]), test_targets=np.array(axes["test"]),
                        training_targets=np.array(axes["training"]),
                        pan_test_targets=np.array(axes["pan_test"]),
                        cross_gain=receipt["cross_gain"],
                        **{a: v for a, v in predictions.items()})
    for arm, model in models.items():
        np.savez_compressed(WORK / f"model-{arm}.npz", **model)
    out = {
        "started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
        "seconds": time.monotonic() - clock,
        "lines": list(LINES), "outer": "gse295214_h1", "seed": SEED,
        "n_training_targets": len(axes["training"]), "n_covered_genes": len(axes["covered"]),
        "n_test_targets_with_full_donor_consensus": len(axes["test"]),
        "evaluation_targets_without_full_donor_consensus": axes["uncovered"],
        "n_evaluation_targets_with_direct_pan_rows_excluded_from_fitting": len(axes["pan_test"]),
        "ridge_rows": receipt["n_rows"], "ridge_features": receipt["n_features"],
        "ridge_lambda_rule": "core341.fit: lambda = number of rows, on RMS-scaled features",
        "held_target_poison_pass": True,
        "subset_linear_max_difference": receipt["subset_linear_max_difference"],
        "normal_equation_max_residual": receipt["normal_equation_max_residual"],
        "same_study_held_target_diagnostic": same_study,
        "same_study_diagnostic_scope": "Pan H1's own held-out responses for the 30 evaluation "
                                       "targets with a direct Pan row; same study, held target; "
                                       "NOT the VCC-2025 H1 reference and not a score",
        "max_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (WORK / "fit.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "same_study_held_target_diagnostic"},
                     indent=2, sort_keys=True))
    print(json.dumps(same_study, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
