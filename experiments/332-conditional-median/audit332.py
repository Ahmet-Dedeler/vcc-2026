"""Independent read-back audit; does not import the experiment's metric code."""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def audit():
    destination = WORK / "audit.json"
    if destination.exists():
        raise FileExistsError("audit exists; do not overwrite")
    execution = json.loads((WORK / "execution.json").read_text())
    if execution["status"] != "COMPLETE":
        raise ValueError("wait for terminal execution, including final input verification")
    result = json.loads((WORK / "results.json").read_text())
    prepared = json.loads((WORK / "prepared.json").read_text())
    registration = json.loads((HERE / "inputs.json").read_text())
    with (WORK / "per-target.csv").open() as f:
        rows = list(csv.DictReader(f))
    saved = {(r["line"], r["target"], r["arm"]): r for r in rows}
    assert len(saved) == len(rows) == 12*82*9
    assert sha(WORK / "per-target.csv") == result["per_target_sha256"]
    for name, expected in registration["sha256"].items():
        assert sha(ROOT / name) == expected, name
    forbidden = {"JVE127", "MG63", "SLR23", "UMRC3"}
    assert not forbidden & {f["line"] for f in prepared["folds"]}
    assert len(execution["completed_splits"]) == 24
    assert len(set(map(tuple, execution["completed_splits"]))) == 24
    errors = {m: 0. for m in ("absolute_error", "squared_error", "cosine", "retrieval")}
    values = {a: {m: [] for m in errors} for a in result["means"]}
    for outer, fold in enumerate(prepared["folds"]):
        line = fold["line"]
        directory = WORK / line
        evaluation = directory / "evaluation.npz"
        assert sha(evaluation) == fold["evaluation_sha256"]
        with np.load(evaluation) as z:
            truth, baseline, scale = z["truth"], z["baseline"], z["scale"]
            targets = z["targets"].tolist()
        assert truth.shape == (82, 512) and np.isfinite(truth).all()
        predictions = {a: np.zeros_like(truth) for a in result["means"]}
        predictions["source_mean"] = baseline
        covered = set()
        for split in fold["splits"]:
            block = split["block"]
            inputs = directory / f"inputs-{block}.npz"
            pred_file = directory / f"prediction-{block}.npz"
            fit = json.loads((directory / f"fit-{block}.json").read_text())
            assert sha(inputs) == fit["input_sha256"] == split["inputs_sha256"]
            assert sha(pred_file) == fit["prediction_sha256"]
            assert fit["outer_truth_read"] is False
            assert not set(split["training_targets"]) & set(split["test_targets"])
            with np.load(inputs) as z:
                indices = z["sample_indices"]
                assert indices.shape == (2048, 3)
                assert outer not in indices[:, 0]
                assert len(set(map(tuple, indices))) == 2048
                assert not set(indices[:, 1]) & set(z["test_indices"])
                np.testing.assert_array_equal(z["x_test"][:, :10], z["x_wrong"][:, :10])
            with np.load(pred_file) as z:
                ix = z["test_indices"]
                assert not covered & set(ix)
                covered.update(ix)
                for a in predictions:
                    if a in ("source_mean", "zero"):
                        continue
                    predictions[a][ix] = z[a].reshape(len(ix), 512)*scale[ix, None]
        assert covered == set(range(82))
        for arm, prediction in predictions.items():
            standardized_difference = (prediction-truth)/scale[:, None]
            # Independent distance implementation. A zero vector has distance 1.
            distances = cdist(prediction, truth, metric="cosine")
            distances[~np.isfinite(distances)] = 1.
            diagonal = np.diag(distances)
            ranks = (distances < diagonal[:, None]).sum(axis=1)
            ranks = ranks + ((distances == diagonal[:, None]).sum(axis=1)-1)/2
            readback = {
                "absolute_error": np.abs(standardized_difference).sum(axis=1)/512,
                "squared_error": (standardized_difference**2).sum(axis=1)/512,
                "cosine": 1-diagonal,
                "retrieval": 1-ranks/81}
            for metric, v in readback.items():
                expected = np.array([float(saved[(line, t, arm)][metric]) for t in targets])
                difference = float(np.max(np.abs(v-expected)))
                errors[metric] = max(errors[metric], difference)
                if difference > 1e-9:
                    raise ArithmeticError(f"{line}/{arm}/{metric}: {difference}")
                values[arm][metric].append(v)
    for arm, metrics in values.items():
        for metric, v in metrics.items():
            assert abs(float(np.mean(v))-result["means"][arm][metric]) < 1e-10
    rng = np.random.default_rng(332)
    target_weights = rng.multinomial(82, np.full(82, 1/82), size=10000)/82
    line_weights = rng.multinomial(12, np.full(12, 1/12), size=10000)/12
    bootstrap_error = 0.
    for comparator, metrics in result["comparisons_primary_minus_control_or_loss_improvement"].items():
        for metric in errors:
            primary = np.array(values["tabicl_median"][metric])
            control = np.array(values[comparator][metric])
            sign = -1 if metric in ("absolute_error", "squared_error") else 1
            delta = sign*(primary-control)
            assert abs(delta.mean()-metrics[metric]["mean"]) < 1e-10
            np.testing.assert_allclose(delta.mean(axis=1), metrics[metric]["line_means"], atol=1e-10)
            # Same registered draws, different contraction from core332's einsum.
            target_draws = (target_weights @ delta.T).mean(axis=1)
            crossed_draws = ((target_weights @ delta.T)*line_weights).sum(axis=1)
            for draws, key in ((target_draws, "target_ci95"), (crossed_draws, "crossed_ci95")):
                ci = np.quantile(draws, [.025, .975])
                difference = float(np.max(np.abs(ci-metrics[metric][key])))
                bootstrap_error = max(bootstrap_error, difference)
                assert difference < 1e-9, (comparator, metric, key, difference)
    depth = json.loads((WORK / "depth.json").read_text())
    assert depth["source_results_sha256"] == sha(WORK / "results.json")
    payload = {"status": "PASS", "source_results_sha256": sha(WORK / "results.json"),
               "depth_sha256": sha(WORK / "depth.json"), "script_sha256": sha(Path(__file__)),
               "rows_recomputed": len(rows), "prediction_splits_verified": 24,
               "maximum_metric_readback_errors": errors,
               "all_frozen_inputs_unchanged": True,
               "bootstrap_intervals_independently_recomputed": list(errors),
               "maximum_bootstrap_endpoint_error": bootstrap_error,
               "note": "Four endpoints and their paired intervals independently recomputed; "
                       "orthogonal-cosine diagnostic is not independently reproduced. "
                       "Not new biological evidence."}
    temp = destination.with_suffix(".tmp")
    temp.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
    temp.replace(destination)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    audit()
