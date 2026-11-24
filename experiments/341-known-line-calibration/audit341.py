"""Independent reconstruction of 341; never imports its fit or endpoint functions."""

import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE/"work"
BANKS = ("data/donors/k562_gwps", "data/donors/orion_hct116",
         "data/donors/orion_hek293t", "data/donors/rpe1",
         "data/donors/nadig_jurkat",
         "experiments/166-donor-acquisition/donors/nadig_hepg2")


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text())


def resample(delta):
    """Alternate counts-matrix reduction with the prespecified random draw order."""
    n_lines, n_targets = delta.shape
    rng = np.random.default_rng(3410909)
    target, crossed = [], []
    for _ in range(100):
        samples = rng.integers(n_targets, size=(100, n_targets))
        lines = rng.integers(n_lines, size=(100, n_lines))
        weights = np.array([np.bincount(s, minlength=n_targets) for s in samples])/n_targets
        byline = weights@delta.T
        target.extend(byline.mean(1))
        crossed.extend(np.take_along_axis(byline, lines, axis=1).mean(1))
    return np.quantile(target, [.025, .975]), np.quantile(crossed, [.025, .975])


def main():
    for name, expected in read(HERE/"inputs.json")["sha256"].items():
        assert sha(ROOT/name) == expected, name
    selection, result = read(HERE/"selection.json"), read(WORK/"results.json")
    execution = read(WORK/"execution.json")
    for name, key in (("execution.json", "execution_sha256"), ("per-target.csv", "per_target_sha256")):
        assert sha(WORK/name) == result[key]
    assert sha(HERE/"inputs.json") == result["inputs_sha256"]
    targets, genes = selection["targets"], selection["genes"]
    raw = []
    for bank in BANKS:
        path = ROOT/bank
        ti = {t: i for i, t in enumerate((path/"perts.txt").read_text().splitlines())}
        gi = {g: i for i, g in enumerate((path/"genes.txt").read_text().splitlines())}
        matrix = np.load(path/"delta.npy", mmap_mode="r")
        indices = [gi[g] for g in genes]
        raw.append(np.stack([np.array(matrix[ti[t], indices], dtype=float) for t in targets]))
    raw = np.array(raw)
    data = raw/np.maximum(np.sqrt(np.einsum("ltg,ltg->lt", raw, raw)/512)[..., None], 1e-12)
    order = sorted(range(1024), key=lambda i: hashlib.sha256(
        ("334-split|"+targets[i]).encode()).hexdigest())
    blocks = [np.array(order[::2]), np.array(order[1::2])]
    records = list(csv.DictReader((WORK/"per-target.csv").open()))
    keyed = {(r["line"], r["target"], r["arm"]): r for r in records}
    assert len(keyed) == len(records) == 36864
    receipt = {(r["outer"], r["block"]): r for r in execution["fits"]}
    max_normal = max_prediction = max_gain = 0.
    max_endpoint = {m: 0. for m in result["means"]["same_line_full"]}
    values = {a: {m: [] for m in max_endpoint} for a in result["means"]}
    model_count = 0
    for outer, bank in enumerate(BANKS):
        line, directory = Path(bank).name, WORK/Path(bank).name
        donors = [j for j in range(6) if j != outer]
        with np.load(directory/"predictions.npz") as archive:
            stored = {a: archive[a] for a in result["means"]}
            truth = archive["truth"]
            np.testing.assert_array_equal(archive["targets"], targets)
            np.testing.assert_array_equal(archive["genes"], genes)
        np.testing.assert_allclose(truth, data[outer], atol=2e-14, rtol=2e-14)
        covered = np.zeros(1024, bool)
        for block in (0, 1):
            record = receipt[(line, block)]
            assert sha(directory/"predictions.npz") == record["predictions_sha256"]
            path = directory/f"models-{block}.npz"
            assert sha(path) == record["models_sha256"]
            with np.load(path) as archive:
                train, test = archive["train_indices"], archive["test_indices"]
                np.testing.assert_array_equal(train, blocks[1-block])
                np.testing.assert_array_equal(test, blocks[block])
                np.testing.assert_array_equal(archive["donor_indices"], donors)
                assert len(train) == len(test) == 512 and not set(train)&set(test)
                assert not covered[test].any()
                covered[test] = True
                means = [data[[j for j in donors if j != r]][:, train].mean(0) for r in donors]
                labels = [data[r, train] for r in donors]
                x, cross_y = np.concatenate(means), np.concatenate(labels)
                query = data[donors][:, test].mean(0)
                np.testing.assert_allclose(query, stored["source_mean"][test], atol=2e-14, rtol=2e-14)
                permutation = archive["label_permutation"]
                cycle = np.random.default_rng(3410909+block+20).permutation(512)
                expected_perm = np.empty(512, dtype=np.int64)
                expected_perm[cycle] = np.roll(cycle, 1)
                np.testing.assert_array_equal(permutation, expected_perm)
                assert np.all(permutation != np.arange(512))
                regimes = {
                    "crossline_full": cross_y,
                    "same_line_full": np.tile(data[outer, train], (5, 1)),
                    "labels_shuffled_full": np.tile(data[outer, train[permutation]], (5, 1)),
                }
                for arm, y in regimes.items():
                    beta, scale = archive[f"{arm}__beta"], archive[f"{arm}__scale"]
                    fitted_gain = archive[f"{arm}__gain"]
                    den = np.einsum("ij,ij->j", x, x)
                    gain = (np.einsum("ij,ij->j", x, y)+np.median(den))/np.maximum(
                        den+np.median(den), 1e-12)
                    error = float(np.max(np.abs(gain-fitted_gain)))
                    assert error < 2e-13, (line, block, arm, error)
                    max_gain = max(max_gain, error)
                    actual_scale = np.maximum(np.linalg.norm(x, axis=0)/np.sqrt(len(x)), 1e-12)
                    np.testing.assert_allclose(actual_scale, scale, atol=2e-13, rtol=2e-13)
                    z = x/scale
                    residual = ((z.T@z+len(x)*np.eye(512))@beta-z.T@(y-x*gain))/len(x)
                    error = float(np.max(np.abs(residual)))
                    assert error < 2e-12, (line, block, arm, error)
                    max_normal = max(max_normal, error)
                    operator = np.diag(fitted_gain)+beta/scale[:, None]
                    error = float(np.max(np.abs(query@operator-stored[arm][test])))
                    assert error < 2e-12, (line, block, arm, error)
                    max_prediction = max(max_prediction, error)
                    if arm != "labels_shuffled_full":
                        diag = "crossline_diag" if arm == "crossline_full" else "same_line_diag"
                        np.testing.assert_allclose(query*gain, stored[diag][test],
                                                   atol=2e-13, rtol=2e-13)
                    model_count += 1
        assert covered.all()
        weight = np.abs(truth)
        # lexsort uses the explicit coordinate as the secondary tie key.
        top = np.lexsort((np.broadcast_to(np.arange(512), weight.shape), -weight), axis=1)[:, :64]
        for arm, prediction in stored.items():
            dot = prediction@truth.T
            denominator = np.outer(np.linalg.norm(prediction, axis=1), np.linalg.norm(truth, axis=1))
            similarity = np.divide(dot, denominator, out=np.zeros_like(dot), where=denominator > 1e-24)
            own = np.diag(similarity)
            correct = prediction*truth > 0
            hit = weight*correct
            endpoints = {
                "cosine": own,
                "decoy_cosine": (similarity.sum(1)-own)/1023,
                "retrieval": (np.diag(rankdata(similarity, axis=1, method="average"))-1)/1023,
                "absolute_error": np.mean(np.abs(prediction-truth), axis=1),
                "squared_error": np.mean(np.square(prediction-truth), axis=1),
                "weighted_sign": hit.sum(1)/weight.sum(1),
                "top64_weighted_sign": np.take_along_axis(hit, top, axis=1).sum(1)/
                                       np.take_along_axis(weight, top, axis=1).sum(1),
            }
            for metric, vector in endpoints.items():
                expected = np.array([float(keyed[(line, t, arm)][metric]) for t in targets])
                error = float(np.max(np.abs(expected-vector)))
                assert error < 2e-12, (line, arm, metric, error)
                max_endpoint[metric] = max(max_endpoint[metric], error)
                values[arm][metric].append(vector)
    max_interval, interval_count = 0., 0
    for comparator in result["comparisons"]:
        for metric in max_endpoint:
            delta = np.array(values["same_line_full"][metric])-np.array(values[comparator][metric])
            if metric in ("absolute_error", "squared_error"):
                delta = -delta
            target, crossed = resample(delta)
            expected = result["comparisons"][comparator][metric]
            error = max(float(np.max(np.abs(target-expected["target_ci95"]))),
                        float(np.max(np.abs(crossed-expected["crossed_ci95"]))),
                        abs(float(delta.mean())-expected["mean"]))
            assert error < 2e-12, (comparator, metric, error)
            np.testing.assert_allclose(delta.mean(1), expected["line_means"], atol=2e-12)
            max_interval = max(max_interval, error)
            interval_count += 1
    for comparator, expected in result["jurkat_weighted_sign"].items():
        delta = (np.array(values["same_line_full"]["weighted_sign"])-
                 np.array(values[comparator]["weighted_sign"]))[4:5]
        target, crossed = resample(delta)
        error = max(float(np.max(np.abs(target-expected["target_ci95"]))),
                    float(np.max(np.abs(crossed-expected["crossed_ci95"]))))
        assert error < 2e-12
        max_interval = max(max_interval, error)
        interval_count += 1
    output = {
        "status": "PASS", "full_ridge_models_audited": model_count,
        "ridge_output_optima_audited": model_count*512, "endpoint_rows_rebuilt": len(records),
        "maximum_normal_equation_residual_per_training_row": max_normal,
        "maximum_dense_prediction_difference": max_prediction,
        "maximum_diagonal_difference": max_gain, "maximum_endpoint_error": max_endpoint,
        "interval_pairs_rebuilt": interval_count, "maximum_interval_error": max_interval,
        "result_sha256": sha(WORK/"results.json"), "inputs_sha256": sha(HERE/"inputs.json"),
        "audit_code_sha256": sha(Path(__file__)),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "same-data reconstruction, not independent biological or cross-study replication",
        "not_recomputed": ["upstream donor extraction", "original held-target poison execution"],
    }
    path = WORK/"independent-audit.json"
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(output, indent=2, sort_keys=True, allow_nan=False)+"\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    with threadpool_limits(limits=6):
        main()
