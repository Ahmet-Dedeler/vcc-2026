"""Independent normal-equation and dense-operator audit of sparse fits."""

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
    rng = np.random.default_rng(3400909)
    target, crossed = [], []
    for _ in range(100):
        samples = rng.integers(1024, size=(100, 1024))
        lines = rng.integers(6, size=(100, 6))
        weights = np.array([np.bincount(s, minlength=1024) for s in samples])/1024
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
    records = list(csv.DictReader((WORK/"per-target.csv").open()))
    keyed = {(r["line"], r["target"], r["arm"]): r for r in records}
    assert len(keyed) == len(records) == 43008
    receipt = {(r["outer"], r["block"]): r for r in execution["fits"]}
    max_normal = max_prediction = max_gain = 0.
    max_endpoint = {m: 0. for m in result["means"]["graph"]}
    values = {a: {m: [] for m in max_endpoint} for a in result["means"]}
    model_count = 0
    for outer, bank in enumerate(BANKS):
        line, directory = Path(bank).name, WORK/Path(bank).name
        donors = [j for j in range(6) if j != outer]
        with np.load(directory/"predictions.npz") as archive:
            stored = {a: archive[a] for a in result["means"]}
            truth = archive["truth"]
        np.testing.assert_allclose(truth, data[outer], atol=2e-14, rtol=2e-14)
        covered = np.zeros(1024, bool)
        for block in (0, 1):
            record = receipt[(line, block)]
            assert sha(directory/"predictions.npz") == record["predictions_sha256"]
            path = directory/f"models-{block}.npz"
            assert sha(path) == record["models_sha256"]
            with np.load(path) as archive:
                train, test = archive["train_indices"], archive["test_indices"]
                assert len(train) == len(test) == 512 and not set(train)&set(test)
                np.testing.assert_array_equal(archive["donor_indices"], donors)
                assert not covered[test].any()
                covered[test] = True
                means = [data[[j for j in donors if j != r]][:, train].mean(0) for r in donors]
                labels = [data[r, train] for r in donors]
                x, y = np.concatenate(means), np.concatenate(labels)
                numerator = sum((p*r).sum(0) for p, r in zip(means, labels))
                denominator = sum((p*p).sum(0) for p in means)
                gain = (numerator+np.median(denominator))/np.maximum(
                    denominator+np.median(denominator), 1e-12)
                query = data[donors][:, test].mean(0)
                np.testing.assert_allclose(query, stored["source_mean"][test], atol=2e-14, rtol=2e-14)
                np.testing.assert_allclose(query*gain, stored["diagonal"][test], atol=2e-13, rtol=2e-13)
                permutation = archive["label_permutation"]
                assert np.array_equal(np.sort(permutation), np.arange(512))
                assert np.all(permutation != np.arange(512))
                ys = np.concatenate([label[permutation] for label in labels])
                centered = x-x.mean(0)
                correlation = np.abs(centered.T@centered)/np.maximum(
                    np.outer(np.linalg.norm(centered, axis=0), np.linalg.norm(centered, axis=0)), 1e-24)
                for arm in ("own_refit", "graph", "random_graph", "correlation_graph", "labels_shuffled"):
                    topology = archive[f"{arm}__neighbors"]
                    beta, scale = archive[f"{arm}__beta"], archive[f"{arm}__scale"]
                    fitted_gain = archive[f"{arm}__gain"]
                    error = float(np.max(np.abs(gain-fitted_gain)))
                    assert error < 2e-13
                    max_gain = max(max_gain, error)
                    if arm in ("graph", "labels_shuffled"):
                        np.testing.assert_array_equal(topology, selection["graph_neighbors"])
                    elif arm == "random_graph":
                        np.testing.assert_array_equal(topology, selection["random_neighbors"])
                    elif arm == "own_refit":
                        np.testing.assert_array_equal(topology, -1)
                    residual = (ys if arm == "labels_shuffled" else y)-x*gain
                    operator = np.diag(fitted_gain)
                    for gene in range(512):
                        valid_neighbors = topology[gene][topology[gene] >= 0]
                        assert gene not in valid_neighbors
                        assert len(set(valid_neighbors)) == len(valid_neighbors)
                        if arm == "correlation_graph" and len(valid_neighbors):
                            remaining = np.setdiff1d(np.arange(512), np.r_[gene, valid_neighbors])
                            assert correlation[gene, valid_neighbors].min() >= (
                                correlation[gene, remaining].max()-1e-12)
                        indices = np.r_[gene, topology[gene]]
                        features = np.zeros((len(x), 9))
                        valid = indices >= 0
                        features[:, valid] = x[:, indices[valid]]
                        actual_scale = np.maximum(np.linalg.norm(features, axis=0)/np.sqrt(len(x)), 1e-12)
                        np.testing.assert_allclose(actual_scale, scale[gene], atol=2e-13, rtol=2e-13)
                        z = features/scale[gene]
                        normal_residual = ((z.T@z+len(x)*np.eye(9))@beta[gene]-
                                           z.T@residual[:, gene])/len(x)
                        error = float(np.max(np.abs(normal_residual)))
                        assert error < 2e-12, (line, block, arm, gene, error)
                        max_normal = max(max_normal, error)
                        for feature, index in enumerate(indices):
                            if index >= 0:
                                operator[index, gene] += beta[gene, feature]/scale[gene, feature]
                    error = float(np.max(np.abs(query@operator-stored[arm][test])))
                    assert error < 2e-12, (line, arm, error)
                    max_prediction = max(max_prediction, error)
                    model_count += 1
        assert covered.all()
        for arm, prediction in stored.items():
            dot = prediction@truth.T
            denominator = np.outer(np.linalg.norm(prediction, axis=1), np.linalg.norm(truth, axis=1))
            similarity = np.divide(dot, denominator, out=np.zeros_like(dot), where=denominator > 1e-24)
            own = np.diag(similarity)
            endpoints = {
                "cosine": own,
                "decoy_cosine": (similarity.sum(1)-own)/1023,
                "retrieval": (np.diag(rankdata(similarity, axis=1, method="average"))-1)/1023,
                "absolute_error": np.mean(np.abs(prediction-truth), axis=1),
                "squared_error": np.mean(np.square(prediction-truth), axis=1),
            }
            for metric, vector in endpoints.items():
                expected = np.array([float(keyed[(line, t, arm)][metric]) for t in targets])
                error = float(np.max(np.abs(expected-vector)))
                assert error < 2e-12, (line, arm, metric, error)
                max_endpoint[metric] = max(max_endpoint[metric], error)
                values[arm][metric].append(expected)
    max_interval = 0.
    for comparator in result["comparisons"]:
        for metric in ("retrieval", "cosine", "squared_error"):
            delta = np.array(values["graph"][metric])-np.array(values[comparator][metric])
            if metric == "squared_error":
                delta = -delta
            target, crossed = resample(delta)
            expected = result["comparisons"][comparator][metric]
            error = max(float(np.max(np.abs(target-expected["target_ci95"]))),
                        float(np.max(np.abs(crossed-expected["crossed_ci95"]))),
                        abs(float(delta.mean())-expected["mean"]))
            assert error < 2e-12
            max_interval = max(max_interval, error)
    output = {
        "status": "PASS", "models_audited": model_count,
        "scalar_ridge_optima_audited": model_count*512, "endpoint_rows_rebuilt": 43008,
        "maximum_normal_equation_residual_per_training_row": max_normal,
        "maximum_dense_prediction_difference": max_prediction,
        "maximum_diagonal_difference": max_gain, "maximum_endpoint_error": max_endpoint,
        "interval_pairs_rebuilt": 18, "maximum_interval_error": max_interval,
        "result_sha256": sha(WORK/"results.json"), "inputs_sha256": sha(HERE/"inputs.json"),
        "audit_code_sha256": sha(Path(__file__)),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "same-data fit/endpoint/interval audit, not independent biological replication",
        "not_recomputed": ["STRING embedding generation", "static metadata graph selection"],
    }
    path = WORK/"independent-audit.json"
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(output, indent=2, sort_keys=True, allow_nan=False)+"\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    with threadpool_limits(limits=6):
        main()
