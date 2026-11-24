"""Independent curve, prediction, endpoint, and bootstrap reconstruction."""

import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
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


def independent_map(query, sorted_values, common):
    # Group tied knots explicitly; interpolate without np.interp.
    knots, inverse = np.unique(sorted_values, return_inverse=True)
    if len(knots) < 2:
        return query.copy()
    heights = np.bincount(inverse, weights=common)/np.bincount(inverse)
    points = np.clip(np.r_[0., query], knots[0], knots[-1])
    upper = np.clip(np.searchsorted(knots, points, side="right"), 1, len(knots)-1)
    lower = upper-1
    fraction = (points-knots[lower])/(knots[upper]-knots[lower])
    mapped = heights[lower]+fraction*(heights[upper]-heights[lower])
    return mapped[1:]-mapped[0]


def resample(delta):
    rng = np.random.default_rng(3390909)
    target, crossed = [], []
    for _ in range(100):
        samples = rng.integers(1024, size=(100, 1024))
        lines = rng.integers(6, size=(100, 6))
        # Counts-matrix reduction differs from the report's indexing reduction.
        weights = np.array([np.bincount(s, minlength=1024) for s in samples])/1024
        perline = weights@delta.T
        target.extend(perline.mean(1))
        crossed.extend(np.take_along_axis(perline, lines, axis=1).mean(1))
    return np.quantile(target, [.025, .975]), np.quantile(crossed, [.025, .975])


def main():
    for name, expected in read(HERE/"inputs.json")["sha256"].items():
        assert sha(ROOT/name) == expected, name
    selection, result = read(HERE/"selection.json"), read(WORK/"results.json")
    execution = read(WORK/"execution.json")
    assert sha(HERE/"inputs.json") == result["inputs_sha256"]
    assert sha(WORK/"execution.json") == result["execution_sha256"]
    assert sha(WORK/"per-target.csv") == result["per_target_sha256"]
    records = list(csv.DictReader((WORK/"per-target.csv").open()))
    keyed = {(r["line"], r["target"], r["arm"]): r for r in records}
    assert len(keyed) == len(records) == 36864
    targets, genes = selection["targets"], selection["genes"]
    raw = []
    for bank in BANKS:
        path = ROOT/bank
        ti = {t: i for i, t in enumerate((path/"perts.txt").read_text().splitlines())}
        gi = {g: i for i, g in enumerate((path/"genes.txt").read_text().splitlines())}
        array = np.load(path/"delta.npy", mmap_mode="r")
        raw.append(np.stack([np.array(array[ti[t], [gi[g] for g in genes]], dtype=float)
                             for t in targets]))
    raw = np.array(raw)
    # Alternate sum-of-squares reduction to verify source reconstruction.
    rms = np.sqrt(np.einsum("ltg,ltg->lt", raw, raw)/512)
    data = raw/np.maximum(rms[..., None], 1e-12)
    permutation = np.array(selection["curve_gene_permutation"])
    receipts = {(r["outer"], r["block"]): r for r in execution["fits"]}
    maximum_prediction_error = maximum_curve_error = 0.
    maximum_endpoint_error = {m: 0. for m in result["means"]["marginal"]}
    values = {a: {m: [] for m in maximum_endpoint_error} for a in result["means"]}
    for outer, bank in enumerate(BANKS):
        line, directory = Path(bank).name, WORK/Path(bank).name
        donors = [i for i in range(6) if i != outer]
        with np.load(directory/"predictions.npz") as archive:
            stored = {a: archive[a] for a in result["means"]}
            truth = archive["truth"]
            np.testing.assert_allclose(truth, data[outer], atol=2e-14, rtol=2e-14)
        covered = np.zeros(1024, bool)
        for block in (0, 1):
            receipt = receipts[(line, block)]
            assert sha(directory/"predictions.npz") == receipt["predictions_sha256"]
            path = directory/f"curves-{block}.npz"
            assert sha(path) == receipt["curves_sha256"]
            with np.load(path) as archive:
                train, test = archive["train_indices"], archive["test_indices"]
                assert len(train) == len(test) == 512 and not set(train)&set(test)
                assert not covered[test].any()
                covered[test] = True
                np.testing.assert_array_equal(archive["donor_indices"], donors)
                training, query = data[donors][:, train], data[donors][:, test]
                quantiles = np.sort(training, axis=1)
                barycenter = quantiles.mean(0)
                sd = np.sqrt(np.mean((training-training.mean(1, keepdims=True))**2, axis=1))
                for name, rebuilt in (("sorted", quantiles), ("barycenter", barycenter), ("sd", sd)):
                    error = float(np.max(np.abs(rebuilt-archive[name])))
                    assert error < 2e-13, (name, error)
                    maximum_curve_error = max(maximum_curve_error, error)
                numerator, denominator = np.zeros(512), np.zeros(512)
                for held in donors:
                    others = [j for j in donors if j != held]
                    p = data[others][:, train].mean(0)
                    numerator += (p*data[held, train]).sum(0)
                    denominator += (p*p).sum(0)
                prior = np.median(denominator)
                gain = (numerator+prior)/np.maximum(denominator+prior, 1e-12)
                np.testing.assert_allclose(gain, archive["gain"], atol=2e-13, rtol=2e-13)
            out = {a: np.zeros((512, 512)) for a in result["means"]}
            out["source_mean"] = query.mean(0)
            out["diagonal"] = query.mean(0)*gain
            for donor in range(5):
                for gene in range(512):
                    x, knots = query[donor, :, gene], quantiles[donor, :, gene]
                    out["marginal"][:, gene] += independent_map(x, knots, barycenter[:, gene])/5
                    other = permutation[gene]
                    out["curves_shuffled"][:, gene] += independent_map(
                        x, quantiles[donor, :, other], barycenter[:, other])/5
                    out["gaussian"][:, gene] += x*(
                        sd[:, gene].mean()/sd[donor, gene] if sd[donor, gene] > 1e-12 else 1.)/5
                    clipped = x if knots[0] == knots[-1] else (
                        np.minimum(np.maximum(x, knots[0]), knots[-1])-
                        min(max(0., knots[0]), knots[-1]))
                    out["clipping"][:, gene] += clipped/5
            for arm in out:
                error = float(np.abs(stored[arm][test]-out[arm]).max())
                assert error < 2e-10, (line, block, arm, error)
                maximum_prediction_error = max(maximum_prediction_error, error)
        assert covered.all()
        for arm, prediction in stored.items():
            dot = prediction@truth.T
            norms = np.linalg.norm(prediction, axis=1)[:, None]*np.linalg.norm(truth, axis=1)[None, :]
            similarity = np.divide(dot, norms, out=np.zeros_like(dot), where=norms > 1e-24)
            own = np.diag(similarity)
            endpoints = {
                "absolute_error": np.abs(prediction-truth).mean(1),
                "squared_error": np.square(prediction-truth).mean(1),
                "cosine": own,
                "decoy_cosine": (similarity.sum(1)-own)/1023,
                "retrieval": np.array([
                    (np.count_nonzero(row < own[i])+
                     .5*(np.count_nonzero(row == own[i])-1))/1023
                    for i, row in enumerate(similarity)]),
            }
            for metric, vector in endpoints.items():
                expected = np.array([float(keyed[(line, t, arm)][metric]) for t in targets])
                error = float(np.max(np.abs(expected-vector)))
                assert error < 2e-12, (line, arm, metric, error)
                maximum_endpoint_error[metric] = max(maximum_endpoint_error[metric], error)
                values[arm][metric].append(expected)
    max_interval_error = 0.
    for comparator in result["comparisons"]:
        for metric in ("retrieval", "cosine", "squared_error"):
            delta = np.array(values["marginal"][metric])-np.array(values[comparator][metric])
            if metric == "squared_error":
                delta = -delta
            target, crossed = resample(delta)
            expected = result["comparisons"][comparator][metric]
            error = max(float(np.max(np.abs(target-expected["target_ci95"]))),
                        float(np.max(np.abs(crossed-expected["crossed_ci95"]))),
                        abs(float(delta.mean())-expected["mean"]))
            assert error < 2e-12, error
            max_interval_error = max(max_interval_error, error)
    assert not all(result["gates"].values())
    assert result["status"] == "DOES_NOT_EARN_FURTHER_DEVELOPMENT"
    output = {
        "status": "PASS", "prediction_splits_rebuilt": 12, "endpoint_rows_rebuilt": 36864,
        "maximum_curve_error": maximum_curve_error,
        "maximum_prediction_error": maximum_prediction_error,
        "maximum_endpoint_error": maximum_endpoint_error,
        "independent_interval_pairs": 15, "maximum_interval_error": max_interval_error,
        "result_sha256": sha(WORK/"results.json"), "inputs_sha256": sha(HERE/"inputs.json"),
        "audit_code_sha256": sha(Path(__file__)),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "reconstruction from same development sources, not independent replication",
    }
    path = WORK/"independent-audit.json"
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(output, indent=2, sort_keys=True, allow_nan=False)+"\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    with threadpool_limits(limits=6):
        main()
