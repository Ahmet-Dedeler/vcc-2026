"""Frozen graph-local response pilot; no candidate or production mutation."""

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
from threadpoolctl import threadpool_limits
import core340 as C

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE/"work"
PRIOR = ROOT/"experiments/334-text-conditioned-transfer"
PRIOR_WORK = PRIOR/"recovery-work"
sys.path.insert(0, str(PRIOR))
import run334 as B
import core334 as P


def inputs():
    own = [HERE/n for n in ("PREREGISTRATION.md", "core340.py", "run340.py",
                            "test340.py", "selection.json")]
    source = [ROOT/b/n for b in B.BANKS for n in
              ("delta.npy", "perts.txt", "genes.txt", "b_ctrl.npy", "meta.json")]
    prior = [PRIOR/n for n in ("core334.py", "run334.py", "selection.json",
                               "inputs.json", "recovery-inputs.json")]
    prior += [PRIOR_WORK/line/n for line in B.LINES for n in
              ("evaluation.npz", "prediction-0.npz", "prediction-1.npz")]
    return own+source+prior+[B.STRING, B.HIST, ROOT/"uv.lock"]


def freeze():
    if (HERE/"inputs.json").exists() or (HERE/"selection.json").exists():
        raise FileExistsError("already frozen")
    selection = json.loads((PRIOR/"selection.json").read_text())
    genes = selection["genes"]
    with np.load(B.STRING) as archive:
        index = {str(g): i for i, g in enumerate(archive["genes"])}
        vectors = np.zeros((len(genes), 64))
        available = np.zeros(len(genes), bool)
        for gene, name in enumerate(genes):
            if name in index and archive["available"][index[name]]:
                vectors[gene] = archive["unit"][index[name]]
                available[gene] = True
    vectors /= np.maximum(np.linalg.norm(vectors, axis=1, keepdims=True), 1e-12)
    graph = C.neighbors_from_similarity(vectors@vectors.T, available=available, positive=True)
    degree = (graph >= 0).sum(1)
    random = C.random_neighbors(degree)
    selection.update({
        "inherited_selection_sha256": B.sha(PRIOR/"selection.json"),
        "string_sha256": B.sha(B.STRING),
        "string_available_genes": int(available.sum()),
        "graph_neighbors": graph.tolist(), "random_neighbors": random.tolist(),
        "graph_degrees": degree.tolist(), "full_degree_genes": int((degree == C.K).sum()),
        "metadata_gate_pass": bool(np.mean(degree == C.K) >= .9),
        "experiment_340_response_values_read": False,
    })
    B.write(HERE/"selection.json", selection)
    B.write(HERE/"inputs.json", {"sha256": {
        str(p.relative_to(ROOT)): B.sha(p) for p in inputs()}})
    print(json.dumps({k: selection[k] for k in
                      ("string_available_genes", "full_degree_genes", "metadata_gate_pass")}))


def verify():
    manifest = json.loads((HERE/"inputs.json").read_text())
    assert set(manifest["sha256"]) == {str(p.relative_to(ROOT)) for p in inputs()}
    for path, expected in manifest["sha256"].items():
        assert B.sha(ROOT/path) == expected, path
    owned = [str(p.relative_to(ROOT)) for p in inputs() if p.parent == HERE]
    owned.append(str((HERE/"inputs.json").relative_to(ROOT)))
    for path in owned:
        subprocess.run(["git", "ls-files", "--error-unmatch", path], cwd=ROOT,
                       stdout=subprocess.DEVNULL, check=True)
    assert not subprocess.check_output(["git", "status", "--porcelain", "--", *owned], cwd=ROOT)


def report(execution, selection):
    records, values = [], {arm: {} for arm in C.ARMS}
    for line in B.LINES:
        with np.load(WORK/line/"predictions.npz") as archive:
            truth = archive["truth"]
            assert list(archive["targets"]) == selection["targets"]
            np.testing.assert_allclose(P.endpoints(truth, truth)["retrieval"], 1.)
            np.testing.assert_allclose(P.endpoints(np.zeros_like(truth), truth)["retrieval"], .5)
            for arm in C.ARMS:
                metrics = P.endpoints(archive[arm], truth)
                for name, vector in metrics.items():
                    values[arm].setdefault(name, []).append(vector)
                records.extend({"line": line, "target": target, "arm": arm,
                                **{m: float(v[i]) for m, v in metrics.items()}}
                               for i, target in enumerate(selection["targets"]))
    path = WORK/"per-target.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    values = {a: {m: np.array(v) for m, v in metrics.items()} for a, metrics in values.items()}
    comparisons = {}
    for comparator in C.ARMS:
        if comparator == "graph":
            continue
        comparisons[comparator] = {}
        for metric, primary in values["graph"].items():
            delta = primary-values[comparator][metric]
            if metric in ("absolute_error", "squared_error"):
                delta = -delta
            comparisons[comparator][metric] = C.interval(delta)
    gates = {f"retrieval_beats_{a}": min(
        comparisons[a]["retrieval"]["target_ci95"][0],
        comparisons[a]["retrieval"]["crossed_ci95"][0]) > 0 for a in comparisons}
    gates["retrieval_point"] = comparisons["source_mean"]["retrieval"]["mean"] >= .005
    gates["cosine_noninferiority"] = all(
        comparisons[a]["cosine"]["mean"] >= -.005 for a in ("source_mean", "diagonal"))
    gates["squared_error"] = all(
        comparisons[a]["squared_error"]["mean"] >= 0 for a in ("source_mean", "diagonal"))
    primary = comparisons["diagonal"]["retrieval"]
    gates["line_robustness"] = primary["positive_lines"] >= 5 and primary["minimum_leave_one_line_out"] > 0
    byline = np.array(primary["line_means"])
    depth = {key: float(byline[selection[key]].mean())
             for key in ("lower_cell_lines", "higher_cell_lines")}
    gates["cell_number_sensitivity"] = all(v > 0 for v in depth.values())
    gates["integrity"] = selection["metadata_gate_pass"] and all(
        r["poison_isolation_pass"] and r["source_mean_max_difference"] == 0 and
        r["diagonal_max_difference"] == 0 and r["subset_linear_max_difference"] < 1e-12
        for r in execution["fits"])
    result = {
        "status": "EARNS_SEPARATE_COUNT_TEST" if all(gates.values()) else
                  "DOES_NOT_EARN_FURTHER_DEVELOPMENT",
        "gates": gates, "comparisons": comparisons,
        "means": {a: {m: float(v.mean()) for m, v in metrics.items()}
                  for a, metrics in values.items()},
        "cell_number_sensitivity": depth,
        "scope": "graph-local ridge on unit-RMS source directions; reused development data",
        "all_six_competition_members_unpriced": True,
        "inputs_sha256": B.sha(HERE/"inputs.json"),
        "execution_sha256": B.sha(WORK/"execution.json"),
        "per_target_sha256": B.sha(path), "reserve_read": False, "submission": False,
    }
    B.write(WORK/"results.json", result)
    print(json.dumps({k: result[k] for k in ("status", "gates", "means")}, indent=2))


def run():
    verify()
    if (WORK/"execution.json").exists() or list(WORK.glob("*/models-*.npz")):
        raise FileExistsError("prior run exists")
    selection = json.loads((HERE/"selection.json").read_text())
    assert selection["metadata_gate_pass"], "metadata coverage failed; no response fitting"
    started, clock = datetime.now(timezone.utc).isoformat(), time.monotonic()
    targets, genes = selection["targets"], selection["genes"]
    checks = B.known()
    loaded = [B.load(bank, targets, genes) for bank in B.BANKS]
    raw, basal = np.array([a for a, _ in loaded]), np.array([b for _, b in loaded])
    direction, rms = P.rms_direction(raw)
    del loaded, raw
    order = P.order(targets, "334-split|")
    blocks = [np.array(order[::2]), np.array(order[1::2])]
    graph = np.array(selection["graph_neighbors"])
    random = np.array(selection["random_neighbors"])
    degrees = np.array(selection["graph_degrees"])
    receipts = []
    for outer, line in enumerate(B.LINES):
        directory = WORK/line
        directory.mkdir(parents=True, exist_ok=True)
        donors = [i for i in range(6) if i != outer]
        baseline = direction[donors].mean(0)
        with np.load(PRIOR_WORK/line/"evaluation.npz") as archive:
            np.testing.assert_array_equal(direction[outer], archive["truth"])
            np.testing.assert_array_equal(baseline, archive["baseline"])
        outputs = {arm: np.full_like(baseline, np.nan) for arm in C.ARMS}
        covered = np.zeros(len(targets), bool)
        for block in (0, 1):
            if time.monotonic()-clock > 3600:
                raise TimeoutError("one-hour budget")
            train, test = blocks[1-block], blocks[block]
            assert not set(train)&set(test) and not covered[test].any()
            covered[test] = True
            inner = P.source_blocks(direction, basal, rms, outer, train)
            gain = P.diagonal(inner)
            x = np.concatenate([features[..., 0] for _, features, _ in inner])
            y = np.concatenate([response for _, _, response in inner])
            permutation = P.derangement(len(train), C.SEED+block+20)
            ys = np.concatenate([response[permutation] for _, _, response in inner])
            correlation = C.correlation_neighbors(x, degrees)
            graphs = {"own_refit": np.full_like(graph, -1), "graph": graph,
                      "random_graph": random, "correlation_graph": correlation,
                      "labels_shuffled": graph}
            models = {arm: C.fit(x, ys if arm == "labels_shuffled" else y, gain, topology)
                      for arm, topology in graphs.items()}
            query = baseline[test]
            predictions = {"source_mean": query, "diagonal": query*gain}
            predictions.update({arm: C.predict(model, query) for arm, model in models.items()})
            with np.load(PRIOR_WORK/line/f"prediction-{block}.npz") as archive:
                np.testing.assert_array_equal(test, archive["test_indices"])
                differences = {arm: float(np.abs(predictions[arm]-archive[arm]).max())
                               for arm in ("source_mean", "diagonal")}
                assert all(v == 0 for v in differences.values())
            poisoned, poisoned_rms = direction.copy(), rms.copy()
            poisoned[outer], poisoned_rms[outer] = 1e6, 1e6
            pin = P.source_blocks(poisoned, basal, poisoned_rms, outer, train)
            px = np.concatenate([features[..., 0] for _, features, _ in pin])
            py = np.concatenate([response for _, _, response in pin])
            np.testing.assert_array_equal(x, px)
            np.testing.assert_array_equal(y, py)
            pgain = P.diagonal(pin)
            np.testing.assert_array_equal(gain, pgain)
            pmodel = C.fit(px, py, pgain, graph)
            for name in pmodel:
                np.testing.assert_array_equal(pmodel[name], models["graph"][name])
            np.testing.assert_array_equal(predictions["graph"],
                C.predict(pmodel, poisoned[donors][:, test].mean(0)))
            subset_predictions = np.mean([
                C.predict(models["graph"], direction[[d for d in donors if d != omit]][:, test].mean(0))
                for omit in donors], axis=0)
            subset_error = float(np.abs(subset_predictions-predictions["graph"]).max())
            assert subset_error < 1e-12
            model_path = directory/f"models-{block}.npz"
            np.savez_compressed(model_path, train_indices=train, test_indices=test,
                                label_permutation=permutation, donor_indices=donors,
                                **{f"{arm}__{name}": value for arm, model in models.items()
                                   for name, value in model.items()})
            for arm, prediction in predictions.items():
                assert np.isfinite(prediction).all()
                outputs[arm][test] = prediction
            receipts.append({
                "outer": line, "block": block, "poison_isolation_pass": True,
                "source_mean_max_difference": differences["source_mean"],
                "diagonal_max_difference": differences["diagonal"],
                "subset_linear_max_difference": subset_error,
                "models_sha256": B.sha(model_path),
            })
            print(f"{line} block {block}: all five models and isolation checks complete", flush=True)
        assert covered.all()
        path = directory/"predictions.npz"
        np.savez_compressed(path, **outputs, truth=direction[outer], targets=targets, genes=genes)
        for receipt in receipts[-2:]:
            receipt["predictions_sha256"] = B.sha(path)
    verify()
    execution = {
        "start_utc": started, "finish_utc": datetime.now(timezone.utc).isoformat(),
        "seconds": time.monotonic()-clock, "fits": receipts, "historical_checks": checks,
        "max_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "inputs_sha256": B.sha(HERE/"inputs.json"),
    }
    B.write(WORK/"execution.json", execution)
    report(execution, selection)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("freeze", "verify", "run"))
    with threadpool_limits(limits=6):
        globals()[parser.parse_args().action]()
