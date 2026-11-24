"""Terminal paired analysis. Never a leaderboard score or submission gate."""

import csv
import json
import numpy as np
import core334 as C
from run334 import HERE, WORK, LINES, sha, write


def report():
    if (WORK/"results.json").exists():
        raise FileExistsError("terminal report exists")
    selection = json.loads((HERE/"selection.json").read_text())
    execution = json.loads((WORK/"execution.json").read_text())
    values = {a:{} for a in C.ARMS}
    table = []
    for line in LINES:
        directory = WORK/line
        with np.load(directory/"evaluation.npz") as z:
            truth = z["truth"]
            assert list(z["targets"]) == selection["targets"]
        predicted = {a:np.full_like(truth, np.nan) for a in C.ARMS}
        coverage = np.zeros(len(truth), bool)
        for block in (0, 1):
            receipt = json.loads((directory/f"fit-{block}.json").read_text())
            assert sha(directory/"evaluation.npz") == receipt["evaluation_sha256"]
            assert sha(directory/f"prediction-{block}.npz") == receipt["prediction_sha256"]
            assert sha(directory/f"models-{block}.pickle") == receipt["models_sha256"]
            assert sha(directory/f"training-{block}.npz") == receipt["training_sha256"]
            with np.load(directory/f"prediction-{block}.npz") as z:
                ix = z["test_indices"]
                assert not coverage[ix].any()
                coverage[ix] = True
                for a in C.ARMS:
                    predicted[a][ix] = z[a]
        assert coverage.all()
        np.testing.assert_allclose(C.endpoints(truth, truth)["retrieval"], 1.)
        np.testing.assert_allclose(C.endpoints(np.zeros_like(truth), truth)["retrieval"], .5)
        for arm, pred in predicted.items():
            endpoints = C.endpoints(pred, truth)
            for metric, v in endpoints.items():
                values[arm].setdefault(metric, []).append(v)
            for i, t in enumerate(selection["targets"]):
                table.append({"line": line, "target": t, "arm": arm,
                              **{k:float(v[i]) for k,v in endpoints.items()}})
    csv_path = WORK/"per-target.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)
    values = {a:{m:np.array(v) for m,v in ms.items()} for a,ms in values.items()}
    comparisons = {}
    for comparator in C.ARMS:
        if comparator == "text":
            continue
        comparisons[comparator] = {}
        for metric, v in values["text"].items():
            delta = v-values[comparator][metric]
            if metric in ("squared_error", "absolute_error"):
                delta = -delta
            comparisons[comparator][metric] = C.interval(delta)
    required = ("source_mean", "source_only", "string", "text_shuffled", "diagonal")
    gates = {f"retrieval_beats_{a}": min(
        comparisons[a]["retrieval"]["target_ci95"][0],
        comparisons[a]["retrieval"]["crossed_ci95"][0]) > 0 for a in required}
    gates["retrieval_point"] = comparisons["source_mean"]["retrieval"]["mean"] >= .005
    gates["cosine_noninferiority"] = all(comparisons[a]["cosine"]["mean"] >= -.005
                                        for a in ("source_mean", "diagonal"))
    gates["squared_error"] = comparisons["source_mean"]["squared_error"]["mean"] >= 0
    contrast = comparisons["source_only"]["retrieval"]
    gates["line_robustness"] = (contrast["positive_lines"] >= 5 and
                               contrast["minimum_leave_one_line_out"] > 0)
    line_delta = np.array(contrast["line_means"])
    depth = {key:float(line_delta[selection[key]].mean())
             for key in ("lower_cell_lines", "higher_cell_lines")}
    gates["cell_number_sensitivity"] = all(v > 0 for v in depth.values())
    gates["integrity"] = all(r["poison_isolation_pass"] for r in execution["fits"])
    result = {"status": "EARNS_SEPARATE_COUNT_TEST" if all(gates.values()) else
                       "DOES_NOT_EARN_FURTHER_DEVELOPMENT",
              "gates": gates, "n_lines": len(LINES), "n_targets": C.N_TARGETS,
              "n_genes": C.N_GENES, "comparisons": comparisons,
              "means": {a:{m:float(v.mean()) for m,v in ms.items()} for a,ms in values.items()},
              "cell_number_sensitivity": depth,
              "all_six_competition_members_unpriced": True,
              "scope": "unit-RMS raw source directions; reused development data",
              "per_target_sha256": sha(csv_path), "inputs_sha256": sha(HERE/"inputs.json"),
              "execution_sha256": sha(WORK/"execution.json"),
              "reserve_read": False, "submission": False}
    write(WORK/"results.json", result)
    print(json.dumps({"status": result["status"], "gates": gates, "means": result["means"]}, indent=2))


if __name__ == "__main__":
    report()
