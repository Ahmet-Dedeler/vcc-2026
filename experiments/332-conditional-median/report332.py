"""Join locked predictions with held responses; all source-stage endpoints."""

import csv
import json
import numpy as np
import core332 as C
from run332 import WORK, HERE, sha, write


def report():
    if (WORK / "results.json").exists():
        raise FileExistsError("terminal analysis exists")
    prepared = json.loads((WORK / "prepared.json").read_text())
    targets = prepared["targets"]
    values = {a: {k: [] for k in ("absolute_error", "squared_error", "cosine",
                                  "retrieval", "orthogonal_cosine")} for a in C.ARMS}
    receipts, diagnostics, table = [], [], []
    for line in C.LINES:
        directory = WORK / line
        prep = json.loads((directory / "prepared.json").read_text())
        evaluation = directory / "evaluation.npz"
        assert sha(evaluation) == prep["evaluation_sha256"]
        with np.load(evaluation) as z:
            truth, baseline, scale = z["truth"], z["baseline"], z["scale"]
            half_a, half_b = z["truth_A"], z["truth_B"]
            assert list(z["targets"]) == targets
        outputs = {a: np.full(truth.shape, np.nan) for a in C.ARMS[2:]}
        covered = np.zeros(len(targets), bool)
        for block in (0, 1):
            receipt = json.loads((directory / f"fit-{block}.json").read_text())
            prediction = directory / f"prediction-{block}.npz"
            assert sha(prediction) == receipt["prediction_sha256"]
            assert receipt["input_sha256"] == prep["splits"][block]["inputs_sha256"]
            receipts.append(receipt)
            with np.load(prediction) as z:
                ix = z["test_indices"]
                assert not covered[ix].any()
                for arm in outputs:
                    outputs[arm][ix] = z[arm].reshape(len(ix), C.N_GENES) * scale[ix, None]
                covered[ix] = True
        assert covered.all()
        outputs["source_mean"], outputs["zero"] = baseline, np.zeros_like(baseline)
        oracle = C.retrieval(truth, truth)
        if not np.allclose(oracle, 1, atol=1e-12):
            raise ArithmeticError("self retrieval instrument failed")
        np.testing.assert_allclose(C.retrieval(np.zeros_like(truth), truth), .5, atol=1e-12)
        diagnostics.append({"line": line, "oracle_retrieval": float(oracle.mean()),
                            "half_cosine": float(C.cosine(half_a, half_b).mean())})
        for arm in C.ARMS:
            pred = outputs[arm]
            assert np.isfinite(pred).all()
            metrics = C.endpoints(pred, truth, baseline, scale)
            for k, v in metrics.items():
                values[arm][k].append(v)
            for i, target in enumerate(targets):
                table.append({"line": line, "target": target, "arm": arm,
                              **{k: float(v[i]) for k, v in metrics.items()}})
    for arm in values:
        for k in values[arm]:
            values[arm][k] = np.array(values[arm][k])
    csv_path = WORK / "per-target.csv"
    with csv_path.open("w", newline="") as stream:
        w = csv.DictWriter(stream, fieldnames=list(table[0]))
        w.writeheader(); w.writerows(table)
    primary = values["tabicl_median"]
    comparisons = {}
    for comparator in C.ARMS:
        if comparator == "tabicl_median":
            continue
        comparisons[comparator] = {}
        for metric in primary:
            delta = primary[metric] - values[comparator][metric]
            if metric in ("absolute_error", "squared_error"):
                delta = -delta
            comparisons[comparator][metric] = C.paired_interval(delta)
    required = ("source_mean", "zero", "no_recipient", "wrong_recipient",
                "shuffled_labels", "absolute_boost")
    gates = {f"absolute_beats_{a}": (
        comparisons[a]["absolute_error"]["target_ci95"][0] > 0
        and comparisons[a]["absolute_error"]["crossed_ci95"][0] > 0) for a in required}
    gates["median_beats_mean"] = comparisons["tabicl_mean"]["absolute_error"]["target_ci95"][0] > 0
    base = comparisons["source_mean"]
    gates["retrieval"] = (base["retrieval"]["mean"] >= .005
                          and base["retrieval"]["target_ci95"][0] >= -.005
                          and base["retrieval"]["crossed_ci95"][0] >= -.005)
    gates["cosine_noninferiority_point"] = base["cosine"]["mean"] >= -.005
    gates["squared_error_noninferiority_point"] = base["squared_error"]["mean"] >= 0
    nc = comparisons["no_recipient"]["absolute_error"]
    gates["line_robustness"] = nc["positive_lines"] >= 8 and nc["minimum_leave_one_line_out"] > 0
    gates["integrity"] = all(f["poison_isolation_pass"] for f in prepared["folds"])
    result = {"status": "EARNS_SEPARATE_EFFECT_STAGE" if all(gates.values()) else
                        "DOES_NOT_EARN_FURTHER_DEVELOPMENT",
              "gates": gates, "comparisons_primary_minus_control_or_loss_improvement": comparisons,
              "means": {a: {k: float(v.mean()) for k, v in m.items()} for a, m in values.items()},
              "n_lines": len(C.LINES), "n_targets": len(targets), "n_genes_per_fold": C.N_GENES,
              "diagnostics": diagnostics, "per_target_sha256": sha(csv_path),
              "inputs_sha256": sha(HERE / "inputs.json"), "scoring": "source-stage proxies only",
              "all_six_competition_members_unpriced": True, "reserve_read": False,
              "submission": False, "fit_seconds": float(sum(r["seconds"] for r in receipts)),
              "largest_fit_max_rss_bytes_macos": max(r["max_rss_bytes_macos"] for r in receipts)}
    write(WORK / "results.json", result)
    print(json.dumps({"status": result["status"], "gates": gates, "means": result["means"]}, indent=2))


if __name__ == "__main__":
    report()
