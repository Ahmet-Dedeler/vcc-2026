"""Preregistered source-only calibration; no production or submission writes."""

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
import core339 as C

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE/"work"
PRIOR = ROOT/"experiments/334-text-conditioned-transfer"
PRIOR_WORK = PRIOR/"recovery-work"
sys.path.insert(0, str(PRIOR))
import run334 as B
import core334 as P


def inputs():
    own = [HERE/n for n in ("PREREGISTRATION.md", "core339.py", "test339.py",
                            "run339.py", "selection.json")]
    inherited = [PRIOR/n for n in ("run334.py", "core334.py", "selection.json",
                                   "inputs.json", "recovery-inputs.json",
                                   "recovery-work/execution.json")]
    inherited += [ROOT/b/n for b in B.BANKS for n in
                  ("delta.npy", "perts.txt", "genes.txt", "b_ctrl.npy", "meta.json")]
    inherited += [PRIOR_WORK/line/n for line in B.LINES for n in
                  ("evaluation.npz", "prediction-0.npz", "prediction-1.npz")]
    return own+inherited+[B.HIST, ROOT/"uv.lock"]


def freeze():
    if (HERE/"inputs.json").exists():
        raise FileExistsError("already frozen")
    selection = json.loads((PRIOR/"selection.json").read_text())
    selection["inherited_selection_sha256"] = B.sha(PRIOR/"selection.json")
    selection["experiment_339_new_response_values_read"] = False
    selection["curve_gene_permutation"] = P.derangement(len(selection["genes"]), C.SEED).tolist()
    if (HERE/"selection.json").exists():
        # A pre-freeze path check failed before any response values were
        # opened. Resume only if that metadata selection is unchanged.
        assert json.loads((HERE/"selection.json").read_text()) == selection
    else:
        B.write(HERE/"selection.json", selection)
    B.write(HERE/"inputs.json", {"sha256": {
        str(p.relative_to(ROOT)): B.sha(p) for p in inputs()}})


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
    rows = []
    values = {arm: {} for arm in C.ARMS}
    for line in B.LINES:
        with np.load(WORK/line/"predictions.npz") as archive:
            truth = archive["truth"]
            assert list(archive["targets"]) == selection["targets"]
            np.testing.assert_allclose(P.endpoints(truth, truth)["retrieval"], 1.)
            np.testing.assert_allclose(P.endpoints(np.zeros_like(truth), truth)["retrieval"], .5)
            for arm in C.ARMS:
                endpoint = P.endpoints(archive[arm], truth)
                for metric, array in endpoint.items():
                    values[arm].setdefault(metric, []).append(array)
                rows.extend({"line": line, "target": target, "arm": arm,
                             **{m: float(v[i]) for m, v in endpoint.items()}}
                            for i, target in enumerate(selection["targets"]))
    path = WORK/"per-target.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    values = {a: {m: np.array(v) for m, v in metrics.items()}
              for a, metrics in values.items()}
    comparisons = {}
    for other in C.ARMS:
        if other == "marginal":
            continue
        comparisons[other] = {}
        for metric, v in values["marginal"].items():
            delta = v-values[other][metric]
            if metric in ("absolute_error", "squared_error"):
                delta = -delta
            comparisons[other][metric] = C.interval(delta)
    gates = {f"retrieval_beats_{other}": min(
        comparisons[other]["retrieval"]["target_ci95"][0],
        comparisons[other]["retrieval"]["crossed_ci95"][0]) > 0
        for other in comparisons}
    primary = comparisons["source_mean"]
    gates["retrieval_point"] = primary["retrieval"]["mean"] >= .005
    gates["cosine_noninferiority"] = all(
        comparisons[a]["cosine"]["mean"] >= -.005 for a in ("source_mean", "diagonal"))
    gates["squared_error"] = primary["squared_error"]["mean"] >= 0
    gates["line_robustness"] = (
        primary["retrieval"]["positive_lines"] >= 5 and
        primary["retrieval"]["minimum_leave_one_line_out"] > 0)
    delta = np.array(primary["retrieval"]["line_means"])
    depth = {key: float(delta[selection[key]].mean())
             for key in ("lower_cell_lines", "higher_cell_lines")}
    gates["cell_number_sensitivity"] = all(d > 0 for d in depth.values())
    gates["integrity"] = all(
        r["poison_isolation_pass"] and r["source_mean_max_difference"] == 0 and
        r["diagonal_max_difference"] == 0 for r in execution["fits"])
    result = {
        "status": "EARNS_SEPARATE_COUNT_TEST" if all(gates.values()) else
                  "DOES_NOT_EARN_FURTHER_DEVELOPMENT",
        "gates": gates, "comparisons": comparisons,
        "means": {a: {m: float(v.mean()) for m, v in ms.items()} for a, ms in values.items()},
        "cell_number_sensitivity": depth,
        "scope": "unit-RMS raw source directions on reused six-line development data",
        "all_six_competition_members_unpriced": True,
        "inputs_sha256": B.sha(HERE/"inputs.json"),
        "per_target_sha256": B.sha(path), "execution_sha256": B.sha(WORK/"execution.json"),
        "reserve_read": False, "submission": False,
    }
    B.write(WORK/"results.json", result)
    print(json.dumps({k: result[k] for k in ("status", "gates", "means")}, indent=2))


def run():
    verify()
    if (WORK/"execution.json").exists() or list(WORK.glob("*/curves-*.npz")):
        raise FileExistsError("prior run exists")
    started, clock = datetime.now(timezone.utc).isoformat(), time.monotonic()
    selection = json.loads((HERE/"selection.json").read_text())
    targets, genes = selection["targets"], selection["genes"]
    checks = B.known()
    loaded = [B.load(bank, targets, genes) for bank in B.BANKS]
    raw = np.array([a for a, _ in loaded])
    basal = np.array([b for _, b in loaded])
    direction, rms = P.rms_direction(raw)
    del loaded, raw
    order = P.order(targets, "334-split|")
    blocks = [np.array(order[::2]), np.array(order[1::2])]
    permutation = np.array(selection["curve_gene_permutation"])
    receipts = []
    for outer, line in enumerate(B.LINES):
        directory = WORK/line
        directory.mkdir(parents=True, exist_ok=True)
        donors = [j for j in range(len(B.LINES)) if j != outer]
        baseline = direction[donors].mean(0)
        with np.load(PRIOR_WORK/line/"evaluation.npz") as archive:
            np.testing.assert_array_equal(direction[outer], archive["truth"])
            np.testing.assert_array_equal(baseline, archive["baseline"])
        predictions = {a: np.full_like(baseline, np.nan) for a in C.ARMS}
        coverage = np.zeros(len(targets), bool)
        for block in (0, 1):
            if time.monotonic()-clock > 2700:
                raise TimeoutError("45-minute total budget")
            train, test = blocks[1-block], blocks[block]
            assert not set(train) & set(test)
            assert not coverage[test].any()
            coverage[test] = True
            training = direction[donors][:, train]
            query = direction[donors][:, test]
            curves = C.fit_curves(training)
            gain = P.diagonal(P.source_blocks(direction, basal, rms, outer, train))
            out = {"source_mean": baseline[test], "diagonal": baseline[test]*gain}
            for arm in ("marginal", "gaussian", "clipping"):
                out[arm] = C.predict(curves, query, arm)
            out["curves_shuffled"] = C.predict(curves, query, gene_permutation=permutation)
            with np.load(PRIOR_WORK/line/f"prediction-{block}.npz") as archive:
                np.testing.assert_array_equal(test, archive["test_indices"])
                differences = {a: float(np.abs(out[a]-archive[a]).max())
                               for a in ("source_mean", "diagonal")}
                assert all(v == 0 for v in differences.values()), differences
            poisoned, poisoned_rms = direction.copy(), rms.copy()
            poisoned[outer], poisoned_rms[outer] = 1e6, 1e6
            pcurves = C.fit_curves(poisoned[donors][:, train])
            for key in curves:
                np.testing.assert_array_equal(curves[key], pcurves[key])
            np.testing.assert_array_equal(out["marginal"],
                C.predict(pcurves, poisoned[donors][:, test]))
            np.testing.assert_array_equal(gain, P.diagonal(
                P.source_blocks(poisoned, basal, poisoned_rms, outer, train)))
            curves_path = directory/f"curves-{block}.npz"
            np.savez_compressed(curves_path, **curves, gain=gain,
                                donor_indices=donors, train_indices=train, test_indices=test)
            for arm, prediction in out.items():
                assert prediction.shape == (len(test), len(genes)) and np.isfinite(prediction).all()
                predictions[arm][test] = prediction
            receipts.append({
                "outer": line, "block": block, "poison_isolation_pass": True,
                "source_mean_max_difference": differences["source_mean"],
                "diagonal_max_difference": differences["diagonal"],
                "curves_sha256": B.sha(curves_path),
            })
            print(f"{line} block {block}: baseline/diagonal exact, calibration complete", flush=True)
        assert coverage.all()
        path = directory/"predictions.npz"
        np.savez_compressed(path, **predictions, truth=direction[outer],
                            targets=targets, genes=genes)
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
