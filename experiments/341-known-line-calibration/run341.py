"""Known-line held-target pilot, deliberately NOT unseen-context validation."""

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
import core341 as C

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE/"work"
PRIOR = ROOT/"experiments/334-text-conditioned-transfer"
OLD_WORK = PRIOR/"recovery-work"
sys.path.insert(0, str(PRIOR))
import run334 as B
import core334 as P


def inputs():
    own = [HERE/n for n in ("PREREGISTRATION.md", "core341.py", "run341.py",
                            "test341.py", "selection.json")]
    source = [ROOT/b/n for b in B.BANKS for n in
              ("delta.npy", "perts.txt", "genes.txt", "b_ctrl.npy", "meta.json")]
    prior = [PRIOR/n for n in ("core334.py", "run334.py", "selection.json",
                               "inputs.json", "recovery-inputs.json")]
    prior += [OLD_WORK/line/n for line in B.LINES for n in
              ("evaluation.npz", "prediction-0.npz", "prediction-1.npz")]
    return own+source+prior+[B.HIST, ROOT/"uv.lock"]


def freeze():
    if (HERE/"inputs.json").exists() or (HERE/"selection.json").exists():
        raise FileExistsError("already frozen")
    selection = json.loads((PRIOR/"selection.json").read_text())
    selection.update(inherited_selection_sha256=B.sha(PRIOR/"selection.json"),
                     information_regime="same-line training targets available",
                     unseen_context_validation=False, cross_study_validation=False)
    B.write(HERE/"selection.json", selection)
    B.write(HERE/"inputs.json", {"sha256": {
        str(path.relative_to(ROOT)): B.sha(path) for path in inputs()}})


def verify():
    manifest = json.loads((HERE/"inputs.json").read_text())
    assert set(manifest["sha256"]) == {str(path.relative_to(ROOT)) for path in inputs()}
    for path, expected in manifest["sha256"].items():
        assert B.sha(ROOT/path) == expected, path
    own = [str(path.relative_to(ROOT)) for path in inputs() if path.parent == HERE]
    own.append(str((HERE/"inputs.json").relative_to(ROOT)))
    for path in own:
        subprocess.run(["git", "ls-files", "--error-unmatch", path], cwd=ROOT,
                       stdout=subprocess.DEVNULL, check=True)
    assert not subprocess.check_output(["git", "status", "--porcelain", "--", *own], cwd=ROOT)


def report(execution, selection):
    records, values = [], {arm: {} for arm in C.ARMS}
    for line in B.LINES:
        with np.load(WORK/line/"predictions.npz") as archive:
            truth = archive["truth"]
            assert list(archive["targets"]) == selection["targets"]
            np.testing.assert_allclose(P.endpoints(truth, truth)["retrieval"], 1.)
            np.testing.assert_allclose(P.endpoints(np.zeros_like(truth), truth)["retrieval"], .5)
            for arm in C.ARMS:
                metrics = {**P.endpoints(archive[arm], truth),
                           **C.sign_endpoints(archive[arm], truth)}
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
        if comparator == "same_line_full":
            continue
        comparisons[comparator] = {}
        for metric, primary in values["same_line_full"].items():
            delta = primary-values[comparator][metric]
            if metric in ("squared_error", "absolute_error"):
                delta = -delta
            comparisons[comparator][metric] = C.interval(delta)
    positive = lambda a, m, floor=0.: min(
        comparisons[a][m]["target_ci95"][0], comparisons[a][m]["crossed_ci95"][0]) > floor
    gates = {f"weighted_sign_beats_{a}": positive(a, "weighted_sign")
             for a in ("crossline_diag", "crossline_full", "same_line_diag", "labels_shuffled_full")}
    gates["weighted_sign_point"] = comparisons["crossline_diag"]["weighted_sign"]["mean"] >= .010
    for a in ("crossline_diag", "crossline_full"):
        gates[f"top64_beats_{a}"] = positive(a, "top64_weighted_sign")
        gates[f"retrieval_noninferior_{a}"] = positive(a, "retrieval", -.005)
    jurkat_index = B.LINES.index("nadig_jurkat")
    jurkat = {a: C.interval((values["same_line_full"]["weighted_sign"]-
                            values[a]["weighted_sign"])[jurkat_index:jurkat_index+1])
              for a in ("crossline_diag", "crossline_full")}
    gates["jurkat_sign"] = all(v["target_ci95"][0] > 0 for v in jurkat.values())
    primary = comparisons["crossline_diag"]["weighted_sign"]
    depth = {key: float(np.array(primary["line_means"])[selection[key]].mean())
             for key in ("lower_cell_lines", "higher_cell_lines")}
    gates["line_and_cell_number_robustness"] = (
        primary["positive_lines"] >= 5 and all(v > 0 for v in depth.values()))
    gates["integrity"] = all(r["held_target_poison_pass"] and
        r["source_mean_max_difference"] == 0 and r["crossline_diag_max_difference"] == 0
        and r["subset_linear_max_difference"] < 1e-12 for r in execution["fits"])
    result = {
        "status": "EARNS_SEPARATE_CROSS_STUDY_COUNT_VALIDATION" if all(gates.values()) else
                  "DOES_NOT_EARN_FURTHER_DEVELOPMENT",
        "gates": gates, "comparisons": comparisons, "jurkat_weighted_sign": jurkat,
        "means": {a: {m: float(v.mean()) for m, v in metrics.items()}
                  for a, metrics in values.items()},
        "cell_number_sensitivity": depth,
        "scope": "known-line held-target development, not unseen-context or cross-study evidence",
        "same_line_training_labels_used": True,
        "all_six_competition_members_unpriced": True,
        "inputs_sha256": B.sha(HERE/"inputs.json"), "execution_sha256": B.sha(WORK/"execution.json"),
        "per_target_sha256": B.sha(path), "reserve_read": False, "submission": False,
    }
    B.write(WORK/"results.json", result)
    print(json.dumps({k: result[k] for k in ("status", "gates", "means")}, indent=2))


def run():
    verify()
    if (WORK/"execution.json").exists() or list(WORK.glob("*/models-*.npz")):
        raise FileExistsError("prior run exists")
    selection = json.loads((HERE/"selection.json").read_text())
    started, clock = datetime.now(timezone.utc).isoformat(), time.monotonic()
    targets, genes = selection["targets"], selection["genes"]
    historical = B.known()
    loaded = [B.load(bank, targets, genes) for bank in B.BANKS]
    raw, basal = np.array([a for a, _ in loaded]), np.array([b for _, b in loaded])
    direction, rms = P.rms_direction(raw)
    del loaded, raw
    order = P.order(targets, "334-split|")
    blocks = [np.array(order[::2]), np.array(order[1::2])]
    receipts = []
    for outer, line in enumerate(B.LINES):
        directory = WORK/line
        directory.mkdir(parents=True, exist_ok=True)
        donors = [i for i in range(6) if i != outer]
        baseline = direction[donors].mean(0)
        with np.load(OLD_WORK/line/"evaluation.npz") as archive:
            np.testing.assert_array_equal(direction[outer], archive["truth"])
            np.testing.assert_array_equal(baseline, archive["baseline"])
        outputs = {arm: np.full_like(baseline, np.nan) for arm in C.ARMS}
        covered = np.zeros(len(targets), bool)
        for block in (0, 1):
            if time.monotonic()-clock > 2700:
                raise TimeoutError("45-minute budget")
            train, test = blocks[1-block], blocks[block]
            assert not set(train)&set(test) and not covered[test].any()
            covered[test] = True
            inner = P.source_blocks(direction, basal, rms, outer, train)
            # Preserve 334's accumulation order for exact baseline replay.
            cross_gain = P.diagonal(inner)
            x = np.concatenate([features[..., 0] for _, features, _ in inner])
            cross_y = np.concatenate([response for _, _, response in inner])
            same_y = np.tile(direction[outer, train], (len(donors), 1))
            permutation = P.derangement(len(train), C.SEED+block+20)
            shuffled_y = np.tile(direction[outer, train[permutation]], (len(donors), 1))
            models = {
                "crossline_full": C.fit(x, cross_y, gain=cross_gain),
                "same_line_full": C.fit(x, same_y),
                "labels_shuffled_full": C.fit(x, shuffled_y),
            }
            query = baseline[test]
            predictions = {"source_mean": query, "crossline_diag": query*cross_gain,
                           "same_line_diag": query*models["same_line_full"]["gain"]}
            predictions.update({a: C.predict(model, query) for a, model in models.items()})
            with np.load(OLD_WORK/line/f"prediction-{block}.npz") as archive:
                np.testing.assert_array_equal(test, archive["test_indices"])
                differences = {
                    "source_mean": float(np.abs(predictions["source_mean"]-archive["source_mean"]).max()),
                    "crossline_diag": float(np.abs(predictions["crossline_diag"]-archive["diagonal"]).max())}
                assert all(v == 0 for v in differences.values())
            poisoned, poisoned_rms = direction.copy(), rms.copy()
            poisoned[outer, test], poisoned_rms[outer, test] = 1e6, 1e6
            pin = P.source_blocks(poisoned, basal, poisoned_rms, outer, train)
            px = np.concatenate([features[..., 0] for _, features, _ in pin])
            py = np.tile(poisoned[outer, train], (len(donors), 1))
            np.testing.assert_array_equal(x, px)
            np.testing.assert_array_equal(same_y, py)
            np.testing.assert_array_equal(cross_gain, P.diagonal(pin))
            pmodel = C.fit(px, py)
            for name in pmodel:
                np.testing.assert_array_equal(pmodel[name], models["same_line_full"][name])
            np.testing.assert_array_equal(predictions["same_line_full"],
                C.predict(pmodel, poisoned[donors][:, test].mean(0)))
            subset = np.mean([C.predict(models["same_line_full"],
                direction[[d for d in donors if d != omit]][:, test].mean(0))
                for omit in donors], axis=0)
            error = float(np.abs(subset-predictions["same_line_full"]).max())
            assert error < 1e-12
            path = directory/f"models-{block}.npz"
            np.savez_compressed(path, train_indices=train, test_indices=test,
                donor_indices=donors, label_permutation=permutation,
                **{f"{arm}__{name}": value for arm, model in models.items()
                   for name, value in model.items()})
            for arm, prediction in predictions.items():
                assert np.isfinite(prediction).all()
                outputs[arm][test] = prediction
            receipts.append({
                "outer": line, "block": block, "held_target_poison_pass": True,
                "source_mean_max_difference": differences["source_mean"],
                "crossline_diag_max_difference": differences["crossline_diag"],
                "subset_linear_max_difference": error, "models_sha256": B.sha(path)})
            print(f"{line} block {block}: matched fits and held-target isolation complete", flush=True)
        assert covered.all()
        path = directory/"predictions.npz"
        np.savez_compressed(path, **outputs, truth=direction[outer], targets=targets, genes=genes)
        for receipt in receipts[-2:]:
            receipt["predictions_sha256"] = B.sha(path)
    verify()
    execution = {
        "start_utc": started, "finish_utc": datetime.now(timezone.utc).isoformat(),
        "seconds": time.monotonic()-clock, "fits": receipts, "historical_checks": historical,
        "max_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "inputs_sha256": B.sha(HERE/"inputs.json"),
        "same_line_training_labels_used": True, "unseen_context_validation": False,
    }
    B.write(WORK/"execution.json", execution)
    report(execution, selection)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("freeze", "verify", "run"))
    with threadpool_limits(limits=6):
        globals()[parser.parse_args().action]()
