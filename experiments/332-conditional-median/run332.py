"""Freeze and serialize the source-stage experiment; no production writes."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
import core332 as C

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
OLD = ROOT / "experiments/268-broad16-context-residuals/work"
RUNTIME = ROOT / "experiments/331-tabicl-runtime"
CHECKPOINT = RUNTIME / "work/tabicl-regressor-v2-20260212.ckpt"
CHECKPOINT_SHA = "0db9cb538f114e79026bf08f45f41ad8dd7ad2de2aaca9a5ca8cd3bd9748ae7a"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temp.replace(path)


def paths():
    own = [HERE / n for n in ("PREREGISTRATION.md", "core332.py", "run332.py",
                              "fit332.py", "report332.py", "test332.py")]
    source = [OLD / "halves" / f"{line}.npz" for line in C.LINES]
    source += [OLD / n for n in ("load268.py", "reserve.json", "inventory.json",
                                 "lolo268_K2.csv")]
    source += [ROOT / "data/challenge/pert_counts.csv",
               ROOT / "data/external/depmap26q1/extracted/metadata/cell_line_metadata.csv",
               ROOT / "experiments/323-current-family-validation/inputs.json",
               ROOT / "uv.lock", RUNTIME / "requirements.lock",
               RUNTIME / "inputs.json", RUNTIME / "work/results.json", CHECKPOINT]
    return own + source


def freeze():
    if (HERE / "inputs.json").exists():
        raise FileExistsError("already frozen")
    assert not set(C.LINES) & C.RESERVED
    write(HERE / "inputs.json", {"sha256": {str(p.relative_to(ROOT)): sha(p) for p in paths()},
                                "lines": C.LINES, "reserved_not_read": sorted(C.RESERVED)})


def verify():
    manifest = json.loads((HERE / "inputs.json").read_text())
    assert set(manifest["sha256"]) == {str(p.relative_to(ROOT)) for p in paths()}
    for name, expected in manifest["sha256"].items():
        if sha(ROOT / name) != expected:
            raise ValueError("input drift: " + name)
    own = [str(p.relative_to(ROOT)) for p in paths() if p.parent == HERE]
    own.append(str((HERE / "inputs.json").relative_to(ROOT)))
    for name in own:
        subprocess.run(["git", "ls-files", "--error-unmatch", name], cwd=ROOT,
                       check=True, stdout=subprocess.DEVNULL)
    if subprocess.check_output(["git", "status", "--porcelain", "--", *own],
                               cwd=ROOT, text=True):
        raise ValueError("commit registered inputs before measurement")
    assert sha(CHECKPOINT) == CHECKPOINT_SHA


def load_old():
    spec = importlib.util.spec_from_file_location("load268", OLD / "load268.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert tuple(mod.WORKING) == C.LINES
    assert not os.environ.get("LINES268")
    return mod


def prepare():
    if (WORK / "prepared.json").exists():
        raise FileExistsError("prepared data already exists")
    mod = load_old()
    h, genes = mod.load_halves(list(C.LINES))
    historical = list(csv.DictReader((OLD / "lolo268_K2.csv").open()))
    targets = sorted({r["target"] for r in historical})
    assert len(targets) == 82
    assert targets == sorted(mod.usable_targets(h, list(C.LINES), 6))
    tblocks = C.target_blocks(targets)
    excluded = set(mod.ALL_TARGETS)
    excluded.update(r["target_gene"] for r in csv.DictReader(
        (ROOT / "data/challenge/pert_counts.csv").open()))
    excluded.update(json.loads((ROOT / "experiments/323-current-family-validation/inputs.json")
                               .read_text())["h1_targets"])
    ctrl = []
    for line in C.LINES:
        ix = [h[line]["targets"].index(t) for t in mod.CONTROLS]
        ctrl.append(np.sum(h[line]["full"][ix], axis=0, dtype=np.float64))
    ctrl = np.array(ctrl)
    measured = np.isfinite(ctrl).all(axis=0) & ~np.isin(genes, list(excluded))
    historical_checks, fold_receipts = [], []
    for outer, line in enumerate(C.LINES):
        training = C.source_indices(len(C.LINES), outer)
        old_r, _, _ = mod.build_responses(h, genes, list(C.LINES), targets, keys=("full",),
                                         min_cpm=5., filter_lines=[C.LINES[i] for i in training])
        base_old = old_r["full"][:, training].mean(axis=1)
        truth_old = old_r["full"][:, outer]
        # Exact historical arithmetic, including float32 cosine accumulation.
        norm = np.maximum(np.linalg.norm(base_old, axis=1) *
                          np.linalg.norm(truth_old, axis=1), 1e-12)
        old_cos = (base_old * truth_old).sum(axis=1) / norm
        old_ret = C.retrieval(base_old, truth_old)
        saved = {r["target"]: r for r in historical if r["held"] == line}
        differences = {"cosine": float(np.max(np.abs(old_cos - [
            float(saved[t]["base_cos_full"]) for t in targets]))),
                       "retrieval": float(np.max(np.abs(old_ret - [
                           float(saved[t]["base_retrieval"]) for t in targets])))}
        if max(differences.values()) > 1e-6:
            raise ArithmeticError(f"historical known-answer failed: {line} {differences}")
        historical_checks.append({"line": line, **differences})
        del old_r, base_old, truth_old

        # Only control values select the new response axis.
        cpm = ctrl[training] / np.nansum(ctrl[training], axis=1, keepdims=True) * 1e6
        mask = measured & (np.mean(cpm, axis=0) >= 5.)
        full_genes = genes[mask].astype(str)
        if len(full_genes) < C.N_GENES:
            raise ValueError("too few eligible genes")
        selected = np.array(C.hash_order(full_genes, "332-gene|")[:C.N_GENES])
        selected_genes = full_genes[selected]
        assert not set(selected_genes) & excluded
        responses = {}
        controls = []
        for key in ("full", "A", "B"):
            rows = []
            for name in C.LINES:
                ix = [h[name]["targets"].index(t) for t in targets]
                ci = [h[name]["targets"].index(t) for t in mod.CONTROLS]
                ctl = mod.logprof(np.sum(h[name][key][ci], axis=0), mask)
                r = mod.logprof(h[name][key][ix], mask) - ctl
                if not np.isfinite(r).all():
                    raise ArithmeticError("nonfinite response; do not drop outcome rows")
                rows.append(r)
                if key == "full":
                    controls.append(ctl[selected])
            responses[key] = np.asarray(rows)
        rf = responses["full"]
        rms = np.sqrt(np.mean(rf**2, axis=2))
        selected_r = rf[..., selected]
        blocks, query, wrong, scale = C.fold_features(selected_r, np.array(controls), rms, outer)
        # Held-response poison audit before any fit.
        poisoned, poisoned_rms = selected_r.copy(), rms.copy()
        poisoned[outer] = 1e6
        poisoned_rms[outer] = 1e6
        bp, qp, wp, sp = C.fold_features(poisoned, np.array(controls), poisoned_rms, outer)
        for a, b in zip(blocks, bp):
            for x, y in zip(a[1:], b[1:]):
                np.testing.assert_array_equal(x, y)
        for a, b in ((query, qp), (wrong, wp), (scale, sp)):
            np.testing.assert_array_equal(a, b)

        fold_dir = WORK / line
        fold_dir.mkdir(parents=True, exist_ok=True)
        baseline = selected_r[training].mean(axis=0)
        eval_file = fold_dir / "evaluation.npz"
        np.savez_compressed(eval_file, targets=np.array(targets), genes=selected_genes,
                            truth=selected_r[outer], baseline=baseline, scale=scale,
                            truth_A=responses["A"][outer][:, selected],
                            truth_B=responses["B"][outer][:, selected])
        # Index the line first to preserve T,G rather than advanced-index G,T.
        with np.load(eval_file) as z:
            assert z["truth_A"].shape == baseline.shape, "half-axis order"
        split_receipts = []
        for block in (0, 1):
            test_targets, train_targets = tblocks[block], tblocks[1-block]
            x, y, yshuf, indices = C.sample_training(
                blocks, train_targets, C.SEED + 100*outer + block)
            assert outer not in indices[:, 0]
            assert not set(indices[:, 1]) & set(test_targets)
            assert len(set(map(tuple, indices))) == C.N_TRAIN
            path = fold_dir / f"inputs-{block}.npz"
            np.savez_compressed(path, x_train=x, y_train=y, y_shuffled=yshuf,
                                x_test=query[test_targets].reshape(-1, len(C.FEATURES)),
                                x_wrong=wrong[test_targets].reshape(-1, len(C.FEATURES)),
                                test_indices=test_targets, sample_indices=indices)
            split_receipts.append({"block": block, "inputs_sha256": sha(path),
                                   "training_targets": [targets[i] for i in train_targets],
                                   "test_targets": [targets[i] for i in test_targets],
                                   "n_train": len(x), "n_test_pairs": len(test_targets)*C.N_GENES})
        receipt = {"line": line, "full_genes": len(full_genes), "selected_genes": selected_genes.tolist(),
                   "evaluation_sha256": sha(eval_file), "poison_isolation_pass": True,
                   "wrong_recipient": C.LINES[(outer+1) % len(C.LINES)], "splits": split_receipts}
        write(fold_dir / "prepared.json", receipt)
        fold_receipts.append(receipt)
        print("prepared", line, "genes", len(full_genes), "KA", differences, flush=True)
    write(WORK / "prepared.json", {"folds": fold_receipts, "historical_checks": historical_checks,
                                   "inputs_sha256": sha(HERE / "inputs.json"),
                                   "targets": targets, "excluded_symbols": sorted(excluded),
                                   "features": C.FEATURES, "reserve_read": False})


def run():
    verify()
    if WORK.exists():
        raise FileExistsError("work directory exists; inspect without overwriting")
    WORK.mkdir()
    started, tick = datetime.now(timezone.utc).isoformat(), time.monotonic()
    receipt = {"started_utc": started, "inputs_sha256": sha(HERE / "inputs.json"),
               "status": "RUNNING", "completed_splits": [], "submissions": False}
    write(WORK / "execution.json", receipt)
    try:
        prepare()
        env = dict(os.environ, HF_HUB_OFFLINE="1", HF_HUB_DISABLE_TELEMETRY="1",
                   DO_NOT_TRACK="1", WANDB_MODE="disabled", OMP_NUM_THREADS="4",
                   OPENBLAS_NUM_THREADS="4", VECLIB_MAXIMUM_THREADS="4")
        for line in C.LINES:
            for block in (0, 1):
                if time.monotonic()-tick > 7200:
                    raise TimeoutError("two-hour biological bound")
                subprocess.run([str(RUNTIME / ".venv/bin/python"), str(HERE / "fit332.py"),
                                line, str(block)], cwd=ROOT, env=env, check=True, timeout=300)
                receipt["completed_splits"].append([line, block])
                write(WORK / "execution.json", receipt)
        if time.monotonic()-tick > 7200:
            raise TimeoutError("two-hour biological bound")
        subprocess.run([sys.executable, str(HERE / "report332.py")], cwd=ROOT, check=True)
        verify()
        receipt["status"] = "COMPLETE"
    except BaseException as exc:
        receipt["status"] = "EXECUTION_STOP"
        receipt["reason"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        receipt["seconds"] = time.monotonic()-tick
        receipt["finished_utc"] = datetime.now(timezone.utc).isoformat()
        receipt["driver_max_rss_bytes_macos"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        write(WORK / "execution.json", receipt)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "run", "verify"))
    globals()[parser.parse_args().mode]()
