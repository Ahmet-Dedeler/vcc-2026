"""Frozen current-family experiment; only own outputs, production imported."""

import argparse
from datetime import datetime, timezone
import gc
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
import pandas as pd

from estimator327 import arms, geometry

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
sys.path.insert(0, str(ROOT / "experiments/326-cap30-current-family"))
import run326 as A

R, E = A.R, A.E
NAMES = ("cap15", "median25", "magnitude25", "power25", "shuffled25")


def input_paths():
    files = [HERE / n for n in ("PREREGISTRATION.md", "estimator327.py", "run327.py",
                                "report327.py", "test327.py")]
    files += [BASE / n for n in ("inputs.json", "evaluation-inputs.json", "anchor-inputs.json",
                                 "run323.py", "evaluate323.py", "report323.py")]
    files += [BASE / "work" / n for n in (
        "stage.npz", "stage.json", "compose.json", "effects.npz",
        "whitener_h1_pruned.npz", "whitener_h1_pruned.json", "reference.npz",
        "reference.json", "reference_de.parquet", "h1-anchors.json",
        "metrics-cap15.json", "stats-cap15.npz", "per-target-cap15.csv")]
    files += [ROOT / "experiments/326-cap30-current-family/run326.py"]
    return files


def freeze():
    path = HERE / "inputs.json"
    if path.exists():
        raise FileExistsError(path)
    R.write_json(path, {
        "inputs_sha256": {str(p.relative_to(ROOT)): R.sha(p) for p in input_paths()},
        "arms": NAMES, "dose": .25, "seed": 20260909,
        "classification": "reused H1 development; no independent confirmation"})


def verify():
    manifest = json.loads((HERE / "inputs.json").read_text())
    expected_paths = {str(p.relative_to(ROOT)) for p in input_paths()}
    if set(manifest["inputs_sha256"]) != expected_paths:
        raise ValueError("input allowlist changed")
    for path, expected in manifest["inputs_sha256"].items():
        if R.sha(ROOT / path) != expected:
            raise ValueError(f"327 input drift: {path}")
    own = [str(p.relative_to(ROOT)) for p in input_paths() if p.parent == HERE]
    own += [str((HERE / "inputs.json").relative_to(ROOT))]
    for path in own:
        subprocess.run(["git", "ls-files", "--error-unmatch", path],
                       cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if subprocess.check_output(["git", "status", "--porcelain", "--", *own], cwd=ROOT, text=True):
        raise ValueError("commit frozen code before measurement")
    E.verify()


def compose():
    verify()
    if (WORK / "effects.npz").exists():
        raise FileExistsError("effects already constructed")
    WORK.mkdir(parents=True, exist_ok=True)
    z = np.load(BASE / "work/stage.npz", allow_pickle=False)
    genes, targets, cols = z["genes"], z["targets"], z["columns"]
    stage = json.loads((BASE / "work/stage.json").read_text())
    comp = json.loads((BASE / "work/compose.json").read_text())
    A.require_same("323 stage", stage["gates_pass"], True)
    A.require_same("323 composition", comp["gates_pass"], True)
    b = R.builder("cap15")
    if b.BETA != 0 or b.ALPHA != .982:
        raise ValueError("unexpected mean-response or dose setting")
    with R.model_env():
        lean = R.load_module("lean327", ROOT / "experiments/170-gamma-mse/stage.py")
    pm = b.EXP061.prepare_perturbmatch(
        ROOT / "data/external/perturbmatch", tuple(genes),
        chunk_rows=128, max_rss_gib=14)
    with lean.lean_norm():
        donors = {n: lean.load_donor_lean(b, ROOT / f"data/donors/{n}", genes)
                  for n in R.DONORS}
        model = b.MultiDonorModel(donors, genes, gene_norm=False, eb=False)
    abc = pd.read_csv(ROOT / "data/challenge/pert_counts.csv").target_gene.tolist()
    lean.trim_model_rows(model, list(targets) + abc)
    del donors
    weights = b.donor_weights(model)
    base = b.build_v9_pre_whiten_effects(model, targets, alpha=b.ALPHA,
                                         beta=b.BETA, weights=weights)
    A.require_same("cap15 pre hash", b.array_sha256(base),
                   stage["arms"]["cap15"]["baseline_pre_sha256"])
    sources = np.zeros((3, len(targets), len(genes)), dtype=np.float64)
    observed = np.zeros_like(sources, dtype=bool)
    for j, name in enumerate(R.DONORS):
        factor = b.ALPHA * model._global * model._renorm[name]
        for i, target in enumerate(targets):
            if target not in model.pert_index[name]:
                continue
            tcols, dcols = model.tcols[name], model.dcols[name]
            row = model.pert_index[name][target]
            sources[j, i, tcols] = factor * model.s[name][row][dcols].astype(float)
            observed[j, i, tcols] = True
    operator_mask = np.zeros(len(genes), bool)
    operator_mask[cols] = True
    staged, eligible = arms(base, sources, observed, operator_mask)
    if list(staged) != list(NAMES):
        raise ValueError("arm order drift")
    if np.any(eligible[:, np.isin(genes, targets)]):
        raise AssertionError("target coordinate was eligible")

    op = b.load_production_operator(
        model, targets, R.OPERATOR, R.OPERATOR.with_suffix(".json"),
        expected_sha256=comp["new_operator_sha256"])
    R.verify_axis(genes, cols, targets, op.axis_sha256)
    if not np.array_equal(cols, op.columns):
        raise ValueError("operator columns differ")
    cd4 = b.EXP062.CD4Prepared(
        mode="gene_norm", mean_effect=z["cd4_mean"], gene_votes=z["cd4_votes"],
        target_condition_count=z["cd4_counts"], condition_details=[])
    positions = np.flatnonzero(np.isin(genes, targets))
    effects, audit = {}, {}
    for name in NAMES:
        pre = staged[name]
        if not np.array_equal(pre[~eligible], base[~eligible]):
            raise AssertionError(f"{name}: outside-support change")
        before = np.linalg.norm(np.where(eligible, base, 0), axis=1)
        after = np.linalg.norm(np.where(eligible, pre, 0), axis=1)
        norm_error = float(np.max(np.abs(before-after)))
        if norm_error > 1e-10:
            raise AssertionError(f"{name}: eligible-norm identity {norm_error}")
        candidate, pm_audit = b.apply_perturbmatch_k562_replacement(
            pre, model, pm, targets, alpha=b.ALPHA, weights=weights)
        if name == "cap15":
            A.require_same("cap15 PM hash", b.array_sha256(candidate),
                           stage["arms"]["cap15"]["candidate_pre_sha256"])
        w, whitening_audit = b.apply_v9_whitening(candidate, op.whitener, cols)
        final, cd4_audit, union_audit = b.apply_frozen_cd4_residual(
            w, cd4, cols, positions)
        if union_audit is not None or b.H1_UNION_G is not None or b.MASS_FIX is not None:
            raise AssertionError("union or mass-fix drift")
        if name == "cap15":
            A.require_same("cap15 final hash", b.array_sha256(final),
                           comp["arms"]["cap15"]["final_sha256"])
        effects[name] = final
        audit[name] = {
            "pre_sha256": b.array_sha256(pre), "pm_sha256": b.array_sha256(candidate),
            "final_sha256": b.array_sha256(final),
            "pre_eligible_norm_max_abs_difference": norm_error,
            "pre_geometry": geometry(pre, base, operator_mask),
            "final_geometry": geometry(final, effects["cap15"], operator_mask),
            "pm": pm_audit, "whitening": whitening_audit, "cd4": cd4_audit,
        }
        print(f"constructed {name}", flush=True)
    np.savez(WORK / "effects.npz", genes=genes, targets=targets, eligible=eligible,
             **effects)
    R.write_json(WORK / "compose.json", {
        "gates_pass": True, "arms": audit, "held_sources_loaded": [],
        "h1_union": "OFF", "dose": .25,
        "eligible_coordinates": int(eligible.sum()),
        "total_target_gene_coordinates": int(eligible.size),
        "eligible_targets": int(eligible.any(1).sum()),
        "operator_sha256": R.sha(R.OPERATOR),
        "effects_sha256": R.sha(WORK / "effects.npz"),
    })


def evaluate(name):
    verify()
    if (WORK / f"metrics-{name}.json").exists():
        raise FileExistsError(f"do not overwrite {name}")
    comp = json.loads((WORK / "compose.json").read_text())
    A.require_same("composition gate", comp["gates_pass"], True)
    A.require_same("effects SHA", R.sha(WORK / "effects.npz"), comp["effects_sha256"])
    if name != "cap15":
        A.require_same("known answer", json.loads((WORK / "known-answer.json").read_text())["pass"], True)
    r, c, emit_c = E.load_reference()
    with np.load(WORK / "effects.npz", allow_pickle=False) as z:
        if not np.array_equal(z["targets"], r["targets"]):
            raise ValueError("target order drift")
        effects = R.map_axis(z[name], z["genes"], r["genes"])
    with A.evaluator_output(WORK):
        result = E.emit_and_score(name, r, c, emit_c, effects)
    if name == "cap15":
        expected = json.loads((BASE / "work/metrics-cap15.json").read_text())
        difference = {m: abs(result["raw"][m]-expected["raw"][m]) for m in E.MEMBERS}
        digests = {k: result[k] == expected[k] for k in
                   ("panel_total_digest_sha256", "emitted_count_digest_sha256")}
        passed = all(digests.values()) and all(difference[m] <= E.TOLERANCE[m] for m in E.MEMBERS)
        R.write_json(WORK / "known-answer.json", {
            "pass": passed, "digests_equal": digests, "raw_absolute_differences": difference})
        if not passed:
            raise AssertionError("cap15 reproduction failed")
    gc.collect()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "compose", *NAMES))
    args = parser.parse_args()
    started = datetime.now(timezone.utc).isoformat()
    timer = time.monotonic()
    if args.mode in ("freeze", "compose"):
        globals()[args.mode]()
    else:
        evaluate(args.mode)
    if args.mode != "freeze":
        R.write_json(WORK / f"execution-{args.mode}.json", {
            "mode": args.mode, "status": "complete", "pid": os.getpid(),
            "started_utc": started,
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.monotonic()-timer,
            "process_max_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        })


if __name__ == "__main__":
    main()
