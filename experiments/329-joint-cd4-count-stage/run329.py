"""One posterior-estimation experiment; imports production, writes only own outputs."""

import argparse
from dataclasses import replace
from datetime import datetime, timezone
import gc
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np

from estimator329 import extractor_delta, geometry, replace_deltas, validate_fit

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
PILOT = ROOT / "experiments/328-joint-cd4-estimation"
sys.path.insert(0, str(ROOT / "experiments/326-cap30-current-family"))
import run326 as A

R, E = A.R, A.E
NAMES = ("cap15", "joint", "diagonal", "shuffled_fit", "shuffled_targets")
FIT_NAMES = {"joint": "joint", "diagonal": "diagonal", "shuffled_fit": "shuffled_joint"}
CONDITIONS = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
OPERATOR_SHA = "a46fd3239a74752d5c72f5c6ae7169d354793f04f7c22b9b3fd39fe6f74b9741"


def input_paths():
    files = [HERE / n for n in ("PREREGISTRATION.md", "estimator329.py", "run329.py",
                                "report329.py", "test329.py")]
    files += [BASE / n for n in ("inputs.json", "evaluation-inputs.json", "anchor-inputs.json",
                                 "run323.py", "evaluate323.py", "report323.py")]
    files += [BASE / "work" / n for n in (
        "stage.npz", "stage.json", "compose.json", "effects.npz",
        "whitener_h1_pruned.npz", "whitener_h1_pruned.json", "reference.npz",
        "reference.json", "reference_de.parquet", "h1-anchors.json",
        "metrics-cap15.json", "stats-cap15.npz", "per-target-cap15.csv")]
    files += [ROOT / "experiments/326-cap30-current-family/run326.py"]
    files += [PILOT / n for n in ("mixture328.py", "inputs.json", "selection.json",
                                  "recovery-inputs.json", "NUMERICAL-RECOVERY.md")]
    files += [PILOT / "work/recovery" / n for n in (
        "results.json", "fit-joint.json", "fit-diagonal.json", "fit-shuffled_joint.json")]
    files += [R.CD4 / n for n in ("gene_map.csv", "target_row_map.csv", "manifest.json")]
    for condition in CONDITIONS:
        files += [R.CD4 / condition / n for n in (
            "log_fc.npy", "lfc_se.npy", "delta.npy", "b_ctrl.npy", "genes.txt",
            "perts.txt", "meta.json")]
    return files


def freeze():
    path = HERE / "inputs.json"
    if path.exists():
        raise FileExistsError(path)
    R.write_json(path, {
        "inputs_sha256": {str(p.relative_to(ROOT)): R.sha(p) for p in input_paths()},
        "arms": NAMES, "batch_size": 8192, "seed": 20260909,
        "classification": "reused H1 development after disclosed source numerical amendment"})


def verify():
    manifest = json.loads((HERE / "inputs.json").read_text())
    if set(manifest["inputs_sha256"]) != {str(p.relative_to(ROOT)) for p in input_paths()}:
        raise ValueError("329 input allowlist drift")
    for path, expected in manifest["inputs_sha256"].items():
        A.require_same("329 input "+path, R.sha(ROOT / path), expected)
    own = [str(p.relative_to(ROOT)) for p in input_paths() if p.parent == HERE]
    own += [str((HERE / "inputs.json").relative_to(ROOT))]
    for path in own:
        subprocess.run(["git", "ls-files", "--error-unmatch", path], cwd=ROOT,
                       check=True, stdout=subprocess.DEVNULL)
    if subprocess.check_output(["git", "status", "--porcelain", "--", *own], cwd=ROOT, text=True):
        raise ValueError("commit 329 freeze before measurement")
    recovery = json.loads((PILOT / "work/recovery/results.json").read_text())
    A.require_same("source recovery gate", recovery["status"], "EARNS_H1_COUNT_STAGE_REGISTRATION")
    if not all(recovery["gates"].values()):
        raise ValueError("source gates failed")
    E.verify()


def load_cd4(b, targets, official_genes):
    """Read only donor coefficients here; reference outcomes are not loaded."""
    banks = [b.EXP062.Bank.load(n, R.CD4 / n, require_cd4_scale=True) for n in CONDITIONS]
    genes = np.array(banks[0].genes)
    if len(set(genes)) != len(genes) or len(set(targets)) != len(targets):
        raise ValueError("ambiguous axes")
    target_pos = {t: i for i, t in enumerate(targets)}
    shape = (len(targets), len(genes), 3)
    y, se = np.full(shape, np.nan, np.float32), np.full(shape, np.nan, np.float32)
    delta = np.zeros(shape, np.float32)
    present = np.zeros((len(targets), 3), bool)
    maps = []
    for j, bank in enumerate(banks):
        if bank.genes != banks[0].genes or not set(bank.perts) <= set(targets):
            raise ValueError("compact-bank axis or target drift")
        rows = np.array([target_pos[t] for t in bank.perts])
        maps.append(rows)
        x = np.load(bank.path / "log_fc.npy", allow_pickle=False)
        s = np.load(bank.path / "lfc_se.npy", allow_pickle=False)
        d = np.load(bank.delta_path, allow_pickle=False)
        if x.shape != bank.shape or s.shape != bank.shape:
            raise ValueError("compact coefficient shape drift")
        if not np.array_equal(extractor_delta(x), d):
            raise AssertionError(f"{bank.name}: identity coefficient conversion failed")
        y[rows, :, j], se[rows, :, j], delta[rows, :, j] = x, s, d
        present[rows, j] = True
    eligible_genes = np.isin(genes, official_genes) & ~np.isin(genes, targets)
    return banks, genes, y, se, delta, present, maps, eligible_genes


def prepare(b, name, cube, banks, maps, targets, genes):
    replaced = []
    for j, (bank, rows) in enumerate(zip(banks, maps)):
        path = WORK / "banks" / name / bank.name / "delta.npy"
        if path.exists():
            raise FileExistsError(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.save(path, cube[rows, :, j], allow_pickle=False)
        replaced.append(replace(bank, delta_path=path))
    prepared = b.EXP062.prepare_cd4_mean(
        replaced, targets, genes, set(targets), b.CD4_MODE,
        chunk_rows=128, max_rss_gib=22)
    cd4 = b.freeze_complete_mean3(prepared)
    if np.count_nonzero(cd4.target_condition_count > 0) != 80:
        raise AssertionError("CD4 carrying coverage drift")
    if np.count_nonzero(cd4.target_condition_count == 3) != 77:
        raise AssertionError("CD4 complete coverage drift")
    incomplete = cd4.target_condition_count != 3
    if np.any(cd4.mean_effect[incomplete]) or np.any(cd4.gene_votes[incomplete]):
        raise AssertionError("incomplete CD4 row entered residual")
    return cd4


def compose():
    verify()
    if (WORK / "effects.npz").exists() or (WORK / "banks").exists():
        raise FileExistsError("do not overwrite construction; inspect any failed attempt")
    WORK.mkdir(parents=True, exist_ok=True)
    with np.load(BASE / "work/stage.npz", allow_pickle=False) as z:
        stage = {k: z[k] for k in ("genes", "targets", "columns", "cap15",
                                  "cd4_mean", "cd4_votes", "cd4_counts")}
    genes, targets, cols = stage["genes"], stage["targets"], stage["columns"]
    comp = json.loads((BASE / "work/compose.json").read_text())
    stage_info = json.loads((BASE / "work/stage.json").read_text())
    A.require_same("323 stage gate", stage_info["gates_pass"], True)
    A.require_same("323 composition gate", comp["gates_pass"], True)
    A.require_same("operator SHA", R.sha(R.OPERATOR), OPERATOR_SHA)
    b = R.builder("cap15")
    if (b.ALPHA != .982 or b.BETA != 0 or b.CD4_DEFLATE_K != 5 or
            b.CD4_TARGET_F != .34 or b.H1_UNION_G is not None or b.MASS_FIX is not None):
        raise ValueError("composition settings drift")
    model = SimpleNamespace(
        names=list(R.DONORS), target_genes=genes,
        tcols={n: np.flatnonzero(np.isin(genes, (
            ROOT / f"data/donors/{n}/genes.txt").read_text().splitlines())) for n in R.DONORS})
    op = b.load_production_operator(
        model, targets, R.OPERATOR, R.OPERATOR.with_suffix(".json"),
        expected_sha256=OPERATOR_SHA)
    R.verify_axis(genes, cols, targets, op.axis_sha256)
    if not np.array_equal(cols, op.columns):
        raise ValueError("operator columns drift")
    w, whitening_audit = b.apply_v9_whitening(stage["cap15"], op.whitener, cols)
    del op, model
    banks, cgenes, y, se, original, present, maps, eligible_genes = load_cd4(b, targets, genes)
    positions = np.flatnonzero(np.isin(genes, targets))
    mask = ~np.isin(genes, targets)
    effects, audits, baseline_cd4 = {}, {}, None
    eligible = None
    for name in NAMES[:-1]:
        if name == "cap15":
            cube = original.copy()
            for j, rows in enumerate(maps):
                cube[rows, :, j] = extractor_delta(y[rows, :, j])
            if not np.array_equal(cube, original):
                raise AssertionError("identity donor cube mismatch")
        else:
            fit = json.loads((PILOT / "work/recovery" /
                              f"fit-{FIT_NAMES[name]}.json").read_text())
            validate_fit(fit, diagonal=name == "diagonal")
            cube, current_eligible = replace_deltas(
                original, y, se, present, eligible_genes, fit, diagonal=name == "diagonal")
            if eligible is not None and not np.array_equal(current_eligible, eligible):
                raise AssertionError("posterior support differs across fits")
            eligible = current_eligible
        cd4 = prepare(b, name, cube, banks, maps, targets, genes)
        if name == "cap15":
            for attr, key in (("mean_effect", "cd4_mean"), ("gene_votes", "cd4_votes"),
                              ("target_condition_count", "cd4_counts")):
                if not np.array_equal(getattr(cd4, attr), stage[key]):
                    raise AssertionError(f"identity CD4 preparation: {attr}")
            baseline_cd4 = cd4.mean_effect.copy()
        else:
            if (not np.array_equal(cd4.gene_votes, stage["cd4_votes"]) or
                    not np.array_equal(cd4.target_condition_count, stage["cd4_counts"])):
                raise AssertionError("candidate CD4 support drift")
        deflated, da = b.deflate_cd4_residual(cd4, 5)
        incomplete = cd4.target_condition_count != 3
        if np.any(deflated.mean_effect[incomplete]) or np.any(deflated.gene_votes[incomplete]):
            raise AssertionError("deflation expanded incomplete-row support")
        if np.any(deflated.mean_effect[:, positions]):
            raise AssertionError("deflation spread excluded target coordinates")
        del deflated
        final, ca, ua = b.apply_frozen_cd4_residual(w, cd4, cols, positions)
        if ca["deflation"] != da:
            raise AssertionError("independent deflation audit mismatch")
        if ua is not None or not np.isfinite(final).all():
            raise AssertionError("union drift or nonfinite effect")
        if np.any(cd4.mean_effect[:, positions]) or np.any(cd4.gene_votes[:, positions]):
            raise AssertionError("panel coordinate entered CD4")
        if name == "cap15":
            A.require_same("identity final hash", b.array_sha256(final),
                           comp["arms"]["cap15"]["final_sha256"])
            A.require_same("identity deflation basis", ca["deflation"]["basis_u_k_sha256"],
                           "31fad06fe6de4f190265903f1846aa2a74dc642af98b52eeb5b70459ed435a4d")
        effects[name] = final
        audits[name] = {
            "final_sha256": b.array_sha256(final),
            "raw_delta_cube_sha256": b.array_sha256(cube),
            "source_geometry": geometry(
                cube.transpose(0, 2, 1).reshape(len(targets), -1),
                original.transpose(0, 2, 1).reshape(len(targets), -1)),
            "prepared_cd4_geometry": geometry(cd4.mean_effect, baseline_cd4, mask),
            "composed_geometry": geometry(final, effects["cap15"], mask),
            "cd4": ca, "condition_details": cd4.condition_details,
        }
        print(f"constructed {name}", flush=True)
    permutation = R.derangement(len(targets), seed=20260909)
    effects["shuffled_targets"] = effects["joint"][permutation]
    np.savez(WORK / "effects.npz", genes=genes, targets=targets, eligible=eligible,
             source_genes=cgenes, **effects, permutation=permutation)
    R.write_json(WORK / "compose.json", {
        "gates_pass": True, "arms": audits, "whitening": whitening_audit,
        "h1_union": "OFF", "recipient_h1_response_reads_in_construction": False,
        "posterior_eligible_triples": int(eligible.sum()),
        "posterior_eligible_targets": int(eligible.any(1).sum()),
        "posterior_eligible_genes": int(eligible.any(0).sum()),
        "raw_fallback_triples": int(eligible.size-eligible.sum()),
        "total_aligned_triples": int(eligible.size),
        "operator_sha256": OPERATOR_SHA,
        "effects_sha256": R.sha(WORK / "effects.npz"),
        "bank_sha256": {str(p.relative_to(ROOT)): R.sha(p)
                        for p in sorted((WORK / "banks").glob("*/*/delta.npy"))},
    })


def evaluate(name):
    verify()
    if (WORK / f"metrics-{name}.json").exists():
        raise FileExistsError(f"refusing to rescore {name}")
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
        differences = {m: abs(result["raw"][m]-expected["raw"][m]) for m in E.MEMBERS}
        digests = {k: result[k] == expected[k] for k in
                   ("panel_total_digest_sha256", "emitted_count_digest_sha256")}
        passed = all(digests.values()) and all(
            differences[m] <= E.TOLERANCE[m] for m in E.MEMBERS)
        R.write_json(WORK / "known-answer.json", {
            "pass": passed, "digests_equal": digests, "raw_absolute_differences": differences})
        if not passed:
            raise AssertionError("identity count-level reproduction failed")


def execute(mode):
    started, timer = datetime.now(timezone.utc).isoformat(), time.monotonic()
    try:
        if mode == "compose":
            compose()
        elif mode == "report":
            import report329
            report329.main()
        else:
            evaluate(mode)
    except BaseException as exc:
        R.write_json(WORK / f"execution-{mode}.json", {
            "status": "failed", "mode": mode, "error": repr(exc),
            "started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.monotonic()-timer})
        raise
    R.write_json(WORK / f"execution-{mode}.json", {
        "status": "complete", "mode": mode, "pid": os.getpid(),
        "started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic()-timer,
        "process_max_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    gc.collect()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "compose", *NAMES, "report", "all"))
    mode = parser.parse_args().mode
    if mode == "freeze":
        freeze()
    elif mode == "all":
        for step in ("compose", *NAMES, "report"):
            execute(step)
    else:
        execute(mode)
