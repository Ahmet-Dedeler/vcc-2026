"""Frozen source-background normalization experiment; no posterior model used."""

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
from types import SimpleNamespace

import h5py
import numpy as np
import pandas as pd

from scales330 import gene_scale, permute_scales, scale_override

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/329-joint-cd4-count-stage"))
import run329 as C
from estimator329 import extractor_delta, geometry

R, E, A = C.R, C.E, C.A
BASE, PILOT = C.BASE, C.PILOT
SOURCE = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"
CONDITIONS = C.CONDITIONS
NAMES = ("cap15", "background1000", "background80", "shuffled_scale", "shuffled_targets")


def input_paths():
    files = [HERE / n for n in (
        "PREREGISTRATION.md", "scales330.py", "run330.py", "report330.py", "test330.py")]
    files += [C.HERE / n for n in ("run329.py", "estimator329.py")]
    files += [BASE / n for n in (
        "inputs.json", "evaluation-inputs.json", "anchor-inputs.json",
        "run323.py", "evaluate323.py", "report323.py")]
    files += [BASE / "work" / n for n in (
        "stage.npz", "stage.json", "compose.json", "effects.npz",
        "whitener_h1_pruned.npz", "whitener_h1_pruned.json",
        "reference.npz", "reference.json", "reference_de.parquet", "h1-anchors.json",
        "metrics-cap15.json")]
    files += [ROOT / "experiments/326-cap30-current-family/run326.py", SOURCE]
    files += [PILOT / n for n in ("inputs.json", "selection.json", "work/source.npz")]
    files += [R.CD4 / "gene_map.csv", R.CD4 / "target_row_map.csv"]
    for name in CONDITIONS:
        files += [R.CD4 / name / n for n in ("delta.npy", "b_ctrl.npy", "genes.txt",
                                            "perts.txt", "meta.json")]
    return files


def selection():
    """Identifiers only; selects no outcome or source value."""
    previous = json.loads((PILOT / "selection.json").read_text())
    targets = previous["targets"][:1000]
    if any(t["split"] != "fit" for t in targets):
        raise ValueError("not the frozen 328 fit cohort")
    excluded = set(previous["excluded_symbols"])
    ids = set(previous["excluded_ids"])
    if any(t["name"].casefold() in excluded or t["id"] in ids for t in targets):
        raise ValueError("panel target in background fit")
    genes = (R.CD4 / CONDITIONS[0] / "genes.txt").read_text().splitlines()
    official = set(pd.read_csv(ROOT / "data/challenge/gene_names.csv").gene_name)
    frame = pd.read_csv(R.CD4 / "gene_map.csv").set_index("gene_name")
    columns = [(i, g, int(frame.loc[g, "source_var_index"])) for i, g in enumerate(genes)
               if g in official and g.casefold() not in excluded]
    if len({g for _, g, _ in columns}) != len(columns):
        raise ValueError("nonunique selected response axis")
    return {
        "targets": targets, "conditions": previous["conditions"],
        "genes": [{"compact_column": i, "name": g, "source_column": j} for i, g, j in columns],
        "excluded_symbols": sorted(excluded),
        "classification": "non-panel source background; reused 328 fit-target identities"}


def freeze():
    path = HERE / "inputs.json"
    if path.exists():
        raise FileExistsError(path)
    R.write_json(path, {
        "inputs_sha256": {str(p.relative_to(ROOT)): R.sha(p) for p in input_paths()},
        "selection": selection(), "arms": NAMES, "seed": 20260909})


def verify():
    manifest = json.loads((HERE / "inputs.json").read_text())
    if set(manifest["inputs_sha256"]) != {str(p.relative_to(ROOT)) for p in input_paths()}:
        raise ValueError("330 input allowlist drift")
    for path, digest in manifest["inputs_sha256"].items():
        A.require_same("330 input "+path, R.sha(ROOT / path), digest)
    if manifest["selection"] != selection():
        raise ValueError("background selection drift")
    own = [str(p.relative_to(ROOT)) for p in input_paths() if p.parent == HERE]
    own += [str((HERE / "inputs.json").relative_to(ROOT))]
    for path in own:
        subprocess.run(["git", "ls-files", "--error-unmatch", path], cwd=ROOT,
                       check=True, stdout=subprocess.DEVNULL)
    if subprocess.check_output(["git", "status", "--porcelain", "--", *own], cwd=ROOT, text=True):
        raise ValueError("commit 330 before measurements")
    E.verify()


def extract():
    verify()
    if (WORK / "scales.npz").exists():
        raise FileExistsError("background already extracted")
    meta = json.loads((HERE / "inputs.json").read_text())["selection"]
    targets = meta["targets"]
    source_cols = np.array([g["source_column"] for g in meta["genes"]])
    compact_cols = np.array([g["compact_column"] for g in meta["genes"]])
    old = json.loads((PILOT / "selection.json").read_text())
    selected_position = {g["name"]: i for i, g in enumerate(meta["genes"])}
    old_positions, new_positions = [], []
    for i, gene in enumerate(old["genes"]):
        if gene["name"] in selected_position:
            old_positions.append(i)
            new_positions.append(selected_position[gene["name"]])
    with np.load(PILOT / "work/source.npz", allow_pickle=False) as z:
        prior_y = z["y"][:1000]
        prior_targets = z["targets"][:1000]
    if not np.array_equal(prior_targets, [t["id"] for t in targets]):
        raise ValueError("prior source target order mismatch")
    scales, short, audit = [], [], {}
    with h5py.File(SOURCE, "r") as source:
        layer = source["layers/log_fc"]
        for j, condition in enumerate(meta["conditions"]):
            y = np.empty((1000, len(source_cols)), np.float32)
            for i, target in enumerate(targets):
                # One full 10,282-column row is small; reorder only after the read.
                y[i] = layer[target["rows"][condition], :][source_cols]
            expected = prior_y[:, old_positions, j].astype(np.float32)
            observed = y[:, new_positions]
            if not np.array_equal(observed, expected, equal_nan=True):
                raise AssertionError("broader extraction failed the 328 known-answer")
            finite = float(np.isfinite(y).mean())
            if finite < .99:
                raise ValueError(f"{condition}: insufficient finite source coefficients")
            x = extractor_delta(y)
            full, info = gene_scale(x)
            small, small_info = gene_scale(x[:80])
            scales.append(full)
            short.append(small)
            audit[condition] = {
                "finite_fraction": finite, "known_answer_genes": len(old_positions),
                "known_answer_coefficients_exact": True, "background1000": info,
                "background80": small_info}
            print(f"extracted and fitted scales: {condition}, {x.shape}", flush=True)
    WORK.mkdir(parents=True, exist_ok=True)
    np.savez(WORK / "scales.npz", columns=compact_cols, genes=[g["name"] for g in meta["genes"]],
             background1000=scales, background80=short, shuffled_scale=permute_scales(scales))
    R.write_json(WORK / "source.json", {
        "gates_pass": True, "conditions": audit, "n_targets": 1000,
        "n_response_genes": len(source_cols), "input_manifest_sha256": R.sha(HERE / "inputs.json"),
        "scales_sha256": R.sha(WORK / "scales.npz"), "recipient_response_reads": False})


def compose():
    verify()
    if (WORK / "effects.npz").exists():
        raise FileExistsError("do not overwrite composition")
    source = json.loads((WORK / "source.json").read_text())
    A.require_same("source gate", source["gates_pass"], True)
    A.require_same("scales SHA", R.sha(WORK / "scales.npz"), source["scales_sha256"])
    with np.load(WORK / "scales.npz", allow_pickle=False) as z:
        columns = z["columns"]
        scale_arms = {name: z[name] for name in NAMES[1:4]}
    with np.load(BASE / "work/stage.npz", allow_pickle=False) as z:
        stage = {k: z[k] for k in ("genes", "targets", "columns", "cap15",
                                  "cd4_mean", "cd4_votes", "cd4_counts")}
    genes, targets, cols = stage["genes"], stage["targets"], stage["columns"]
    old = json.loads((BASE / "work/compose.json").read_text())
    A.require_same("323 composition gate", old["gates_pass"], True)
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
        expected_sha256=C.OPERATOR_SHA)
    R.verify_axis(genes, cols, targets, op.axis_sha256)
    if not np.array_equal(cols, op.columns):
        raise ValueError("operator columns drift")
    w, wa = b.apply_v9_whitening(stage["cap15"], op.whitener, cols)
    del op, model
    positions = np.flatnonzero(np.isin(genes, targets))
    mask = ~np.isin(genes, targets)
    effects, audit = {}, {}
    for name in NAMES[:-1]:
        replacement = None if name == "cap15" else dict(zip(CONDITIONS, scale_arms[name]))
        changes = {}
        with scale_override(b.EXP062, replacement, columns, changes):
            _, cd4 = b.prepare_frozen_cd4(
                R.CD4, targets, genes, targets, chunk_rows=128, max_rss_gib=22)
        for attr, key in (("gene_votes", "cd4_votes"), ("target_condition_count", "cd4_counts")):
            if not np.array_equal(getattr(cd4, attr), stage[key]):
                raise AssertionError("CD4 support drift: "+attr)
        if (np.count_nonzero(cd4.target_condition_count > 0) != 80 or
                np.count_nonzero(cd4.target_condition_count == 3) != 77):
            raise AssertionError("coverage drift")
        if name == "cap15" and not np.array_equal(cd4.mean_effect, stage["cd4_mean"]):
            raise AssertionError("baseline CD4 preparation mismatch")
        deflated, da = b.deflate_cd4_residual(cd4, 5)
        incomplete = cd4.target_condition_count != 3
        if (np.any(cd4.mean_effect[incomplete]) or np.any(deflated.mean_effect[incomplete]) or
                np.any(cd4.mean_effect[:, positions]) or np.any(deflated.mean_effect[:, positions])):
            raise AssertionError("excluded support expanded")
        del deflated
        final, ca, ua = b.apply_frozen_cd4_residual(w, cd4, cols, positions)
        if ua is not None or ca["deflation"] != da or not np.isfinite(final).all():
            raise AssertionError("unexpected union, deflation or finite-state failure")
        if name == "cap15":
            A.require_same("baseline final SHA", b.array_sha256(final),
                           old["arms"]["cap15"]["final_sha256"])
            A.require_same("baseline U5 SHA", da["basis_u_k_sha256"],
                           "31fad06fe6de4f190265903f1846aa2a74dc642af98b52eeb5b70459ed435a4d")
        effects[name] = final
        audit[name] = {
            "scale_override": changes, "cd4": ca, "condition_details": cd4.condition_details,
            "prepared_cd4_geometry": geometry(cd4.mean_effect, stage["cd4_mean"], mask),
            "composed_geometry": geometry(final, effects["cap15"], mask),
            "final_sha256": b.array_sha256(final)}
        print("constructed "+name, flush=True)
    permutation = R.derangement(len(targets), seed=20260909)
    effects["shuffled_targets"] = effects["background1000"][permutation]
    np.savez(WORK / "effects.npz", genes=genes, targets=targets, **effects, permutation=permutation)
    R.write_json(WORK / "compose.json", {
        "gates_pass": True, "arms": audit, "whitening": wa, "h1_union": "OFF",
        "recipient_response_reads_in_construction": False,
        "operator_sha256": C.OPERATOR_SHA, "scales_sha256": R.sha(WORK / "scales.npz"),
        "effects_sha256": R.sha(WORK / "effects.npz")})


def evaluate(name):
    verify()
    if (WORK / f"metrics-{name}.json").exists():
        raise FileExistsError("refusing to overwrite a score")
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
        passed = all(digests.values()) and all(differences[m] <= E.TOLERANCE[m] for m in E.MEMBERS)
        R.write_json(WORK / "known-answer.json", {
            "pass": passed, "digests_equal": digests, "raw_absolute_differences": differences})
        if not passed:
            raise AssertionError("baseline count reproduction failed")


def execute(mode):
    started, timer = datetime.now(timezone.utc).isoformat(), time.monotonic()
    try:
        if mode in ("extract", "compose"):
            globals()[mode]()
        elif mode == "report":
            import report330
            report330.main()
        else:
            evaluate(mode)
    except BaseException as exc:
        R.write_json(WORK / f"execution-{mode}.json", {
            "status": "failed", "mode": mode, "error": repr(exc), "started_utc": started,
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.monotonic()-timer})
        raise
    R.write_json(WORK / f"execution-{mode}.json", {
        "status": "complete", "mode": mode, "pid": os.getpid(), "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic()-timer,
        "process_max_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    gc.collect()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "extract", "compose", *NAMES, "report", "all"))
    mode = parser.parse_args().mode
    if mode == "freeze":
        freeze()
    elif mode == "all":
        for step in ("extract", "compose", *NAMES, "report"):
            execute(step)
    else:
        execute(mode)
