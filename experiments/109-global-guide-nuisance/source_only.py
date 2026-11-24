#!/usr/bin/env python3
"""Fail-closed source-only orchestrator for experiment 109.

This command can preflight or generate experiment-083 compact K562 caches,
run the preregistered experiment-109 tune/test/null procedure, and atomically
publish a source-only bundle. It never imports a held-context evaluator,
opens held-context outcomes, generates counts, or submits predictions.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import fit_freeze
import model
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXP083 = ROOT / "experiments/083-guide-latent-decontamination"
DEFAULT_SOURCE_H5AD = (
    ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
)
DEFAULT_COMPACT_K562 = ROOT / "data/donors/k562_gwps"
DEFAULT_PANEL_CSV = ROOT / "data/challenge/pert_counts.csv"
DEFAULT_PM_PERTS = ROOT / "data/external/perturbmatch/perts.txt"
DEFAULT_PM_GENES = ROOT / "data/external/perturbmatch/genes.txt"
DEFAULT_WORK = HERE / "work"
DEFAULT_GUIDE_CACHE = DEFAULT_WORK / "guide-cache"
DEFAULT_FEATURE_CACHE = DEFAULT_WORK / "feature-cache"
DEFAULT_OUTPUT = DEFAULT_WORK / "frozen-source"
CACHE_METADATA = {
    "guide": ("metadata.json", *fit_freeze.GUIDE_CACHE_OUTPUTS),
    "feature": ("seed-feature-provenance.json", *fit_freeze.FEATURE_CACHE_OUTPUTS),
}


def _write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _savez(path: Path, **values: object) -> None:
    with path.open("wb") as handle:
        np.savez_compressed(handle, **values)


def _cache_state(path: Path, kind: str) -> dict[str, object]:
    required = CACHE_METADATA[kind]
    if not path.exists():
        return {"state": "absent", "path": str(path), "missing": list(required)}
    if not path.is_dir():
        return {
            "state": "invalid",
            "path": str(path),
            "missing": list(required),
            "reason": "cache path is not a directory",
        }
    missing = [name for name in required if not (path / name).is_file()]
    state = "complete" if not missing else "partial"
    return {"state": state, "path": str(path), "missing": missing}


def _verify_hash(path: Path, expected: str, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} is missing: {path}")
    observed = fit_freeze.sha256_file(path)
    if observed != expected:
        raise ValueError(f"{label} provenance mismatch: {observed} != {expected}")


def _verify_generation_inputs(
    *,
    source_h5ad: Path,
    compact_k562: Path,
    panel_csv: Path,
    perturbmatch_perts: Path,
    perturbmatch_genes: Path,
    seed_table: Path | None,
    allow_seed_download: bool,
    hash_h5ad: bool,
) -> dict[str, object]:
    source_axes = fit_freeze.load_canonical_source_axes()
    source = source_axes["source"]
    if not source_h5ad.is_file():
        raise FileNotFoundError(f"K562 source H5AD is missing: {source_h5ad}")
    observed_bytes = source_h5ad.stat().st_size
    expected_bytes = int(source["h5ad_bytes"])
    if observed_bytes != expected_bytes:
        raise ValueError(
            f"K562 source byte count changed: {observed_bytes} != {expected_bytes}"
        )
    h5ad_sha256 = None
    if hash_h5ad:
        h5ad_sha256 = fit_freeze.sha256_file(source_h5ad)
        if h5ad_sha256 != source["h5ad_sha256_from_pinned_source_audit"]:
            raise ValueError("K562 source H5AD hash differs from the pinned audit")

    _verify_hash(panel_csv, str(source["panel_csv_sha256"]), "panel CSV")
    _verify_hash(
        perturbmatch_perts,
        str(source["perturbmatch_perts_sha256"]),
        "PerturbMatch target axis",
    )
    _verify_hash(
        perturbmatch_genes,
        str(source["perturbmatch_genes_sha256"]),
        "PerturbMatch gene axis",
    )

    artifact_manifest = json.loads(
        (EXP083 / "artifact_manifest.json").read_text()
    )
    compact_hashes = artifact_manifest["donors"]["k562_gwps"]
    for filename, expected in compact_hashes.items():
        _verify_hash(
            compact_k562 / filename,
            str(expected),
            f"compact K562 {filename}",
        )

    seed_receipt: dict[str, object]
    if seed_table is not None:
        expected_seed = source["seed_metadata"]
        if seed_table.stat().st_size != int(expected_seed["compressed_bytes"]):
            raise ValueError("local seed table byte count differs from pinned source")
        _verify_hash(
            seed_table,
            str(expected_seed["compressed_sha256"]),
            "local seed table",
        )
        seed_receipt = {
            "mode": "local",
            "path": str(seed_table),
            "sha256": expected_seed["compressed_sha256"],
        }
    elif allow_seed_download:
        seed_receipt = {
            "mode": "pinned_https_download",
            "source": source["seed_metadata"]["source"],
            "sha256": source["seed_metadata"]["compressed_sha256"],
        }
    else:
        raise FileNotFoundError(
            "feature-cache generation needs --seed-table or "
            "--allow-seed-download"
        )
    return {
        "source_h5ad": {
            "path": str(source_h5ad),
            "bytes": observed_bytes,
            "sha256": h5ad_sha256,
            "full_hash_checked": hash_h5ad,
        },
        "compact_k562": {
            "path": str(compact_k562),
            "files": dict(compact_hashes),
        },
        "seed_table": seed_receipt,
    }


def _extract_command(
    args: argparse.Namespace,
    *,
    output: Path,
) -> list[str]:
    return [
        sys.executable,
        str(EXP083 / "extract_guide_effects.py"),
        "--h5ad",
        str(args.source_h5ad),
        "--compact-k562",
        str(args.compact_k562),
        "--panel-csv",
        str(args.panel_csv),
        "--output",
        str(output),
        "--max-rss-gib",
        str(args.max_rss_gib),
        "--guarded",
    ]


def _feature_command(
    args: argparse.Namespace,
    *,
    guide_cache: Path,
    output: Path,
) -> list[str]:
    command = [
        sys.executable,
        str(EXP083 / "build_seed_features.py"),
        "--guide-ids",
        str(guide_cache / "guide_ids.txt"),
        "--output",
        str(output),
    ]
    if args.seed_table is not None:
        command.extend(("--source", str(args.seed_table)))
    return command


def preflight(args: argparse.Namespace, *, hash_h5ad: bool = False) -> dict[str, object]:
    protocol = fit_freeze.validate_preregistration()
    guide_state = _cache_state(args.guide_cache, "guide")
    feature_state = _cache_state(args.feature_cache, "feature")
    errors: list[str] = []
    commands: list[list[str]] = []
    generation_inputs: dict[str, object] | None = None
    managed_paths = {
        "guide_cache": args.guide_cache,
        "feature_cache": args.feature_cache,
        "output": args.output,
    }
    for left_name, left in managed_paths.items():
        for right_name, right in managed_paths.items():
            if left_name >= right_name:
                continue
            if left == right or left in right.parents or right in left.parents:
                errors.append(
                    f"{left_name} and {right_name} paths overlap: {left} / {right}"
                )

    for state in (guide_state, feature_state):
        if state["state"] in {"partial", "invalid"}:
            errors.append(
                f"{state['path']}: {state['state']} cache; missing={state['missing']}"
            )

    caches_complete = (
        guide_state["state"] == "complete"
        and feature_state["state"] == "complete"
    )
    cache_manifest: dict[str, object] | None = None
    if caches_complete:
        try:
            cache_manifest, _ = fit_freeze.build_source_input_manifest(
                args.guide_cache,
                args.feature_cache,
            )
        except (OSError, ValueError) as exc:
            errors.append(f"compact-cache verification failed: {exc}")
    elif args.cache_mode == "verify":
        errors.append("cache mode is verify but both authenticated caches are absent")
    elif not errors:
        try:
            generation_inputs = _verify_generation_inputs(
                source_h5ad=args.source_h5ad,
                compact_k562=args.compact_k562,
                panel_csv=args.panel_csv,
                perturbmatch_perts=args.perturbmatch_perts,
                perturbmatch_genes=args.perturbmatch_genes,
                seed_table=args.seed_table,
                allow_seed_download=args.allow_seed_download,
                hash_h5ad=hash_h5ad,
            )
            if guide_state["state"] == "absent":
                commands.append(_extract_command(args, output=args.guide_cache))
            if feature_state["state"] == "absent":
                commands.append(
                    _feature_command(
                        args,
                        guide_cache=args.guide_cache,
                        output=args.feature_cache,
                    )
                )
        except (OSError, ValueError) as exc:
            errors.append(str(exc))

    if args.output.exists():
        errors.append(f"source output already exists; refusing overwrite: {args.output}")

    return {
        "experiment": fit_freeze.EXPERIMENT,
        "status": "READY" if not errors else "BLOCKED",
        "protocol_status": protocol["status"],
        "protocol_sha256": fit_freeze.sha256_file(fit_freeze.PROTOCOL_PATH),
        "cache_mode": args.cache_mode,
        "guide_cache": guide_state,
        "feature_cache": feature_state,
        "caches_authenticated": cache_manifest is not None,
        "cache_manifest_preview": cache_manifest,
        "generation_inputs": generation_inputs,
        "planned_commands": commands,
        "output": str(args.output),
        "held_evaluation_enabled": False,
        "held_outcomes_opened": False,
        "errors": errors,
    }


def _run_command(command: Sequence[str]) -> None:
    print(json.dumps({"source_only_command": list(command)}), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def _start_memory_guard(
    args: argparse.Namespace,
    *,
    state_dir: Path,
) -> tuple[subprocess.Popen[bytes], dict[str, Path]]:
    """Start the pinned guard before any cache generation or source fitting."""
    guard = fit_freeze.validate_preregistration()["source_only_runner"][
        "memory_guard"
    ]
    if fit_freeze.sha256_file(fit_freeze.MEMORY_GUARD_PATH) != (
        fit_freeze.MEMORY_GUARD_SHA256
    ):
        raise ValueError("pinned memory-guard script changed")
    state_dir.mkdir(parents=True, exist_ok=False)
    paths = {
        "ready": state_dir / "ready.json",
        "stop": state_dir / "stop",
        "report": state_dir / "report.json",
    }
    command = [
        sys.executable,
        str(fit_freeze.MEMORY_GUARD_PATH),
        "--monitor-pid",
        str(os.getpid()),
        "--ready-file",
        str(paths["ready"]),
        "--stop-file",
        str(paths["stop"]),
        "--report-file",
        str(paths["report"]),
        "--max-rss-gib",
        str(guard["max_rss_gib"]),
        "--min-free-percent",
        str(guard["min_free_percent"]),
        "--poll-seconds",
        str(guard["poll_seconds"]),
        "--quiet",
    ]
    monitor = subprocess.Popen(command, cwd=ROOT)
    deadline = time.monotonic() + 15.0
    while not paths["ready"].is_file():
        if monitor.poll() is not None:
            raise RuntimeError("memory guard exited before becoming ready")
        if time.monotonic() >= deadline:
            monitor.terminate()
            raise RuntimeError("memory guard did not become ready")
        time.sleep(0.05)
    ready = json.loads(paths["ready"].read_text())
    if (
        ready.get("status") != "active"
        or int(ready.get("root_pid", -1)) != os.getpid()
        or int(ready.get("monitor_pid", -1)) != monitor.pid
    ):
        monitor.terminate()
        raise RuntimeError("memory guard ready receipt is invalid")
    return monitor, paths


def _finish_memory_guard(
    monitor: subprocess.Popen[bytes],
    paths: Mapping[str, Path],
    *,
    destination: Path,
    source_manifest_content_sha256: str,
    source_cache_snapshot_sha256: str,
) -> dict[str, object]:
    paths["stop"].touch(exist_ok=False)
    returncode = monitor.wait(timeout=30)
    if returncode != 0:
        raise RuntimeError(f"memory guard exited with status {returncode}")
    report = json.loads(paths["report"].read_text())
    receipt = {
        "schema_version": fit_freeze.MEMORY_GUARD_RECEIPT_SCHEMA_VERSION,
        "experiment": fit_freeze.EXPERIMENT,
        "status": report.get("status"),
        "scope": fit_freeze.validate_preregistration()["source_only_runner"][
            "memory_guard"
        ]["scope"],
        "guard_script": str(fit_freeze.MEMORY_GUARD_PATH.relative_to(ROOT)),
        "guard_script_sha256": fit_freeze.MEMORY_GUARD_SHA256,
        "root_pid": report.get("root_pid"),
        "monitor_pid": report.get("monitor_pid"),
        "max_rss_gib": report.get("max_rss_gib"),
        "min_free_percent": report.get("min_free_percent"),
        "poll_seconds": report.get("poll_seconds"),
        "peak_rss_gib": report.get("peak_rss_gib"),
        "samples": report.get("samples"),
        "last_process_count": report.get("last_process_count"),
        "minimum_free_percent": report.get("minimum_free_percent"),
        "elapsed_seconds": report.get("elapsed_seconds"),
        "source_manifest_content_sha256": source_manifest_content_sha256,
        "source_cache_snapshot_sha256": source_cache_snapshot_sha256,
    }
    _write_json(destination, receipt)
    return fit_freeze.validate_memory_guard_receipt(
        destination,
        expected_source_manifest_content_sha256=source_manifest_content_sha256,
        expected_source_cache_snapshot_sha256=source_cache_snapshot_sha256,
    )


def _snapshot_source_caches(
    *,
    guide_cache: Path,
    feature_cache: Path,
    staging: Path,
) -> tuple[
    dict[str, object],
    dict[str, np.ndarray],
    dict[str, object],
    dict[str, dict[str, str]],
]:
    """Copy authenticated inputs into the bundle and fit only from that copy."""
    _, _, original_manifest = fit_freeze.load_authenticated_source_cache(
        guide_cache,
        feature_cache,
    )
    before = fit_freeze.source_cache_member_hashes(guide_cache, feature_cache)
    snapshot_root = staging / "source-cache"
    for kind, source in (("guide", guide_cache), ("feature", feature_cache)):
        destination = snapshot_root / kind
        destination.mkdir(parents=True, exist_ok=False)
        for filename in fit_freeze.SOURCE_CACHE_FILES[kind]:
            shutil.copyfile(source / filename, destination / filename)
    copied = fit_freeze.source_cache_member_hashes(
        snapshot_root / "guide",
        snapshot_root / "feature",
        exact_members=True,
    )
    if copied != before:
        raise ValueError("source cache changed while its frozen snapshot was copied")
    cache, derived, snapshot_manifest = fit_freeze.load_authenticated_source_cache(
        snapshot_root / "guide",
        snapshot_root / "feature",
    )
    if fit_freeze.source_manifest_content(snapshot_manifest) != (
        fit_freeze.source_manifest_content(original_manifest)
    ):
        raise ValueError("source cache snapshot content differs from authenticated input")
    return cache, derived, original_manifest, before


def _finalize_source_cache_snapshot(
    *,
    guide_cache: Path,
    feature_cache: Path,
    staging: Path,
    original_manifest: Mapping[str, object],
    before: Mapping[str, Mapping[str, str]],
) -> dict[str, object]:
    """Rehash original and copied caches after computation, then freeze receipt."""
    after = fit_freeze.source_cache_member_hashes(guide_cache, feature_cache)
    snapshot = fit_freeze.source_cache_member_hashes(
        staging / "source-cache/guide",
        staging / "source-cache/feature",
        exact_members=True,
    )
    if dict(before) != after or dict(before) != snapshot:
        raise ValueError("source cache changed during source computation")
    _, _, after_manifest = fit_freeze.load_authenticated_source_cache(
        guide_cache,
        feature_cache,
    )
    _, _, snapshot_manifest = fit_freeze.load_authenticated_source_cache(
        staging / "source-cache/guide",
        staging / "source-cache/feature",
    )
    content = fit_freeze.source_manifest_content(original_manifest)
    if (
        fit_freeze.source_manifest_content(after_manifest) != content
        or fit_freeze.source_manifest_content(snapshot_manifest) != content
    ):
        raise ValueError("source manifest changed during source computation")
    receipt = {
        "schema_version": fit_freeze.SOURCE_CACHE_SNAPSHOT_SCHEMA_VERSION,
        "experiment": fit_freeze.EXPERIMENT,
        "guide_cache": "source-cache/guide",
        "feature_cache": "source-cache/feature",
        "members": snapshot,
        "source_manifest_content_sha256": fit_freeze.sha256_json(content),
        "original_member_hashes_before": dict(before),
        "original_member_hashes_after": after,
    }
    _write_json(staging / "source-cache-snapshot.json", receipt)
    fit_freeze.validate_source_cache_snapshot(
        staging,
        expected_source_manifest=original_manifest,
    )
    return receipt


def _generate_cache(
    destination: Path,
    *,
    command_builder,
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{destination.name}.",
            suffix=".tmp",
            dir=destination.parent,
        )
    )
    try:
        command = command_builder(temporary)
        _run_command(command)
        if destination.exists():
            raise FileExistsError(f"cache appeared during generation: {destination}")
        os.replace(temporary, destination)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def generate_or_verify_caches(args: argparse.Namespace) -> None:
    guide_state = _cache_state(args.guide_cache, "guide")["state"]
    feature_state = _cache_state(args.feature_cache, "feature")["state"]
    if guide_state == "absent":
        _verify_generation_inputs(
            source_h5ad=args.source_h5ad,
            compact_k562=args.compact_k562,
            panel_csv=args.panel_csv,
            perturbmatch_perts=args.perturbmatch_perts,
            perturbmatch_genes=args.perturbmatch_genes,
            seed_table=args.seed_table,
            allow_seed_download=args.allow_seed_download,
            hash_h5ad=True,
        )
        _generate_cache(
            args.guide_cache,
            command_builder=lambda output: _extract_command(args, output=output),
        )
    elif guide_state != "complete":
        raise RuntimeError("guide cache is partial or invalid")

    if feature_state == "absent":
        _generate_cache(
            args.feature_cache,
            command_builder=lambda output: _feature_command(
                args,
                guide_cache=args.guide_cache,
                output=output,
            ),
        )
    elif feature_state != "complete":
        raise RuntimeError("feature cache is partial or invalid")

    fit_freeze.load_authenticated_source_cache(
        args.guide_cache,
        args.feature_cache,
    )


def _candidate_vector(candidate: model.Candidate) -> np.ndarray:
    return np.asarray(
        [
            candidate.response_rank,
            candidate.ridge,
            candidate.eta,
            candidate.weight_power,
            candidate.huber_c,
        ],
        dtype=np.float64,
    )


def _run_source_nulls(
    cache: Mapping[str, object],
    *,
    response_axis: model.AuthenticatedResponseAxis,
    splits: Mapping[str, Sequence[str]],
    singleton_features: np.ndarray,
    singleton_target_ids: Sequence[str],
    singleton_axis_contract: model.AxisContract,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    dict[str, object],
]:
    protocol = fit_freeze.validate_preregistration()
    within_repetitions = int(
        protocol["nulls"]["within_target_guide_feature_permutation"][
            "repetitions"
        ]
    )
    cross_repetitions = int(
        protocol["nulls"]["cross_target_guide_row_permutation"][
            "repetitions"
        ]
    )
    development_targets = tuple(
        sorted(target for name in fit_freeze.REQUIRED_SPLITS for target in splits[name])
    )
    identifiable = np.asarray(cache["identifiable"], dtype=bool)
    strata, strata_audit = fit_freeze.build_split_local_cross_target_strata(
        cache["features"],
        cache["guide_targets"],
        identifiable,
        splits,
    )
    results: dict[str, list[dict[str, object]]] = {
        "within_target": [],
        "cross_target": [],
    }
    values: dict[str, np.ndarray] = {
        "within_target": np.empty(within_repetitions, dtype=np.float64),
        "cross_target": np.empty(cross_repetitions, dtype=np.float64),
    }
    assignments = {
        "within_target": np.empty(
            (within_repetitions, len(cache["guide_targets"])),
            dtype=np.int64,
        ),
        "cross_target": np.empty(
            (cross_repetitions, len(cache["guide_targets"])),
            dtype=np.int64,
        ),
    }
    for family, repetitions in (
        ("within_target", within_repetitions),
        ("cross_target", cross_repetitions),
    ):
        for replicate in range(repetitions):
            seed = model.EXP083.stable_seed(
                model.SPLIT_SEED,
                f"exp109-{family}-null",
                replicate,
            )
            if family == "within_target":
                permutation = fit_freeze.within_target_permutation_indices(
                    cache["guide_targets"],
                    seed=seed,
                    row_mask=identifiable,
                    targets=development_targets,
                )
            else:
                permutation = fit_freeze.cross_target_permutation_indices(
                    cache["guide_targets"],
                    strata,
                    seed=seed,
                    row_mask=identifiable,
                )
            assignments[family][replicate] = permutation
            permuted = np.asarray(cache["features"])[permutation]
            null_cache = dict(cache)
            null_cache["features"] = permuted
            tune_records, selected = fit_freeze.recompute_tune_records(
                null_cache,
                response_axis=response_axis,
                splits=splits,
                singleton_features=singleton_features,
                singleton_target_ids=singleton_target_ids,
                singleton_axis_contract=singleton_axis_contract,
            )
            _, holdouts = fit_freeze.recompute_source_test_holdouts(
                null_cache,
                response_axis=response_axis,
                fit_tune_targets=tuple(
                    sorted((*splits["fit"], *splits["tune"]))
                ),
                test_targets=splits["test"],
                singleton_features=singleton_features,
                singleton_target_ids=singleton_target_ids,
                singleton_axis_contract=singleton_axis_contract,
                candidate=selected,
            )
            values[family][replicate] = float(
                holdouts["metrics"]["cosine_delta"]
            )
            selection_evidence = fit_freeze.tune_selection_payload(
                tune_records,
                selected,
            )
            selected_metrics = next(
                dict(record["metrics"])
                for record in selection_evidence["candidate_records"]
                if fit_freeze.candidate_from_payload(record["candidate"]) == selected
            )
            results[family].append(
                {
                    "replicate": replicate,
                    "seed": seed,
                    "permutation_indices_sha256": fit_freeze.sha256_array(
                        permutation
                    ),
                    "selected_candidate": fit_freeze.candidate_payload(selected),
                    "selected_tune_metrics": selected_metrics,
                    "candidate_records": selection_evidence["candidate_records"],
                    "source_test_metrics": holdouts["metrics"],
                    "permuted_features_sha256": fit_freeze.sha256_array(permuted),
                }
            )
            if (replicate + 1) % 10 == 0 or replicate + 1 == repetitions:
                print(
                    json.dumps(
                        {
                            "null_family": family,
                            "completed": replicate + 1,
                            "repetitions": repetitions,
                        }
                    ),
                    flush=True,
                )
    audit = {
        "schema_version": 2,
        "within_target": results["within_target"],
        "cross_target": results["cross_target"],
        "cross_target_strata": strata_audit,
    }
    return (
        values["within_target"],
        values["cross_target"],
        assignments["within_target"],
        assignments["cross_target"],
        audit,
    )


def _save_source_test_holdouts(
    path: Path,
    holdouts: Mapping[str, object],
    *,
    source_test_model_sha256: str,
    source_input_manifest_sha256: str,
    response_axis_sha256: str,
) -> None:
    _savez(
        path,
        schema_version=np.asarray(fit_freeze.SOURCE_TEST_HOLDOUT_SCHEMA_VERSION),
        raw_residual=np.asarray(holdouts["truth"]),
        recovered_nuisance=np.asarray(holdouts["predicted"]),
        shrink_recovered_nuisance=np.asarray(holdouts["shrink_predicted"]),
        target_labels=np.asarray(holdouts["labels"]),
        held_guide_rows=np.asarray(holdouts["held_guide_rows"], dtype=np.int64),
        source_test_model_sha256=np.asarray(source_test_model_sha256),
        source_input_manifest_sha256=np.asarray(source_input_manifest_sha256),
        response_axis_sha256=np.asarray(response_axis_sha256),
    )


def _save_singleton_corrections(
    path: Path,
    *,
    fitted: model.StableSubspaceModel,
    source_derived: Mapping[str, np.ndarray],
    source_input_manifest: Mapping[str, object],
    source_test_model_sha256: str,
    selected_model_sha256: str,
    source_input_manifest_sha256: str,
) -> dict[str, object]:
    source_axes = fit_freeze.load_canonical_source_axes()
    singleton_targets, singleton_contract = fit_freeze.canonical_axis(
        "singleton_centering_targets",
        source_axes,
    )
    panel_targets, _ = fit_freeze.canonical_axis("panel_targets", source_axes)
    genes, _ = fit_freeze.canonical_axis("gene_axis", source_axes)
    pm_targets, _ = fit_freeze.canonical_axis("perturbmatch_targets", source_axes)
    pm_genes, _ = fit_freeze.canonical_axis("perturbmatch_genes", source_axes)
    singleton_features = np.asarray(source_derived["singleton_features"])
    panel_features = np.asarray(source_derived["panel_singleton_features"])
    development_features = np.asarray(source_derived["development_features"])
    development_ids = np.asarray(source_derived["development_target_ids"])
    panel_index = np.asarray(
        source_derived["panel_singleton_source_index"],
        dtype=np.int64,
    )
    final_invariants = fit_freeze.validate_preregistered_final_invariants(
        fitted,
        source_derived=source_derived,
        source_input_manifest=source_input_manifest,
    )
    centered = model.centered_singleton_corrections(
        fitted.nuisance(singleton_features),
        singleton_targets,
        np.ones(len(singleton_targets), dtype=bool),
        eta=fitted.candidate.eta,
        centering_axis_contract=singleton_contract,
    )
    valid = panel_index >= 0
    response_corrections = np.zeros(
        (len(panel_targets), fitted.response_width),
        dtype=np.float64,
    )
    response_corrections[valid] = centered[panel_index[valid]]
    full = model.expand_response_corrections(fitted, response_corrections)
    slot_corrections = model.apply_k562_slot_retention(
        full,
        np.asarray([target in set(pm_targets) for target in panel_targets]),
        np.asarray([gene in set(pm_genes) for gene in genes]),
    )
    diagnostics = model.stable_score_extrapolation_diagnostics(
        fitted,
        development_features,
        singleton_features,
    )
    _savez(
        path,
        schema_version=np.asarray(fit_freeze.SINGLETON_CORRECTION_SCHEMA_VERSION),
        singleton_target_ids=np.asarray(singleton_targets),
        singleton_features=singleton_features,
        panel_singleton_features=panel_features,
        development_features=development_features,
        development_target_ids=development_ids,
        panel_target_ids=np.asarray(panel_targets),
        panel_singleton_source_index=panel_index,
        response_corrections=response_corrections,
        source_slot_corrections=slot_corrections,
        source_test_model_sha256=np.asarray(source_test_model_sha256),
        selected_model_sha256=np.asarray(selected_model_sha256),
        source_input_manifest_sha256=np.asarray(source_input_manifest_sha256),
        response_axis_sha256=np.asarray(fitted.response_axis_sha256),
        centering_axis_sha256=np.asarray(singleton_contract.sha256),
        extrapolation_metric_names=np.asarray(tuple(diagnostics)),
        extrapolation_metric_values=np.asarray(tuple(diagnostics.values())),
    )
    return final_invariants


def _save_held_label_nulls(
    path: Path,
    *,
    panel_support: np.ndarray,
    singleton_corrections_sha256: str,
) -> None:
    source_axes = fit_freeze.load_canonical_source_axes()
    panel_targets, _ = fit_freeze.canonical_axis("panel_targets", source_axes)
    pm_targets, _ = fit_freeze.canonical_axis("perturbmatch_targets", source_axes)
    support = np.asarray(panel_support, dtype=bool)
    strata = np.zeros(len(panel_targets), dtype=np.int64)
    strata[support] = 1
    strata[
        support
        & np.asarray([target not in set(pm_targets) for target in panel_targets])
    ] = 2
    repetitions = int(
        fit_freeze.validate_preregistration()["nulls"][
            "held_singleton_correction_label_permutation"
        ]["repetitions"]
    )
    permutations = fit_freeze.deterministic_stratified_permutations(
        strata,
        support,
        repetitions=repetitions,
        seed=fit_freeze.HELD_LABEL_NULL_SEED,
    )
    _savez(
        path,
        schema_version=np.asarray(fit_freeze.HELD_LABEL_NULL_SCHEMA_VERSION),
        seed=np.asarray(fit_freeze.HELD_LABEL_NULL_SEED),
        permutations=permutations,
        strata=strata,
        movable=support,
        panel_targets_sha256=np.asarray(fit_freeze.sha256_lines(panel_targets)),
        singleton_corrections_sha256=np.asarray(singleton_corrections_sha256),
    )


def _publish(staging: Path, output: Path) -> None:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite source output: {output}")
    os.replace(staging, output)


def run_source_only(args: argparse.Namespace) -> dict[str, object]:
    output = args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    lock = output.parent / f".{output.name}.source-only.lock"
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise RuntimeError(f"source-only runner is already active: {lock}") from exc
    os.write(descriptor, f"{os.getpid()}\n".encode())
    os.fsync(descriptor)

    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{output.name}.",
            suffix=".tmp",
            dir=output.parent,
        )
    )
    guard_state = output.parent / f".{output.name}.memory-guard"
    monitor: subprocess.Popen[bytes] | None = None
    guard_paths: dict[str, Path] | None = None
    try:
        monitor, guard_paths = _start_memory_guard(args, state_dir=guard_state)
        initial = preflight(args, hash_h5ad=False)
        if initial["errors"]:
            raise RuntimeError("; ".join(map(str, initial["errors"])))
        if args.cache_mode == "generate":
            generate_or_verify_caches(args)
        cache, source_derived, source_manifest, source_member_hashes = (
            _snapshot_source_caches(
                guide_cache=args.guide_cache,
                feature_cache=args.feature_cache,
                staging=staging,
            )
        )
        fit_freeze.write_lifecycle(staging, state="source_open")
        _write_json(staging / "source-input-manifest.json", source_manifest)
        source_input_sha256 = fit_freeze.sha256_file(
            staging / "source-input-manifest.json"
        )
        response_axis = fit_freeze.canonical_response_axis()
        singleton_targets, singleton_contract = fit_freeze.canonical_axis(
            "singleton_centering_targets"
        )
        splits = fit_freeze.development_split(
            response_axis.development_targets,
            response_axis.panel_targets,
        )
        fit_freeze.write_axis_payload(
            staging / "axes.json",
            response_axis,
            singleton_targets,
            singleton_contract,
        )

        tune_records, selected = fit_freeze.recompute_tune_records(
            cache,
            response_axis=response_axis,
            splits=splits,
            singleton_features=source_derived["singleton_features"],
            singleton_target_ids=singleton_targets,
            singleton_axis_contract=singleton_contract,
        )
        _write_json(
            staging / "tune-selection.json",
            fit_freeze.tune_selection_payload(tune_records, selected),
        )
        fit_tune_targets = tuple(sorted((*splits["fit"], *splits["tune"])))
        source_test_model, holdouts = fit_freeze.recompute_source_test_holdouts(
            cache,
            response_axis=response_axis,
            fit_tune_targets=fit_tune_targets,
            test_targets=splits["test"],
            singleton_features=source_derived["singleton_features"],
            singleton_target_ids=singleton_targets,
            singleton_axis_contract=singleton_contract,
            candidate=selected,
        )
        fit_freeze.save_frozen_model(
            staging / "source-test-model.npz",
            source_test_model,
            role="source_test_model",
        )
        source_test_model_sha256 = fit_freeze.sha256_file(
            staging / "source-test-model.npz"
        )
        _save_source_test_holdouts(
            staging / "source-test-holdouts.npz",
            holdouts,
            source_test_model_sha256=source_test_model_sha256,
            source_input_manifest_sha256=source_input_sha256,
            response_axis_sha256=source_test_model.response_axis_sha256,
        )
        validated_holdouts = fit_freeze.validate_source_test_holdouts(
            staging / "source-test-holdouts.npz",
            expected_test_targets=splits["test"],
            expected_guide_targets=cache["guide_targets"],
            expected_model_sha256=source_test_model_sha256,
            expected_source_input_sha256=source_input_sha256,
            expected_response_axis_sha256=source_test_model.response_axis_sha256,
            expected_response_width=source_test_model.response_width,
            recomputed=holdouts,
        )
        protocol = fit_freeze.validate_preregistration()
        bootstrap = fit_freeze.target_cluster_bootstrap(
            validated_holdouts["per_guide_cosine_delta"],
            validated_holdouts["labels"],
            repetitions=int(
                protocol["source_validation"][
                    "test_target_cluster_bootstrap_repetitions"
                ]
            ),
            seed=model.SPLIT_SEED,
        )
        _savez(
            staging / "source-bootstrap.npz",
            schema_version=np.asarray(1),
            seed=np.asarray(model.SPLIT_SEED),
            cosine_delta=bootstrap,
            source_test_holdouts_sha256=np.asarray(
                fit_freeze.sha256_file(staging / "source-test-holdouts.npz")
            ),
        )
        (
            within,
            cross,
            within_permutations,
            cross_permutations,
            null_audit,
        ) = _run_source_nulls(
            cache,
            response_axis=response_axis,
            splits=splits,
            singleton_features=source_derived["singleton_features"],
            singleton_target_ids=singleton_targets,
            singleton_axis_contract=singleton_contract,
        )
        shrink_delta = float(
            validated_holdouts["shrink_metrics"]["cosine_delta"]
        )
        _savez(
            staging / "source-nulls.npz",
            schema_version=np.asarray(2),
            within_target_cosine_delta=within,
            cross_target_cosine_delta=cross,
            within_target_permutation_indices=within_permutations,
            cross_target_permutation_indices=cross_permutations,
            direction_preserving_shrink_delta=np.asarray([shrink_delta]),
        )
        _write_json(staging / "source-null-selection.json", null_audit)
        source_test = fit_freeze.source_test_decision(
            validated_holdouts["metrics"],
            bootstrap_samples=bootstrap,
            within_target_null=within,
            cross_target_null=cross,
            direction_preserving_shrink_delta=shrink_delta,
        )
        source_test["selected_candidate"] = fit_freeze.candidate_payload(selected)
        if source_test["passes"] is not True:
            snapshot_receipt = _finalize_source_cache_snapshot(
                guide_cache=args.guide_cache,
                feature_cache=args.feature_cache,
                staging=staging,
                original_manifest=source_manifest,
                before=source_member_hashes,
            )
            snapshot_sha256 = fit_freeze.sha256_file(
                staging / "source-cache-snapshot.json"
            )
            assert monitor is not None and guard_paths is not None
            _finish_memory_guard(
                monitor,
                guard_paths,
                destination=staging / "memory-guard.json",
                source_manifest_content_sha256=str(
                    snapshot_receipt["source_manifest_content_sha256"]
                ),
                source_cache_snapshot_sha256=snapshot_sha256,
            )
            monitor = None
            source_test["held_evaluation_authorized"] = False
            _write_json(staging / "source-test.json", source_test)
            _write_json(
                staging / "source-run.json",
                {
                    "experiment": fit_freeze.EXPERIMENT,
                    "status": "source_rejected",
                    "decision": "REJECT",
                    "source_test": source_test,
                    "held_evaluation_enabled": False,
                    "held_outcomes_opened": False,
                },
            )
            fit_freeze.write_lifecycle(
                staging,
                state="source_rejected",
                source_input_manifest_sha256=source_input_sha256,
            )
            _publish(staging, output)
            return json.loads((output / "source-run.json").read_text())

        deployment_model = model.fit_stable_subspace_model(
            cache["guide_effects"],
            cache["features"],
            cache["guide_targets"],
            cache["identifiable"],
            response_axis.development_targets,
            response_axis,
            selected,
        )
        fit_freeze.save_frozen_model(
            staging / "selected-model.npz",
            deployment_model,
            role="deployment_model",
        )
        selected_model_sha256 = fit_freeze.sha256_file(
            staging / "selected-model.npz"
        )
        final_invariants = _save_singleton_corrections(
            staging / "singleton-corrections.npz",
            fitted=deployment_model,
            source_derived=source_derived,
            source_input_manifest=source_manifest,
            source_test_model_sha256=source_test_model_sha256,
            selected_model_sha256=selected_model_sha256,
            source_input_manifest_sha256=source_input_sha256,
        )
        panel_index = np.asarray(
            source_derived["panel_singleton_source_index"],
            dtype=np.int64,
        )
        _save_held_label_nulls(
            staging / "held-label-nulls.npz",
            panel_support=panel_index >= 0,
            singleton_corrections_sha256=fit_freeze.sha256_file(
                staging / "singleton-corrections.npz"
            ),
        )
        snapshot_receipt = _finalize_source_cache_snapshot(
            guide_cache=args.guide_cache,
            feature_cache=args.feature_cache,
            staging=staging,
            original_manifest=source_manifest,
            before=source_member_hashes,
        )
        snapshot_sha256 = fit_freeze.sha256_file(
            staging / "source-cache-snapshot.json"
        )
        assert monitor is not None and guard_paths is not None
        _finish_memory_guard(
            monitor,
            guard_paths,
            destination=staging / "memory-guard.json",
            source_manifest_content_sha256=str(
                snapshot_receipt["source_manifest_content_sha256"]
            ),
            source_cache_snapshot_sha256=snapshot_sha256,
        )
        monitor = None
        source_test["final_model_invariants"] = final_invariants
        source_test["held_evaluation_authorized"] = False
        source_test["bindings"] = {
            "source_input_manifest_sha256": source_input_sha256,
            "source_cache_snapshot_sha256": snapshot_sha256,
            "memory_guard_sha256": fit_freeze.sha256_file(
                staging / "memory-guard.json"
            ),
            "source_test_model_sha256": source_test_model_sha256,
            "source_test_holdouts_sha256": fit_freeze.sha256_file(
                staging / "source-test-holdouts.npz"
            ),
            "selected_model_sha256": selected_model_sha256,
            "tune_selection_sha256": fit_freeze.sha256_file(
                staging / "tune-selection.json"
            ),
            "source_bootstrap_sha256": fit_freeze.sha256_file(
                staging / "source-bootstrap.npz"
            ),
            "source_nulls_sha256": fit_freeze.sha256_file(
                staging / "source-nulls.npz"
            ),
            "source_null_selection_sha256": fit_freeze.sha256_file(
                staging / "source-null-selection.json"
            ),
            "axes_sha256": fit_freeze.sha256_file(staging / "axes.json"),
            "singleton_corrections_sha256": fit_freeze.sha256_file(
                staging / "singleton-corrections.npz"
            ),
            "held_label_nulls_sha256": fit_freeze.sha256_file(
                staging / "held-label-nulls.npz"
            ),
            "selected_candidate": fit_freeze.candidate_payload(selected),
            "split_hashes": {
                name: fit_freeze.sha256_lines(splits[name])
                for name in fit_freeze.REQUIRED_SPLITS
            },
        }
        _write_json(staging / "source-test.json", source_test)
        manifest = fit_freeze.build_freeze_manifest(
            staging,
            splits=splits,
            frozen_files=fit_freeze.REQUIRED_FROZEN_ARTIFACTS,
            source_test=source_test,
            response_axis=response_axis,
            centering_target_ids=singleton_targets,
            centering_axis_contract=singleton_contract,
        )
        _write_json(staging / "freeze-manifest.json", manifest)
        _write_json(
            staging / "held-evaluation-blocked.json",
            {
                "experiment": fit_freeze.EXPERIMENT,
                "blocked": True,
                "reason": (
                    "held evaluation remains a separately reviewed phase; "
                    "evaluate.py intentionally refuses authorization"
                ),
                "freeze_manifest_sha256": fit_freeze.sha256_file(
                    staging / "freeze-manifest.json"
                ),
                "held_outcomes_opened": False,
            },
        )
        _write_json(
            staging / "source-run.json",
            {
                "experiment": fit_freeze.EXPERIMENT,
                "status": "source_frozen",
                "decision": "SOURCE_ONLY_FREEZE",
                "source_test": source_test,
                "freeze_manifest_sha256": fit_freeze.sha256_file(
                    staging / "freeze-manifest.json"
                ),
                "held_evaluation_enabled": False,
                "held_outcomes_opened": False,
            },
        )
        fit_freeze.write_lifecycle(
            staging,
            state="source_frozen",
            source_input_manifest_sha256=source_input_sha256,
        )
        _publish(staging, output)
        return json.loads((output / "source-run.json").read_text())
    except BaseException:
        if monitor is not None and monitor.poll() is None:
            monitor.terminate()
            try:
                monitor.wait(timeout=5)
            except subprocess.TimeoutExpired:
                monitor.kill()
        shutil.rmtree(staging, ignore_errors=True)
        raise
    finally:
        shutil.rmtree(guard_state, ignore_errors=True)
        os.close(descriptor)
        lock.unlink(missing_ok=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument(
        "--preflight",
        action="store_true",
        help="validate paths/provenance and print the source-only execution plan",
    )
    action.add_argument(
        "--run-source",
        action="store_true",
        help="run source-only cache, tune/test/null, and atomic freeze stages",
    )
    parser.add_argument(
        "--cache-mode",
        choices=("verify", "generate"),
        default="verify",
        help="verify existing caches or generate only absent caches",
    )
    parser.add_argument("--guide-cache", type=Path, default=DEFAULT_GUIDE_CACHE)
    parser.add_argument("--feature-cache", type=Path, default=DEFAULT_FEATURE_CACHE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--source-h5ad", type=Path, default=DEFAULT_SOURCE_H5AD)
    parser.add_argument("--compact-k562", type=Path, default=DEFAULT_COMPACT_K562)
    parser.add_argument("--panel-csv", type=Path, default=DEFAULT_PANEL_CSV)
    parser.add_argument("--perturbmatch-perts", type=Path, default=DEFAULT_PM_PERTS)
    parser.add_argument("--perturbmatch-genes", type=Path, default=DEFAULT_PM_GENES)
    parser.add_argument("--seed-table", type=Path)
    parser.add_argument("--allow-seed-download", action="store_true")
    parser.add_argument("--full-hash-preflight", action="store_true")
    parser.add_argument("--max-rss-gib", type=float, default=48.0)
    parser.add_argument(
        "--guarded",
        action="store_true",
        help="required acknowledgment for the heavy source-only run",
    )
    return parser


def _resolve_paths(args: argparse.Namespace) -> None:
    for name in (
        "guide_cache",
        "feature_cache",
        "output",
        "source_h5ad",
        "compact_k562",
        "panel_csv",
        "perturbmatch_perts",
        "perturbmatch_genes",
        "seed_table",
    ):
        value = getattr(args, name)
        if value is not None:
            setattr(args, name, value.expanduser().resolve())


def main() -> int:
    args = build_parser().parse_args()
    _resolve_paths(args)
    if args.max_rss_gib <= 0:
        raise SystemExit("--max-rss-gib must be positive")
    if args.preflight:
        result = preflight(args, hash_h5ad=args.full_hash_preflight)
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0 if not result["errors"] else 2
    if not args.guarded:
        raise SystemExit("refusing heavy source-only execution without --guarded")
    started = time.time()
    result = run_source_only(args)
    result["elapsed_seconds"] = time.time() - started
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["decision"] == "SOURCE_ONLY_FREEZE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
