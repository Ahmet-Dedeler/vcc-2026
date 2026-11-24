#!/usr/bin/env python3
"""Fit and freeze Experiment 107 conditional contrasts before held truth is opened."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

import numpy as np
from scipy import sparse

import model_input
import modeling
import source_artifacts


HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
DEFAULT_INPUT = WORK / "source-model-input.npz"
DEFAULT_OUTPUT = WORK / "frozen-source-model"
MAX_RSS_GIB = 24.0
MIN_FREE_PERCENT = 5
MEMORY_POLL_SECONDS = 0.25
EXPECTED_TARGET_COUNT = 18_465
EXPECTED_GUIDE_COUNT = 83_401
EXPECTED_ORIGINAL_RESPONSE_COUNT = 374
EXPECTED_EXPORTABLE_RESPONSE_COUNT = 349


class SourceFitError(RuntimeError):
    """The frozen source-fit contract was violated."""


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_bytes), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise SourceFitError(f"cannot read {path}: {error}") from error
    if not isinstance(value, dict):
        raise SourceFitError(f"{path} must contain a JSON object")
    return value


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def guard_memory(stage: str, limit_gib: float) -> float:
    if not np.isfinite(limit_gib) or not 0 < limit_gib <= MAX_RSS_GIB:
        raise SourceFitError(f"memory limit must lie in (0, {MAX_RSS_GIB}] GiB")
    peak = peak_rss_gib()
    if peak >= limit_gib:
        raise SourceFitError(
            f"{stage}: peak RSS {peak:.3f} GiB reached limit {limit_gib:.3f} GiB"
        )
    return peak


@contextmanager
def active_process_tree_memory_monitor(
    *,
    limit_gib: float,
    scratch_parent: Path,
) -> Iterator[dict[str, Any]]:
    """Run a separate watchdog that can terminate this fit during excess use."""

    guard_memory("memory-monitor startup", limit_gib)
    monitor_script = HERE.parents[1] / "scripts" / "run_with_memory_guard.py"
    if not monitor_script.is_file():
        raise SourceFitError(f"memory monitor is missing: {monitor_script}")
    scratch_parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".source-fit-memory-monitor.",
        dir=scratch_parent,
    ) as directory:
        monitor_directory = Path(directory)
        ready_path = monitor_directory / "ready.json"
        stop_path = monitor_directory / "stop"
        report_path = monitor_directory / "report.json"
        command = [
            sys.executable,
            str(monitor_script),
            "--monitor-pid",
            str(os.getpid()),
            "--max-rss-gib",
            str(limit_gib),
            "--min-free-percent",
            str(MIN_FREE_PERCENT),
            "--poll-seconds",
            str(MEMORY_POLL_SECONDS),
            "--ready-file",
            str(ready_path),
            "--stop-file",
            str(stop_path),
            "--report-file",
            str(report_path),
            "--quiet",
        ]
        monitor = subprocess.Popen(
            command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            start_new_session=True,
        )
        deadline = time.monotonic() + 10.0
        while not ready_path.exists():
            if monitor.poll() is not None:
                raise SourceFitError("active memory monitor exited before readiness")
            if time.monotonic() >= deadline:
                monitor.terminate()
                raise SourceFitError("active memory monitor did not become ready")
            time.sleep(0.02)
        ready = load_json(ready_path)
        if (
            ready.get("status") != "active"
            or ready.get("root_pid") != os.getpid()
            or ready.get("max_rss_gib") != limit_gib
        ):
            monitor.terminate()
            raise SourceFitError("active memory monitor readiness receipt is invalid")

        state: dict[str, Any] = {"ready": ready}
        failed = False
        try:
            yield state
        except BaseException:
            failed = True
            raise
        finally:
            stop_path.touch()
            try:
                return_code = monitor.wait(timeout=10)
            except subprocess.TimeoutExpired:
                monitor.terminate()
                try:
                    monitor.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    monitor.kill()
                    monitor.wait(timeout=5)
                return_code = monitor.returncode
            if report_path.exists():
                state["report"] = load_json(report_path)
            if not failed:
                report = state.get("report")
                if return_code != 0 or not isinstance(report, dict):
                    raise SourceFitError(
                        "active memory monitor did not produce a passing receipt"
                    )
                if report.get("status") != "complete":
                    raise SourceFitError(
                        "active memory monitor did not complete cleanly"
                    )


def atomic_savez(path: Path, **arrays: np.ndarray) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            np.savez_compressed(handle, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def _binary_csr(
    indptr: np.ndarray,
    indices: np.ndarray,
    *,
    rows: int,
    columns: int,
    name: str,
) -> sparse.csr_matrix:
    try:
        return model_input.csr_from_components(
            indptr,
            indices,
            rows=rows,
            columns=columns,
        )
    except model_input.ModelInputError as error:
        raise SourceFitError(f"{name}: {error}") from error


def load_frozen_input(
    input_path: Path,
    *,
    source_verification: dict[str, Any] | None = None,
) -> tuple[
    np.ndarray,
    sparse.csr_matrix,
    np.ndarray,
    np.ndarray,
    modeling.GuideEvidence,
    dict[str, Any],
    dict[str, Any],
    Any,
]:
    """Load and revalidate the source-only input without held-data access."""

    input_path = input_path.resolve()
    report_path = input_path.with_suffix(".json")
    report = load_json(report_path)
    if report.get("status") != "complete":
        raise SourceFitError("source-model input is incomplete")
    if report.get("primary_condition") != "untreated":
        raise SourceFitError(
            "source-model input is not the frozen untreated population"
        )
    if report.get("nadig_or_challenge_outcome_opened") is not False:
        raise SourceFitError(
            "source-model input does not certify the held-data boundary"
        )
    if report.get("implementation_sha256") != sha256_file(HERE / "model_input.py"):
        raise SourceFitError("source-model input implementation changed")
    if report.get("output_path") != str(input_path.relative_to(HERE)):
        raise SourceFitError("source-model input path changed")
    if report.get("output_sha256") != sha256_file(input_path):
        raise SourceFitError("source-model input payload changed")
    verified_source = (
        source_artifacts.verify_source_artifacts()
        if source_verification is None
        else source_verification
    )
    try:
        model_input.validate_model_input_lineage(
            report,
            source_verification=verified_source,
        )
    except model_input.ModelInputError as error:
        raise SourceFitError(str(error)) from error

    with np.load(input_path, allow_pickle=False) as values:
        required = {
            "response",
            "nuisance",
            "channels",
            "barcodes",
            "target_indptr",
            "target_indices",
            "guide_indptr",
            "guide_indices",
            "target_symbols",
            "guide_ids",
            "guide_to_target",
            "response_symbols",
            "response_vcc_indices",
            "response_exportable",
            "original_response_symbols",
            "original_response_vcc_indices",
            "original_response_exportable",
            "response_original_indices",
            "precomputed_training_channels",
            "precomputed_distinct_called_guides",
            "precomputed_guides_with_three_untreated_cells",
            "precomputed_untreated_carrier_cells",
            "precomputed_training_channel_count",
            "precomputed_eligible",
        }
        if set(values.files) != required:
            raise SourceFitError("source-model input keys changed")
        arrays = {name: values[name] for name in required}

    response = np.asarray(arrays["response"], dtype=np.float64)
    nuisance = np.asarray(arrays["nuisance"], dtype=np.float64)
    channels = np.asarray(arrays["channels"])
    target_symbols = np.asarray(arrays["target_symbols"]).astype(str)
    guide_ids = np.asarray(arrays["guide_ids"]).astype(str)
    if (
        response.ndim != 2
        or response.shape[1] != EXPECTED_EXPORTABLE_RESPONSE_COUNT
        or nuisance.shape != (response.shape[0], 3)
    ):
        raise SourceFitError("source response and nuisance axes changed")
    if channels.shape != (response.shape[0],):
        raise SourceFitError("source channel axis changed")
    if len(set(channels.tolist())) != 16 or set(channels.tolist()) != set(
        modeling.FROZEN_CHANNELS
    ):
        raise SourceFitError("source input does not contain exactly channels 1..16")
    if target_symbols.shape != (EXPECTED_TARGET_COUNT,) or guide_ids.shape != (
        EXPECTED_GUIDE_COUNT,
    ):
        raise SourceFitError("source target or guide axis changed")

    target_design = _binary_csr(
        arrays["target_indptr"],
        arrays["target_indices"],
        rows=response.shape[0],
        columns=target_symbols.size,
        name="target design",
    )
    guide_presence = _binary_csr(
        arrays["guide_indptr"],
        arrays["guide_indices"],
        rows=response.shape[0],
        columns=guide_ids.size,
        name="guide presence",
    )
    guide_to_target = np.asarray(arrays["guide_to_target"], dtype=np.int64)
    evidence = modeling.GuideEvidence(
        guide_presence=guide_presence,
        guide_target_indices=guide_to_target,
    )
    precomputed_support = modeling.TargetSupportEvidence(
        training_channels=np.asarray(
            arrays["precomputed_training_channels"],
            dtype=np.int64,
        ),
        distinct_called_guides=np.asarray(
            arrays["precomputed_distinct_called_guides"],
            dtype=np.int64,
        ),
        guides_with_three_untreated_cells=np.asarray(
            arrays["precomputed_guides_with_three_untreated_cells"],
            dtype=np.int64,
        ),
        untreated_carrier_cells=np.asarray(
            arrays["precomputed_untreated_carrier_cells"],
            dtype=np.int64,
        ),
        training_channel_count=np.asarray(
            arrays["precomputed_training_channel_count"],
            dtype=np.int64,
        ),
        supported=np.asarray(arrays["precomputed_eligible"], dtype=bool),
    )
    observed_support = modeling.derive_target_support(
        evidence,
        target_design,
        channels,
        training_channels=modeling.FROZEN_CHANNELS,
    )
    try:
        model_input.assert_target_support_equal(
            precomputed_support,
            observed_support,
            label="frozen untreated target support",
        )
    except model_input.ModelInputError as error:
        raise SourceFitError(str(error)) from error
    response_symbols = np.asarray(arrays["response_symbols"]).astype(str)
    response_vcc_indices = np.asarray(
        arrays["response_vcc_indices"],
        dtype=np.int32,
    )
    response_exportable = np.asarray(
        arrays["response_exportable"],
        dtype=bool,
    )
    original_response_symbols = np.asarray(arrays["original_response_symbols"]).astype(
        str
    )
    original_response_vcc_indices = np.asarray(
        arrays["original_response_vcc_indices"],
        dtype=np.int32,
    )
    original_response_exportable = np.asarray(
        arrays["original_response_exportable"],
        dtype=bool,
    )
    response_original_indices = np.asarray(
        arrays["response_original_indices"],
        dtype=np.int64,
    )
    axes = {
        "target_symbols": target_symbols,
        "guide_ids": guide_ids,
        "guide_to_target": guide_to_target,
        "response_symbols": response_symbols,
        "response_vcc_indices": response_vcc_indices,
        "response_exportable": response_exportable,
        "original_response_symbols": original_response_symbols,
        "original_response_vcc_indices": original_response_vcc_indices,
        "original_response_exportable": original_response_exportable,
        "response_original_indices": response_original_indices,
        "barcodes": np.asarray(arrays["barcodes"]).astype(str),
        "precomputed_support": precomputed_support,
    }
    if (
        response_symbols.shape != (EXPECTED_EXPORTABLE_RESPONSE_COUNT,)
        or response_vcc_indices.shape != (EXPECTED_EXPORTABLE_RESPONSE_COUNT,)
        or response_exportable.shape != (EXPECTED_EXPORTABLE_RESPONSE_COUNT,)
        or not np.all(response_exportable)
        or np.any(response_vcc_indices < 0)
        or original_response_symbols.shape != (EXPECTED_ORIGINAL_RESPONSE_COUNT,)
        or original_response_vcc_indices.shape != (EXPECTED_ORIGINAL_RESPONSE_COUNT,)
        or original_response_exportable.shape != (EXPECTED_ORIGINAL_RESPONSE_COUNT,)
        or np.count_nonzero(original_response_exportable)
        != EXPECTED_EXPORTABLE_RESPONSE_COUNT
        or response_original_indices.shape != (EXPECTED_EXPORTABLE_RESPONSE_COUNT,)
        or not np.array_equal(
            response_original_indices,
            np.flatnonzero(original_response_exportable),
        )
        or not np.array_equal(
            response_symbols,
            original_response_symbols[response_original_indices],
        )
        or not np.array_equal(
            response_vcc_indices,
            original_response_vcc_indices[response_original_indices],
        )
        or axes["barcodes"].shape != (response.shape[0],)
    ):
        raise SourceFitError("source response or cell metadata axis changed")
    if (
        report.get("response_genes") != EXPECTED_EXPORTABLE_RESPONSE_COUNT
        or report.get("modeled_response_genes") != EXPECTED_EXPORTABLE_RESPONSE_COUNT
        or report.get("exportable_response_genes") != EXPECTED_EXPORTABLE_RESPONSE_COUNT
        or report.get("original_response_genes") != EXPECTED_ORIGINAL_RESPONSE_COUNT
        or report.get("original_vcc_mapped_response_genes")
        != model_input.EXPECTED_VCC_MAPPED_RESPONSE_COUNT
        or report.get("excluded_response_genes")
        != EXPECTED_ORIGINAL_RESPONSE_COUNT - EXPECTED_EXPORTABLE_RESPONSE_COUNT
        or report.get("original_response_axis_sha256")
        != model_input.sha256_arrays(
            symbols=original_response_symbols,
            vcc_indices=original_response_vcc_indices,
            exportable=original_response_exportable,
        )
        or report.get("modeled_response_axis_sha256")
        != model_input.sha256_arrays(
            original_indices=response_original_indices.astype(np.int32),
            symbols=response_symbols,
            vcc_indices=response_vcc_indices,
        )
    ):
        raise SourceFitError("source response provenance changed")
    if not np.all(np.isfinite(response)) or not np.all(np.isfinite(nuisance)):
        raise SourceFitError("source-model input contains non-finite values")
    design_audit_path, _ = model_input.design_audit_paths(input_path)
    try:
        raw_design_audit = model_input.load_production_design_audit(
            design_audit_path,
            target_design=target_design,
            nuisance=nuisance,
            channels=channels,
            guide_evidence=evidence,
        )
        validated_design_audit = modeling.validate_production_design_audit(
            raw_design_audit,
            target_design,
            nuisance,
            channels,
            evidence,
        )
    except (model_input.ModelInputError, modeling.ModelingError) as error:
        raise SourceFitError(str(error)) from error
    return (
        response,
        target_design,
        nuisance,
        channels,
        evidence,
        axes,
        report,
        validated_design_audit,
    )


def _solver_report(
    diagnostics: modeling.RidgeSolveDiagnostics | None,
) -> dict[str, Any]:
    if diagnostics is None:
        raise SourceFitError("production estimator lacks ridge diagnostics")
    return {
        "ridge_lambda": diagnostics.ridge_lambda,
        "iterations": diagnostics.iterations,
        "operator_applications": diagnostics.operator_applications,
        "warm_started": diagnostics.warm_started,
        "all_components_converged": bool(np.all(diagnostics.converged)),
        "maximum_relative_residual": float(
            np.max(diagnostics.relative_residual_norms, initial=0.0)
        ),
        "relative_tolerance": diagnostics.relative_tolerance,
        "absolute_tolerance": diagnostics.absolute_tolerance,
    }


def _pack_fold_vectors(
    values: tuple[np.ndarray, ...],
) -> tuple[np.ndarray, np.ndarray]:
    indptr = np.zeros(len(values) + 1, dtype=np.int64)
    for index, value in enumerate(values, start=1):
        indptr[index] = indptr[index - 1] + np.asarray(value).size
    data = (
        np.concatenate([np.asarray(value, dtype=np.float64) for value in values])
        if indptr[-1]
        else np.zeros(0, dtype=np.float64)
    )
    return indptr, data


def _fold_stability(
    fold_effects: np.ndarray,
    fold_transferable: np.ndarray,
    response_exportable: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return pairwise fold agreement for identifiable conditional contrasts."""

    effects = np.asarray(fold_effects, dtype=np.float64)
    transferable = np.asarray(fold_transferable, dtype=bool)
    exportable = np.asarray(response_exportable, dtype=bool)
    if effects.ndim != 3 or transferable.shape != effects.shape[:2]:
        raise SourceFitError("fold stability axes changed")
    if exportable.shape != (effects.shape[2],) or not np.any(exportable):
        raise SourceFitError("fold stability response mask changed")
    pairs = np.asarray(
        [
            (first, second)
            for first in range(effects.shape[0])
            for second in range(first + 1, effects.shape[0])
        ],
        dtype=np.int8,
    )
    cosine = np.full((pairs.shape[0], effects.shape[1]), np.nan, dtype=np.float64)
    sign_agreement = np.full_like(cosine, np.nan)
    for pair_index, (first, second) in enumerate(pairs):
        left = effects[first][:, exportable]
        right = effects[second][:, exportable]
        left_norm = np.linalg.norm(left, axis=1)
        right_norm = np.linalg.norm(right, axis=1)
        valid = (
            transferable[first]
            & transferable[second]
            & np.isfinite(left_norm)
            & np.isfinite(right_norm)
            & (left_norm > 0.0)
            & (right_norm > 0.0)
        )
        cosine[pair_index, valid] = np.sum(left[valid] * right[valid], axis=1) / (
            left_norm[valid] * right_norm[valid]
        )
        sign_agreement[pair_index, valid] = np.mean(
            np.signbit(left[valid]) == np.signbit(right[valid]),
            axis=1,
        )
    return pairs, cosine, sign_agreement


def _restore_original_response_axis(
    values: np.ndarray,
    *,
    original_response_count: int,
    response_original_indices: np.ndarray,
) -> np.ndarray:
    """Zero-pad compact modeled responses onto the original source axis."""

    array = np.asarray(values)
    indices = np.asarray(response_original_indices, dtype=np.int64)
    if array.ndim < 1 or array.shape[-1] != indices.size:
        raise SourceFitError("compact response axis changed before zero padding")
    if (
        indices.shape != (EXPECTED_EXPORTABLE_RESPONSE_COUNT,)
        or np.any(indices < 0)
        or np.any(indices >= original_response_count)
        or np.any(np.diff(indices) <= 0)
    ):
        raise SourceFitError("original response index mapping changed")
    expanded = np.zeros(
        (*array.shape[:-1], original_response_count),
        dtype=array.dtype,
    )
    expanded[..., indices] = array
    return expanded


def fit_and_freeze(
    input_path: Path = DEFAULT_INPUT,
    output_directory: Path = DEFAULT_OUTPUT,
    *,
    memory_limit_gib: float = MAX_RSS_GIB,
) -> dict[str, Any]:
    """Fit and publish minimum-norm high-MOI conditional target contrasts."""

    input_path = input_path.resolve()
    if HERE not in input_path.parents:
        raise SourceFitError(
            f"refusing to read source-model input outside experiment: {input_path}"
        )
    output_directory = output_directory.resolve()
    if HERE not in output_directory.parents:
        raise SourceFitError(
            f"refusing to write outside the experiment directory: {output_directory}"
        )
    if output_directory.exists():
        raise SourceFitError(f"refusing to overwrite {output_directory}")
    source_verification = source_artifacts.verify_source_artifacts()
    (
        response,
        targets,
        nuisance,
        channels,
        evidence,
        axes,
        input_report,
        validated_design_audit,
    ) = load_frozen_input(
        input_path,
        source_verification=source_verification,
    )
    if input_report.get("source_raw_archive_sha256") != source_verification.get(
        "raw_archive_sha256"
    ):
        raise SourceFitError("source-model input used a different source archive")
    if input_report.get("source_oracle_receipt_sha256") != source_verification.get(
        "oracle_receipt_sha256"
    ):
        raise SourceFitError("source-model input used a different R-oracle receipt")
    for field in ("preflight_receipt_sha256", "numpy_version", "scipy_version"):
        input_field = (
            "source_preflight_receipt_sha256"
            if field == "preflight_receipt_sha256"
            else field
        )
        if input_report.get(input_field) != source_verification.get(field):
            raise SourceFitError(f"source-model input lineage changed at {input_field}")

    output_directory.parent.mkdir(parents=True, exist_ok=True)
    monitor_state: dict[str, Any]
    with active_process_tree_memory_monitor(
        limit_gib=memory_limit_gib,
        scratch_parent=output_directory.parent,
    ) as monitor_state:
        result = modeling.select_production_estimator(
            response,
            targets,
            nuisance,
            channels,
            evidence,
            design_audit=validated_design_audit,
        )
    memory_report = monitor_state["report"]
    estimator = result.estimator
    try:
        model_input.assert_target_support_equal(
            axes["precomputed_support"],
            result.final_support,
            label="fitted all-channel untreated target support",
        )
    except model_input.ModelInputError as error:
        raise SourceFitError(str(error)) from error
    if not np.any(estimator.diagnostics.transferable):
        raise SourceFitError("source estimator exported no transferable targets")
    if np.any(estimator.exported_effects[~estimator.diagnostics.transferable] != 0.0):
        raise SourceFitError("unsupported source targets are not exact zero fallback")
    if not np.all(np.isfinite(estimator.exported_effects)):
        raise SourceFitError("source estimator exported non-finite effects")
    for name, error in {
        "nullspace annihilation": (estimator.diagnostics.nullspace_annihilation_error),
        "projection idempotence": (estimator.diagnostics.projection_idempotence_error),
        "gauge-shift invariance": (estimator.diagnostics.gauge_shift_invariance_error),
    }.items():
        if error > modeling.PROJECTION_POSTCONDITION_TOLERANCE:
            raise SourceFitError(f"source estimator failed {name}: {error:.3e}")
    fold_geometry = result.fold_geometry_diagnostics
    if len(fold_geometry) != len(modeling.FROZEN_HELD_CHANNEL_FOLDS):
        raise SourceFitError("source selector did not preserve every fold diagnostic")
    fold_zero_variance = np.stack(
        [diagnostics.zero_variance for diagnostics in fold_geometry]
    )
    fold_nuisance_aliased = np.stack(
        [diagnostics.nuisance_aliased for diagnostics in fold_geometry]
    )
    fold_projection_diagonal = np.stack(
        [diagnostics.projection_diagonal for diagnostics in fold_geometry]
    )
    fold_transferable = np.stack(
        [diagnostics.transferable for diagnostics in fold_geometry]
    )
    fold_support_training_channels = np.stack(
        [support.training_channels for support in result.fold_support]
    )
    fold_support_distinct_called_guides = np.stack(
        [support.distinct_called_guides for support in result.fold_support]
    )
    fold_support_guides_with_three = np.stack(
        [support.guides_with_three_untreated_cells for support in result.fold_support]
    )
    fold_support_untreated_carriers = np.stack(
        [support.untreated_carrier_cells for support in result.fold_support]
    )
    fold_support_training_channel_count = np.stack(
        [support.training_channel_count for support in result.fold_support]
    )
    fold_support_supported = np.stack(
        [support.supported for support in result.fold_support]
    )
    unmodeled_indptr, unmodeled_values = _pack_fold_vectors(
        tuple(diagnostics.held_unmodeled_target_burden for diagnostics in fold_geometry)
    )
    fold_pairs, fold_effect_cosine, fold_effect_sign_agreement = _fold_stability(
        result.fold_selected_effects,
        fold_transferable,
        axes["response_exportable"],
    )
    original_response_count = axes["original_response_symbols"].size
    response_original_indices = axes["response_original_indices"]
    expanded_exported_effects = _restore_original_response_axis(
        estimator.exported_effects,
        original_response_count=original_response_count,
        response_original_indices=response_original_indices,
    )
    expanded_coefficients = _restore_original_response_axis(
        estimator.coefficients_binary_units,
        original_response_count=original_response_count,
        response_original_indices=response_original_indices,
    )
    expanded_projected = _restore_original_response_axis(
        estimator.projected_coefficients,
        original_response_count=original_response_count,
        response_original_indices=response_original_indices,
    )
    expanded_response_slopes = _restore_original_response_axis(
        estimator.nuisance_response_slopes,
        original_response_count=original_response_count,
        response_original_indices=response_original_indices,
    )
    expanded_response_basis = _restore_original_response_axis(
        estimator.response_basis,
        original_response_count=original_response_count,
        response_original_indices=response_original_indices,
    )
    expanded_fold_effects = _restore_original_response_axis(
        result.fold_selected_effects,
        original_response_count=original_response_count,
        response_original_indices=response_original_indices,
    )
    forbidden_response = ~axes["original_response_exportable"]
    for name, values in {
        "exported effects": expanded_exported_effects,
        "binary-unit coefficients": expanded_coefficients,
        "projected coefficients": expanded_projected,
        "nuisance response slopes": expanded_response_slopes,
        "response basis": expanded_response_basis,
        "fold selected effects": expanded_fold_effects,
    }.items():
        if np.any(values[..., forbidden_response] != 0.0):
            raise SourceFitError(
                f"{name} contain a nonzero forbidden response coordinate"
            )

    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{output_directory.name}.",
            dir=output_directory.parent,
        )
    )
    try:
        model_path = staging / "source-model.npz"
        atomic_savez(
            model_path,
            exported_effects=expanded_exported_effects,
            coefficients_binary_units=expanded_coefficients,
            projected_coefficients=expanded_projected,
            nuisance_response_slopes=expanded_response_slopes,
            nuisance_target_slopes=estimator.nuisance_target_slopes,
            response_basis=expanded_response_basis,
            nullspace_basis=estimator.nullspace_basis,
            evidence_support=estimator.diagnostics.evidence_support,
            zero_variance=estimator.diagnostics.zero_variance,
            nuisance_aliased=estimator.diagnostics.nuisance_aliased,
            fitted=estimator.diagnostics.fitted,
            projection_diagonal=estimator.diagnostics.projection_diagonal,
            transferable=estimator.diagnostics.transferable,
            target_scales=estimator.diagnostics.target_scales,
            singular_values=estimator.diagnostics.singular_values,
            retained_singular_directions=(
                estimator.diagnostics.retained_singular_directions
            ),
            target_symbols=axes["target_symbols"],
            response_symbols=axes["original_response_symbols"],
            response_vcc_indices=axes["original_response_vcc_indices"],
            response_exportable=axes["original_response_exportable"],
            selected_lambda=np.asarray(result.selected_lambda),
            lambda_grid=result.lambda_grid,
            held_cell_counts=result.held_cell_counts,
            fold_mean_gains=result.fold_mean_gains,
            mean_held_cell_gains=result.mean_held_cell_gains,
            fold_selected_effects=expanded_fold_effects,
            fold_zero_variance=fold_zero_variance,
            fold_nuisance_aliased=fold_nuisance_aliased,
            fold_projection_diagonal=fold_projection_diagonal,
            fold_transferable=fold_transferable,
            fold_support_training_channels=fold_support_training_channels,
            fold_support_distinct_called_guides=(fold_support_distinct_called_guides),
            fold_support_guides_with_three_untreated_cells=(
                fold_support_guides_with_three
            ),
            fold_support_untreated_carrier_cells=(fold_support_untreated_carriers),
            fold_support_training_channel_count=(fold_support_training_channel_count),
            fold_support_supported=fold_support_supported,
            final_support_training_channels=result.final_support.training_channels,
            final_support_distinct_called_guides=(
                result.final_support.distinct_called_guides
            ),
            final_support_guides_with_three_untreated_cells=(
                result.final_support.guides_with_three_untreated_cells
            ),
            final_support_untreated_carrier_cells=(
                result.final_support.untreated_carrier_cells
            ),
            final_support_training_channel_count=(
                result.final_support.training_channel_count
            ),
            final_support_supported=result.final_support.supported,
            fold_unmodeled_indptr=unmodeled_indptr,
            fold_unmodeled_target_burden=unmodeled_values,
            fold_stability_pairs=fold_pairs,
            fold_effect_cosine=fold_effect_cosine,
            fold_effect_sign_agreement=fold_effect_sign_agreement,
            projection_nullspace_annihilation_error=np.asarray(
                estimator.diagnostics.nullspace_annihilation_error
            ),
            projection_idempotence_error=np.asarray(
                estimator.diagnostics.projection_idempotence_error
            ),
            projection_gauge_shift_invariance_error=np.asarray(
                estimator.diagnostics.gauge_shift_invariance_error
            ),
        )
        with np.load(model_path, allow_pickle=False) as frozen:
            frozen_projection_errors = modeling.validate_production_projection(
                frozen["coefficients_binary_units"],
                frozen["projected_coefficients"],
                frozen["nullspace_basis"],
            )
        fold_support = [
            {
                "held_channels": list(result.held_channel_folds[index]),
                "training_channels": support.training_channels.tolist(),
                "supported_targets": int(np.count_nonzero(support.supported)),
                "distinct_guide_gate": int(
                    np.count_nonzero(
                        support.distinct_called_guides
                        >= modeling.MIN_DISTINCT_CALLED_GUIDES
                    )
                ),
                "three_cell_guide_gate": int(
                    np.count_nonzero(
                        support.guides_with_three_untreated_cells
                        >= modeling.MIN_GUIDES_WITH_THREE_UNTREATED_CELLS
                    )
                ),
                "carrier_gate": int(
                    np.count_nonzero(
                        support.untreated_carrier_cells
                        >= modeling.MIN_UNTREATED_CARRIER_CELLS
                    )
                ),
                "channel_gate": int(
                    np.count_nonzero(
                        support.training_channel_count >= modeling.MIN_TRAINING_CHANNELS
                    )
                ),
                "solver": [
                    _solver_report(item)
                    for item in result.fold_solver_diagnostics[index]
                ],
                "zero_variance_targets": int(
                    np.count_nonzero(fold_geometry[index].zero_variance)
                ),
                "nuisance_aliased_targets": int(
                    np.count_nonzero(fold_geometry[index].nuisance_aliased)
                ),
                "transferable_targets": int(
                    np.count_nonzero(fold_geometry[index].transferable)
                ),
                "effective_rank": fold_geometry[index].effective_rank,
                "condition_number": fold_geometry[index].condition_number,
                "held_unmodeled_target_burden": {
                    "nonzero_cells": int(
                        np.count_nonzero(
                            fold_geometry[index].held_unmodeled_target_burden
                        )
                    ),
                    "mean": float(
                        np.mean(fold_geometry[index].held_unmodeled_target_burden)
                    ),
                    "maximum": float(
                        np.max(
                            fold_geometry[index].held_unmodeled_target_burden,
                            initial=0.0,
                        )
                    ),
                },
            }
            for index, support in enumerate(result.fold_support)
        ]
        final_solver = _solver_report(estimator.ridge_solve_diagnostics)
        manifest = {
            "status": "frozen_before_held_outcome_access",
            "experiment": "107-jurkat-high-moi",
            "primary_condition": "untreated",
            "held_truth_opened": False,
            "nadig_expression_or_scores_opened": False,
            "challenge_outcomes_opened": False,
            "effect_semantics": (
                "identifiable_high_moi_conditional_target_contrast_"
                "minimum_norm_row_space_component_not_singleton_marginal_effect"
            ),
            "input_path": str(input_path.relative_to(HERE)),
            "input_sha256": sha256_file(input_path),
            "input_report_sha256": sha256_file(input_path.with_suffix(".json")),
            "model_input_implementation_sha256": input_report["implementation_sha256"],
            "modeling_implementation_sha256": sha256_file(HERE / "modeling.py"),
            "fit_implementation_sha256": sha256_file(Path(__file__)),
            "source_artifacts_implementation_sha256": sha256_file(
                HERE / "source_artifacts.py"
            ),
            "source_raw_archive_sha256": source_verification["raw_archive_sha256"],
            "source_oracle_receipt_sha256": source_verification[
                "oracle_receipt_sha256"
            ],
            "source_preflight_receipt_sha256": source_verification[
                "preflight_receipt_sha256"
            ],
            "design_audit_path": input_report["design_audit_path"],
            "design_audit_report_sha256": input_report["design_audit_report_sha256"],
            "design_audit_output_sha256": input_report["design_audit_output_sha256"],
            "design_audit_validated_before_response_fit": True,
            "numpy_version": source_verification["numpy_version"],
            "scipy_version": source_verification["scipy_version"],
            "cells": int(response.shape[0]),
            "targets": int(targets.shape[1]),
            "response_genes": int(original_response_count),
            "modeled_response_genes": int(response.shape[1]),
            "excluded_response_genes": int(original_response_count - response.shape[1]),
            "target_design_nnz": int(targets.nnz),
            "guide_design_nnz": int(evidence.guide_presence.nnz),
            "selected_lambda": result.selected_lambda,
            "lambda_grid": result.lambda_grid.tolist(),
            "held_channel_folds": [list(fold) for fold in result.held_channel_folds],
            "fold_mean_gains": result.fold_mean_gains.tolist(),
            "mean_held_cell_gains": result.mean_held_cell_gains.tolist(),
            "held_cell_counts": result.held_cell_counts.tolist(),
            "folds": fold_support,
            "channel_stability": {
                "definition": (
                    "pairwise selected-lambda source-effect cosine and sign "
                    "agreement computed only on the compact 349-coordinate "
                    "modeled response axis"
                ),
                "pair_count": int(fold_pairs.shape[0]),
                "targets_with_any_pair": int(
                    np.count_nonzero(np.any(np.isfinite(fold_effect_cosine), axis=0))
                ),
                "median_pairwise_cosine": (
                    float(np.nanmedian(fold_effect_cosine))
                    if np.any(np.isfinite(fold_effect_cosine))
                    else None
                ),
                "median_pairwise_sign_agreement": (
                    float(np.nanmedian(fold_effect_sign_agreement))
                    if np.any(np.isfinite(fold_effect_sign_agreement))
                    else None
                ),
                "full_arrays_in_model_file": True,
            },
            "projection_postconditions": {
                "tolerance": modeling.PROJECTION_POSTCONDITION_TOLERANCE,
                "nullspace_annihilation_error": frozen_projection_errors[0],
                "projection_idempotence_error": frozen_projection_errors[1],
                "gauge_shift_invariance_error": frozen_projection_errors[2],
                "validated_after_serialization": True,
            },
            "final_support": {
                "evidence_supported_targets": int(
                    np.count_nonzero(result.final_support.supported)
                ),
                "fitted_targets": int(np.count_nonzero(estimator.diagnostics.fitted)),
                "transferable_targets": int(
                    np.count_nonzero(estimator.diagnostics.transferable)
                ),
                "zero_variance_targets": int(
                    np.count_nonzero(estimator.diagnostics.zero_variance)
                ),
                "nuisance_aliased_targets": int(
                    np.count_nonzero(estimator.diagnostics.nuisance_aliased)
                ),
                "effective_rank": estimator.diagnostics.effective_rank,
                "condition_number": estimator.diagnostics.condition_number,
                "minimum_transferable_projection_diagonal": float(
                    np.min(
                        estimator.diagnostics.projection_diagonal[
                            estimator.diagnostics.transferable
                        ]
                    )
                ),
                "full_per_target_fold_and_final_arrays_in_model_file": True,
            },
            "final_solver": final_solver,
            "response_basis_rank": estimator.effective_response_rank,
            "exportable_response_genes": int(
                np.count_nonzero(axes["original_response_exportable"])
            ),
            "active_memory_monitor": memory_report,
            "model_file": "source-model.npz",
            "model_sha256": sha256_file(model_path),
        }
        atomic_write_json(staging / "manifest.json", manifest)
        os.replace(staging, output_directory)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    return manifest


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    result.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT)
    result.add_argument("--memory-limit-gib", type=float, default=MAX_RSS_GIB)
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        report = fit_and_freeze(
            args.input,
            args.output_directory,
            memory_limit_gib=args.memory_limit_gib,
        )
    except (
        SourceFitError,
        model_input.ModelInputError,
        modeling.ModelingError,
        source_artifacts.SourceArtifactError,
        OSError,
        ValueError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
