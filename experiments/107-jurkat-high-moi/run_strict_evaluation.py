#!/usr/bin/env python3
"""One-shot, authenticated execution driver for Experiment 107.

Importing this module is outcome-safe. The real Nadig and held-context outcomes
are touched only by ``run_real_evaluation`` after the execution amendment has
authenticated the frozen plan, validator, artifact binding, runner, and held
input bindings.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import math
import os
import shutil
import uuid
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np

import strict_nadig_validation as STRICT


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXPERIMENT = "107-jurkat-high-moi"
AMENDMENT_PATH = HERE / "STRICT_EVALUATION_EXECUTION_AMENDMENT.json"
RUNNER_PATH = Path(__file__).resolve()
RESULT_DIRECTORY = HERE / "work/strict-evaluation-results"
LOCK_PATH = HERE / "work/.strict-evaluation-run.lock"
EXPECTED_ARTIFACT_BINDING = (
    "b11deb9572f8a57a73c7518ed2e4ba26c078ac7cd0202fa077a4fa9987d3ecca"
)
RESULT_ARRAYS_NAME = "evaluation-arrays.npz"
RESULT_SUMMARY_NAME = "evaluation-summary.json"
RESULT_MANIFEST_NAME = "result-manifest.json"
TS = 50_000.0


class EvaluationDriverError(RuntimeError):
    """Raised when an execution or result invariant is violated."""


@dataclass(frozen=True)
class GalleryInput:
    incumbent: np.ndarray
    source: np.ndarray
    truth: np.ndarray
    transferable: np.ndarray
    cluster_basis: np.ndarray


@dataclass(frozen=True)
class ReferenceDE:
    genes: np.ndarray
    gate: np.ndarray
    gate_indices: np.ndarray
    control_mean_cpm: np.ndarray
    targets: np.ndarray
    lfc: np.ndarray
    significant: tuple[np.ndarray, ...]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise EvaluationDriverError(message)


def _load_json(path: Path) -> dict[str, Any]:
    path = Path(path)
    _require(path.is_file() and not path.is_symlink(), f"unsafe JSON file: {path}")
    try:
        with path.open(encoding="utf-8") as handle:
            value = json.load(
                handle,
                object_pairs_hook=STRICT._reject_duplicate_keys,  # noqa: SLF001
                parse_constant=STRICT._reject_constant,  # noqa: SLF001
            )
    except (OSError, ValueError, TypeError) as error:
        raise EvaluationDriverError(f"cannot load {path}: {error}") from error
    _require(isinstance(value, dict), f"{path} is not a JSON object")
    return value


def _confined_regular_file(root: Path, relative: str, label: str) -> Path:
    _require(
        isinstance(relative, str)
        and relative
        and relative == Path(relative).as_posix()
        and not Path(relative).is_absolute()
        and ".." not in Path(relative).parts,
        f"{label} path is not canonical repo-relative POSIX",
    )
    resolved_root = Path(root).resolve()
    candidate = resolved_root / relative
    _require(
        candidate.is_file() and not candidate.is_symlink(),
        f"{label} is missing or unsafe",
    )
    resolved = candidate.resolve()
    _require(
        resolved == resolved_root or resolved_root in resolved.parents,
        f"{label} escapes repository root",
    )
    return resolved


def _verify_binding(root: Path, binding: Mapping[str, Any], label: str) -> Path:
    _require(
        isinstance(binding, Mapping)
        and set(binding) == {"path", "bytes", "sha256"},
        f"{label} binding schema changed",
    )
    path = _confined_regular_file(root, str(binding["path"]), label)
    _require(
        path.stat().st_size == int(binding["bytes"]),
        f"{label} byte count changed",
    )
    _require(
        STRICT.sha256_file(path) == str(binding["sha256"]),
        f"{label} SHA-256 changed",
    )
    return path


def authenticate_execution_amendment(
    amendment_path: Path = AMENDMENT_PATH,
    *,
    root: Path = ROOT,
    runner_path: Path = RUNNER_PATH,
    plan_path: Path = STRICT.PRETRUTH_PLAN_PATH,
    validator_path: Path | None = None,
    expected_artifact_binding: str = EXPECTED_ARTIFACT_BINDING,
    verify_held_inputs: bool = True,
) -> dict[str, Any]:
    """Authenticate every execution binding without parsing outcome arrays."""

    validator = (
        Path(STRICT.__file__).resolve()
        if validator_path is None
        else Path(validator_path).resolve()
    )
    amendment = _load_json(amendment_path)
    _require(
        set(amendment)
        == {
            "schema_version",
            "experiment",
            "status",
            "bindings",
            "execution_contract",
        },
        "execution amendment top-level schema changed",
    )
    _require(amendment["schema_version"] == 1, "execution amendment version changed")
    _require(amendment["experiment"] == EXPERIMENT, "execution experiment changed")
    _require(
        amendment["status"] == "frozen_before_experiment_107_outcome_execution",
        "execution amendment is not frozen",
    )
    bindings = amendment["bindings"]
    _require(
        isinstance(bindings, Mapping)
        and set(bindings)
        == {
            "pretruth_plan",
            "strict_validator",
            "strict_artifact_binding_sha256",
            "runner",
            "final_provenance_receipt",
            "held_inputs",
        },
        "execution amendment binding schema changed",
    )
    bound_plan = _verify_binding(root, bindings["pretruth_plan"], "pretruth plan")
    bound_validator = _verify_binding(
        root, bindings["strict_validator"], "strict validator"
    )
    bound_runner = _verify_binding(root, bindings["runner"], "evaluation runner")
    _require(
        bound_plan == Path(plan_path).resolve(),
        "execution amendment binds a different pretruth plan",
    )
    _require(
        bound_validator == validator,
        "execution amendment binds a different strict validator",
    )
    _require(
        bound_runner == Path(runner_path).resolve(),
        "execution amendment binds a different runner",
    )
    artifact_binding = str(bindings["strict_artifact_binding_sha256"])
    _require(
        artifact_binding == expected_artifact_binding,
        "execution amendment artifact binding changed",
    )
    receipt_path = _verify_binding(
        root,
        bindings["final_provenance_receipt"],
        "final provenance receipt",
    )
    receipt = _load_json(receipt_path)
    _require(
        receipt.get("status")
        == "btbd11_validation_and_provenance_repair_complete",
        "final provenance receipt status changed",
    )
    _require(
        receipt.get("artifact_binding_sha256") == artifact_binding,
        "final provenance receipt artifact binding disagrees with amendment",
    )
    plan = _load_json(bound_plan)
    _require(plan.get("status") == "pretruth_complete", "pretruth plan is not final")
    _require(
        plan.get("nadig_outcome_opened") is False,
        "pretruth plan records prior Nadig outcome access",
    )
    _require(
        plan.get("artifact_bindings", {}).get("binding_sha256") == artifact_binding,
        "pretruth plan artifact binding disagrees with amendment",
    )
    _require(
        plan.get("validation_implementation") == dict(bindings["strict_validator"]),
        "pretruth plan validator binding disagrees with amendment",
    )
    held_inputs = bindings["held_inputs"]
    _require(
        isinstance(held_inputs, Mapping)
        and set(held_inputs) == set(STRICT.AUTHORIZING_HELD_SURFACES),
        "execution amendment held-input set changed",
    )
    held_paths: dict[str, dict[str, Path]] = {}
    for surface in STRICT.AUTHORIZING_HELD_SURFACES:
        surface_bindings = held_inputs[surface]
        _require(
            isinstance(surface_bindings, Mapping)
            and set(surface_bindings) == {"harness", "reference_de"},
            f"{surface} held-input schema changed",
        )
        if verify_held_inputs:
            held_paths[surface] = {
                name: _verify_binding(
                    root, surface_bindings[name], f"{surface} {name}"
                )
                for name in ("harness", "reference_de")
            }
    execution = amendment["execution_contract"]
    _require(
        execution
        == {
            "fractions": list(STRICT.FRACTIONS),
            "gallery_count": 6,
            "targets_per_gallery": 300,
            "permutation_replicates": STRICT.PERMUTATION_REPLICATES,
            "cluster_bootstrap_replicates": STRICT.BOOTSTRAP_REPLICATES,
            "held_surfaces": list(STRICT.AUTHORIZING_HELD_SURFACES),
            "outcome_loader_must_follow_authentication": True,
            "immutable_result_directory": (
                "experiments/107-jurkat-high-moi/work/"
                "strict-evaluation-results"
            ),
        },
        "execution contract changed",
    )
    return {
        "amendment": amendment,
        "amendment_path": Path(amendment_path).resolve(),
        "plan": plan,
        "plan_path": bound_plan,
        "validator_path": bound_validator,
        "runner_path": bound_runner,
        "artifact_binding_sha256": artifact_binding,
        "held_paths": held_paths,
    }


def authorize_then_load_nadig(
    amendment_path: Path = AMENDMENT_PATH,
    *,
    outcome_loader: Callable[..., tuple[Any, dict[str, Any], np.ndarray]]
    | None = None,
    authentication_kwargs: Mapping[str, Any] | None = None,
) -> tuple[dict[str, Any], Any, dict[str, Any], np.ndarray]:
    """Authenticate the execution amendment before invoking the outcome loader."""

    authenticated = authenticate_execution_amendment(
        amendment_path, **dict(authentication_kwargs or {})
    )
    loader = (
        STRICT.load_nadig_outcome_after_authorization
        if outcome_loader is None
        else outcome_loader
    )
    source, plan, outcome = loader(
        STRICT.SOURCE_DIRECTORY,
        STRICT.PRETRUTH_PLAN_PATH,
        STRICT.NADIG_DIRECTORY,
        artifact_manifest_path=STRICT.ARTIFACT_MANIFEST_PATH,
        receipt_path=STRICT.OUTCOME_ACCESS_RECEIPT_PATH,
    )
    return authenticated, source, plan, outcome


def _fraction_suffix(fraction: float) -> str:
    _require(float(fraction) in STRICT.FRACTIONS, "fraction is outside frozen grid")
    return f"f{int(round(float(fraction) * 1000)):03d}"


def _nested_payload(nested: Any) -> dict[str, Any]:
    return {
        "fractions_by_held_gallery": list(
            map(float, nested.fractions_by_held_gallery)
        ),
        "baseline_gallery_means": np.asarray(
            nested.baseline_gallery_means, dtype=np.float64
        ).tolist(),
        "selected_gallery_means": np.asarray(
            nested.selected_gallery_means, dtype=np.float64
        ).tolist(),
        "held_gallery_gains": np.asarray(
            nested.held_gallery_gains, dtype=np.float64
        ).tolist(),
        "macro_baseline": float(nested.macro_baseline),
        "macro_selected": float(nested.macro_selected),
        "macro_gain": float(nested.macro_gain),
    }


def cosine_per_target(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    """Own-target cosine, with zero-norm rows reported as zero."""

    predicted = np.asarray(prediction, dtype=np.float64)
    observed = np.asarray(truth, dtype=np.float64)
    _require(
        predicted.ndim == 2 and predicted.shape == observed.shape,
        "cosine arrays do not align",
    )
    numerator = np.einsum("ij,ij->i", predicted, observed)
    denominator = np.linalg.norm(predicted, axis=1) * np.linalg.norm(
        observed, axis=1
    )
    return np.divide(
        numerator,
        denominator,
        out=np.zeros(predicted.shape[0], dtype=np.float64),
        where=denominator > 0,
    )


def evaluate_nadig_statistics(
    galleries: Sequence[GalleryInput],
    *,
    fallback_evidence: Mapping[str, Any],
    held_context_metrics: Mapping[str, Any],
    strict_module: Any = STRICT,
) -> dict[str, Any]:
    """Run the complete frozen statistical protocol on six in-memory galleries."""

    _require(len(galleries) == 6, "evaluation requires exactly six galleries")
    incumbents = tuple(
        np.asarray(item.incumbent, dtype=np.float64) for item in galleries
    )
    sources = tuple(np.asarray(item.source, dtype=np.float64) for item in galleries)
    truths = tuple(np.asarray(item.truth, dtype=np.float64) for item in galleries)
    transferable = tuple(
        np.asarray(item.transferable, dtype=bool) for item in galleries
    )
    for index, (base, source, truth, support, item) in enumerate(
        zip(incumbents, sources, truths, transferable, galleries, strict=True),
        start=1,
    ):
        _require(
            base.shape == source.shape == truth.shape,
            f"gallery {index} score arrays do not align",
        )
        _require(
            base.ndim == 2
            and base.shape[0] == 300
            and support.shape == (300,),
            f"gallery {index} does not contain the frozen 300 targets",
        )
        _require(
            np.asarray(item.cluster_basis).shape == base.shape,
            f"gallery {index} cluster basis does not align",
        )
    pds_grid = strict_module.pds_fraction_grid(
        incumbents,
        sources,
        truths,
        transferable,
        fractions=strict_module.FRACTIONS,
    )
    nested = strict_module.nested_blend_selection(
        pds_grid, fractions=strict_module.FRACTIONS
    )
    selected_fraction = strict_module.final_fraction(
        pds_grid, fractions=strict_module.FRACTIONS
    )
    permutation_gains = strict_module.permutation_null_distribution(
        incumbents,
        sources,
        truths,
        transferable,
        repetitions=strict_module.PERMUTATION_REPLICATES,
        fractions=strict_module.FRACTIONS,
    )
    cluster_labels = tuple(
        strict_module.reciprocal_knn_clusters(
            np.asarray(item.cluster_basis, dtype=np.float64),
            gallery_one_based=index,
        )
        for index, item in enumerate(galleries, start=1)
    )
    bootstrap = strict_module.gallery_stratified_cluster_bootstrap(
        pds_grid,
        cluster_labels,
        repetitions=strict_module.BOOTSTRAP_REPLICATES,
        fractions=strict_module.FRACTIONS,
    )
    selected_deltas = tuple(
        np.asarray(pds_grid[fraction][index], dtype=np.float64)
        - np.asarray(pds_grid[0.0][index], dtype=np.float64)
        for index, fraction in enumerate(nested.fractions_by_held_gallery)
    )
    concentration = strict_module.gain_concentration(
        selected_deltas,
        cluster_labels,
        nested.macro_gain,
    )
    pure_source_per_target = tuple(
        strict_module.pds_per_target(source, truth, support)
        for source, truth, support in zip(
            sources, truths, transferable, strict=True
        )
    )
    pure_source_gallery_pds = np.asarray(
        [values.mean() for values in pure_source_per_target], dtype=np.float64
    )
    pure_source_per_target_cosine = tuple(
        cosine_per_target(source, truth)
        for source, truth in zip(sources, truths, strict=True)
    )
    pure_source_gallery_cosine = np.asarray(
        [values.mean() for values in pure_source_per_target_cosine],
        dtype=np.float64,
    )
    gate_report = strict_module.promotion_gate_report(
        pure_source_gallery_pds,
        nested,
        bootstrap,
        permutation_gains,
        concentration,
        held_context_metrics,
        fallback_evidence=fallback_evidence,
    )
    return {
        "pds_grid": pds_grid,
        "nested": nested,
        "selected_fraction": float(selected_fraction),
        "permutation_gains": np.asarray(permutation_gains, dtype=np.float64),
        "cluster_labels": cluster_labels,
        "bootstrap": bootstrap,
        "selected_deltas": selected_deltas,
        "concentration": concentration,
        "pure_source_per_target": pure_source_per_target,
        "pure_source_gallery_pds": pure_source_gallery_pds,
        "pure_source_per_target_cosine": pure_source_per_target_cosine,
        "pure_source_gallery_cosine": pure_source_gallery_cosine,
        "gate_report": gate_report,
    }


def _load_integration(surface: str) -> dict[str, np.ndarray]:
    path = HERE / "work/strict-validation" / surface / "integration.npz"
    _require(path.is_file() and not path.is_symlink(), f"missing integration: {surface}")
    try:
        with np.load(path, allow_pickle=False) as archive:
            return {name: np.asarray(archive[name]) for name in archive.files}
    except (OSError, ValueError, KeyError) as error:
        raise EvaluationDriverError(
            f"cannot load frozen integration {surface}: {error}"
        ) from error


def _validate_integration_candidates(
    arrays: Mapping[str, np.ndarray],
    *,
    surface: str,
) -> None:
    required = {
        "targets",
        "genes",
        "score_genes",
        "score_indices",
        "incumbent_effects",
        "source_effects",
        "transferable",
        "coordinate_supported",
        *(
            name
            for fraction in STRICT.FRACTIONS
            for name in (
                f"candidate_{_fraction_suffix(fraction)}",
                f"effective_{_fraction_suffix(fraction)}",
            )
        ),
    }
    _require(set(arrays) == required, f"{surface} integration schema changed")
    base = np.asarray(arrays["incumbent_effects"])
    source = np.asarray(arrays["source_effects"])
    transferable = np.asarray(arrays["transferable"], dtype=bool)
    coordinate_supported = np.asarray(arrays["coordinate_supported"], dtype=bool)
    for fraction in STRICT.FRACTIONS:
        suffix = _fraction_suffix(fraction)
        expected, effective = STRICT.blend_with_source(
            base,
            source,
            fraction,
            transferable,
            coordinate_supported,
        )
        _require(
            np.array_equal(expected, arrays[f"candidate_{suffix}"]),
            f"{surface} frozen candidate {fraction} disagrees with strict blend",
        )
        _require(
            np.array_equal(effective, arrays[f"effective_{suffix}"]),
            f"{surface} frozen effective mask {fraction} changed",
        )


def preflight_frozen_integrations() -> None:
    """Fully load and validate all eight source-only integrations."""

    surfaces = (
        *(f"nadig/gallery_{number}" for number in range(1, 7)),
        *STRICT.AUTHORIZING_HELD_SURFACES,
    )
    for surface in surfaces:
        arrays = _load_integration(surface)
        _validate_integration_candidates(arrays, surface=surface)
        del arrays


def _load_nadig_galleries(
    plan: Mapping[str, Any],
    outcome: np.ndarray,
) -> tuple[list[GalleryInput], list[dict[str, np.ndarray]]]:
    _require(len(plan.get("galleries", ())) == 6, "pretruth plan lost six galleries")
    galleries: list[GalleryInput] = []
    full: list[dict[str, np.ndarray]] = []
    for number, axes in enumerate(plan["galleries"], start=1):
        _require(int(axes["gallery"]) == number, "pretruth gallery order changed")
        arrays = _load_integration(f"nadig/gallery_{number}")
        _validate_integration_candidates(arrays, surface=f"nadig.gallery_{number}")
        targets = np.asarray(arrays["targets"]).astype(str)
        score_genes = np.asarray(arrays["score_genes"]).astype(str)
        score_indices = np.asarray(arrays["score_indices"], dtype=np.int64)
        _require(
            targets.tolist() == list(axes["targets"]),
            f"gallery {number} target axis differs from pretruth plan",
        )
        _require(
            score_genes.tolist() == list(axes["score_genes"])
            and score_indices.tolist() == list(axes["score_indices_in_official"]),
            f"gallery {number} score axis differs from pretruth plan",
        )
        truth = np.asarray(
            outcome[
                np.ix_(
                    np.asarray(axes["target_indices_in_nadig"], dtype=np.int64),
                    np.asarray(axes["score_indices_in_nadig"], dtype=np.int64),
                )
            ],
            dtype=np.float64,
        )
        incumbent = np.asarray(
            arrays["incumbent_effects"][:, score_indices], dtype=np.float64
        )
        source = np.asarray(
            arrays["source_effects"][:, score_indices], dtype=np.float64
        )
        transferable = np.asarray(arrays["transferable"], dtype=bool)
        _require(
            transferable.tolist() == list(axes["source_transferable"]),
            f"gallery {number} transferability differs from pretruth plan",
        )
        galleries.append(
            GalleryInput(
                incumbent=incumbent,
                source=source,
                truth=truth,
                transferable=transferable,
                cluster_basis=incumbent,
            )
        )
        full.append(arrays)
    return galleries, full


def _fallback_evidence(
    integrations: Sequence[Mapping[str, np.ndarray]],
    fractions_by_gallery: Sequence[float],
) -> dict[str, tuple[np.ndarray, ...]]:
    _require(
        len(integrations) == len(fractions_by_gallery) == 6,
        "fallback evidence requires six selected galleries",
    )
    candidates = []
    effective = []
    for arrays, fraction in zip(integrations, fractions_by_gallery, strict=True):
        suffix = _fraction_suffix(float(fraction))
        candidates.append(np.asarray(arrays[f"candidate_{suffix}"]))
        effective.append(np.asarray(arrays[f"effective_{suffix}"], dtype=bool))
    return {
        "incumbents": tuple(
            np.asarray(arrays["incumbent_effects"]) for arrays in integrations
        ),
        "candidates": tuple(candidates),
        "transferable": tuple(effective),
        "coordinate_supported": tuple(
            np.asarray(arrays["coordinate_supported"], dtype=bool)
            for arrays in integrations
        ),
    }


def _sum_rows(matrix: Any, rows: np.ndarray, genes: int, chunk_rows: int = 512) -> np.ndarray:
    output = np.zeros(genes, dtype=np.float64)
    ordered = np.sort(np.asarray(rows, dtype=np.int64))
    _require(ordered.size > 0, "held outcome group has no cells")
    for start in range(0, len(ordered), chunk_rows):
        block = matrix[ordered[start : start + chunk_rows]]
        output += np.asarray(block.sum(axis=0), dtype=np.float64).ravel()
    _require(np.all(np.isfinite(output)) and output.sum() > 0, "invalid held sum")
    return output


def _load_held_effects(
    harness_path: Path,
    expected_targets: Sequence[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read one held H5AD after authorization and return genes, control, effects."""

    import anndata as ad

    handle = ad.read_h5ad(harness_path, backed="r")
    try:
        genes = np.asarray(handle.var_names).astype(str)
        labels = np.asarray(handle.obs["target_gene"].astype(str))
        targets = tuple(map(str, expected_targets))
        observed = tuple(sorted(set(labels) - {"non-targeting"}))
        _require(observed == targets, f"{harness_path.name} target axis changed")
        control_sum = _sum_rows(
            handle.X, np.flatnonzero(labels == "non-targeting"), len(genes)
        )
        control = np.log1p(TS * control_sum / control_sum.sum())
        effects = np.empty((len(targets), len(genes)), dtype=np.float64)
        for row, target in enumerate(targets):
            total = _sum_rows(handle.X, np.flatnonzero(labels == target), len(genes))
            effects[row] = np.log1p(TS * total / total.sum()) - control
    finally:
        handle.file.close()
    _require(np.all(np.isfinite(effects)), "held effect matrix is non-finite")
    return genes, control, effects


def _load_reference_de(path: Path) -> ReferenceDE:
    try:
        with np.load(path, allow_pickle=False) as archive:
            genes = archive["genes"].astype(str)
            gate = archive["gate"].astype(bool)
            control = archive["ref_lin"].astype(np.float64)
            targets = archive["perts"].astype(str)
            lfc = archive["lfc"].astype(np.float64)
            pointers = archive["sig_ptr"].astype(np.int64)
            indices = archive["sig_idx"].astype(np.int64)
    except (OSError, ValueError, KeyError) as error:
        raise EvaluationDriverError(f"cannot load reference DE {path}: {error}") from error
    gate_indices = np.flatnonzero(gate)
    _require(len(genes) == len(gate), "reference-DE gene axis is malformed")
    _require(len(gate_indices) == len(control), "reference-DE control axis is malformed")
    _require(
        lfc.shape == (len(targets), len(gate_indices)),
        "reference-DE LFC matrix is malformed",
    )
    _require(
        pointers.shape == (len(targets) + 1,)
        and pointers[0] == 0
        and pointers[-1] == len(indices)
        and np.all(np.diff(pointers) >= 0),
        "reference-DE significance pointers are malformed",
    )
    significant = tuple(
        indices[pointers[row] : pointers[row + 1]] for row in range(len(targets))
    )
    return ReferenceDE(
        genes=genes,
        gate=gate,
        gate_indices=gate_indices,
        control_mean_cpm=control,
        targets=targets,
        lfc=lfc,
        significant=significant,
    )


def directional_credit_per_target(
    effect: np.ndarray,
    control_log_profile: np.ndarray,
    targets: Sequence[str],
    reference: ReferenceDE,
) -> np.ndarray:
    """Existing Experiment-106 directional-credit definition, kept exact."""

    values = np.asarray(effect, dtype=np.float64)
    control = np.asarray(control_log_profile, dtype=np.float64)
    _require(
        values.ndim == 2 and values.shape[1] == len(control),
        "directional effect and control axes do not align",
    )
    profile = np.expm1(control[None, :] + values) / TS
    np.clip(profile, 0.0, None, out=profile)
    predicted_lfc = np.log2(
        (1e6 * profile[:, reference.gate] + 1e-9)
        / (reference.control_mean_cpm[None, :] + 1e-9)
    )
    reference_rows = {
        str(target): index for index, target in enumerate(reference.targets)
    }
    result = np.full(len(targets), np.nan, dtype=np.float64)
    for row, target in enumerate(map(str, targets)):
        _require(target in reference_rows, f"reference DE is missing target {target}")
        reference_row = reference_rows[target]
        significant = np.asarray(
            reference.significant[reference_row], dtype=np.int64
        )
        full_indices = reference.gate_indices[significant]
        truth = reference.lfc[reference_row, significant]
        keep = (
            (reference.genes[full_indices] != target)
            & np.isfinite(truth)
            & (truth != 0)
        )
        significant = significant[keep]
        truth = truth[keep]
        _require(significant.size > 0, f"{target} has no directional-credit genes")
        result[row] = float(
            np.sign(truth) @ predicted_lfc[row, significant]
        ) / max(float(np.abs(truth).sum()), 1e-300)
    _require(np.all(np.isfinite(result)), "directional credit is non-finite")
    return result


def evaluate_held_surface(
    surface: str,
    fraction: float,
    *,
    harness_path: Path,
    reference_de_path: Path,
) -> dict[str, Any]:
    arrays = _load_integration(surface)
    _validate_integration_candidates(arrays, surface=surface)
    targets = tuple(np.asarray(arrays["targets"]).astype(str))
    official_genes = np.asarray(arrays["genes"]).astype(str)
    score_genes = np.asarray(arrays["score_genes"]).astype(str)
    score_indices = np.asarray(arrays["score_indices"], dtype=np.int64)
    held_genes, control, held_truth = _load_held_effects(harness_path, targets)
    reference = _load_reference_de(reference_de_path)
    _require(
        np.array_equal(reference.genes, held_genes),
        f"{surface} reference-DE gene axis differs from held outcome",
    )
    official_index = {gene: index for index, gene in enumerate(official_genes)}
    _require(
        len(official_index) == len(official_genes),
        f"{surface} official gene axis contains duplicates",
    )
    missing = [gene for gene in held_genes if gene not in official_index]
    _require(not missing, f"{surface} held genes are absent from official axis")
    held_from_official = np.asarray(
        [official_index[gene] for gene in held_genes], dtype=np.int64
    )
    held_index = {gene: index for index, gene in enumerate(held_genes)}
    _require(
        score_genes.tolist() == official_genes[score_indices].tolist(),
        f"{surface} integration score symbols changed",
    )
    truth_score_indices = np.asarray(
        [held_index[gene] for gene in score_genes], dtype=np.int64
    )
    suffix = _fraction_suffix(fraction)
    incumbent_full = np.asarray(arrays["incumbent_effects"], dtype=np.float64)
    candidate_full = np.asarray(arrays[f"candidate_{suffix}"], dtype=np.float64)
    incumbent_pds = STRICT.pds_per_target(
        incumbent_full[:, score_indices],
        held_truth[:, truth_score_indices],
    )
    candidate_pds = STRICT.pds_per_target(
        candidate_full[:, score_indices],
        held_truth[:, truth_score_indices],
    )
    incumbent_credit = directional_credit_per_target(
        incumbent_full[:, held_from_official],
        control,
        targets,
        reference,
    )
    candidate_credit = directional_credit_per_target(
        candidate_full[:, held_from_official],
        control,
        targets,
        reference,
    )
    return {
        "metrics": {
            "pds_delta": float(np.mean(candidate_pds - incumbent_pds)),
            "directional_credit_delta": float(
                np.mean(candidate_credit - incumbent_credit)
            ),
        },
        "incumbent_pds": incumbent_pds,
        "candidate_pds": candidate_pds,
        "incumbent_directional_credit": incumbent_credit,
        "candidate_directional_credit": candidate_credit,
        "targets": np.asarray(targets, dtype="U"),
        "score_genes": score_genes,
    }


def _json_safe(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float):
        _require(math.isfinite(value), "result JSON contains a non-finite float")
    return value


def _write_immutable_npz(path: Path, arrays: Mapping[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    _require(not path.exists() and not path.is_symlink(), f"replay refused: {path}")
    try:
        with path.open("xb") as handle:
            np.savez_compressed(
                handle,
                **{
                    name: np.ascontiguousarray(np.asarray(value))
                    for name, value in sorted(arrays.items())
                },
            )
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(path, 0o444)
    except BaseException:
        try:
            path.chmod(0o644)
        except OSError:
            pass
        path.unlink(missing_ok=True)
        raise


def _logical_binding(actual_path: Path, final_path: Path) -> dict[str, Any]:
    return {
        "path": final_path.relative_to(ROOT).as_posix(),
        "bytes": actual_path.stat().st_size,
        "sha256": STRICT.sha256_file(actual_path),
    }


def write_immutable_results(
    result_directory: Path,
    *,
    summary: Mapping[str, Any],
    arrays: Mapping[str, np.ndarray],
    input_bindings: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Atomically install one immutable result directory and refuse replay."""

    destination = Path(result_directory).resolve()
    _require(
        HERE.resolve() in destination.parents,
        "result directory must remain below Experiment 107",
    )
    _require(
        not destination.exists() and not destination.is_symlink(),
        f"immutable result directory already exists; replay refused: {destination}",
    )
    staging = destination.with_name(f".{destination.name}.{uuid.uuid4().hex}.tmp")
    _require(not staging.exists(), f"temporary result directory exists: {staging}")
    staging.mkdir(parents=True, mode=0o700)
    try:
        arrays_path = staging / RESULT_ARRAYS_NAME
        summary_path = staging / RESULT_SUMMARY_NAME
        manifest_path = staging / RESULT_MANIFEST_NAME
        _write_immutable_npz(arrays_path, arrays)
        STRICT.write_immutable_json(summary_path, _json_safe(summary))
        manifest = {
            "schema_version": 1,
            "experiment": EXPERIMENT,
            "status": "immutable_strict_evaluation_complete",
            "outputs": {
                "arrays": _logical_binding(
                    arrays_path, destination / RESULT_ARRAYS_NAME
                ),
                "summary": _logical_binding(
                    summary_path, destination / RESULT_SUMMARY_NAME
                ),
            },
            "inputs": list(input_bindings),
        }
        manifest["binding_sha256"] = STRICT.sha256_json_payload(manifest)
        STRICT.write_immutable_json(manifest_path, manifest)
        os.chmod(staging, 0o555)
        os.replace(staging, destination)
        parent_fd = os.open(destination.parent, os.O_RDONLY)
        try:
            os.fsync(parent_fd)
        finally:
            os.close(parent_fd)
        return manifest
    except BaseException:
        if staging.exists():
            os.chmod(staging, 0o700)
            for path in staging.iterdir():
                try:
                    path.chmod(0o600)
                except OSError:
                    pass
            shutil.rmtree(staging)
        raise


@contextmanager
def _execution_lock() -> Iterable[None]:
    LOCK_PATH.parent.mkdir(parents=True, exist_ok=True)
    with LOCK_PATH.open("a+") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise EvaluationDriverError(
                "another Experiment-107 evaluation process holds the lock"
            ) from error
        yield


def _result_arrays(
    evaluation: Mapping[str, Any],
    held: Mapping[str, Mapping[str, Any]],
) -> dict[str, np.ndarray]:
    fractions = tuple(map(float, STRICT.FRACTIONS))
    arrays: dict[str, np.ndarray] = {
        "fractions": np.asarray(fractions, dtype=np.float64),
        "pds_grid": np.stack(
            [
                np.stack(evaluation["pds_grid"][fraction])
                for fraction in fractions
            ]
        ),
        "nested_selected_fractions": np.asarray(
            evaluation["nested"].fractions_by_held_gallery, dtype=np.float64
        ),
        "nested_baseline_gallery_means": np.asarray(
            evaluation["nested"].baseline_gallery_means, dtype=np.float64
        ),
        "nested_selected_gallery_means": np.asarray(
            evaluation["nested"].selected_gallery_means, dtype=np.float64
        ),
        "nested_held_gallery_gains": np.asarray(
            evaluation["nested"].held_gallery_gains, dtype=np.float64
        ),
        "selected_per_target_deltas": np.stack(evaluation["selected_deltas"]),
        "pure_source_per_target_pds": np.stack(
            evaluation["pure_source_per_target"]
        ),
        "pure_source_gallery_pds": np.asarray(
            evaluation["pure_source_gallery_pds"], dtype=np.float64
        ),
        "pure_source_per_target_cosine": np.stack(
            evaluation["pure_source_per_target_cosine"]
        ),
        "pure_source_gallery_cosine": np.asarray(
            evaluation["pure_source_gallery_cosine"], dtype=np.float64
        ),
        "permutation_null_gains": np.asarray(
            evaluation["permutation_gains"], dtype=np.float64
        ),
        "cluster_bootstrap_draws": np.asarray(
            evaluation["bootstrap"]["draws"], dtype=np.float64
        ),
    }
    for index, labels in enumerate(evaluation["cluster_labels"], start=1):
        arrays[f"gallery_{index}_cluster_labels"] = np.asarray(labels).astype("U")
    for surface, values in held.items():
        arrays[f"{surface}_targets"] = np.asarray(values["targets"]).astype("U")
        arrays[f"{surface}_score_genes"] = np.asarray(
            values["score_genes"]
        ).astype("U")
        for name in (
            "incumbent_pds",
            "candidate_pds",
            "incumbent_directional_credit",
            "candidate_directional_credit",
        ):
            arrays[f"{surface}_{name}"] = np.asarray(values[name], dtype=np.float64)
    return arrays


def run_real_evaluation() -> dict[str, Any]:
    """Consume all frozen outcomes exactly once and write immutable results."""

    with _execution_lock():
        _require(
            not RESULT_DIRECTORY.exists() and not RESULT_DIRECTORY.is_symlink(),
            "strict evaluation result already exists; replay refused",
        )
        authenticated = authenticate_execution_amendment()
        preflight_frozen_integrations()
        _source, plan, outcome = STRICT.load_nadig_outcome_after_authorization(
            STRICT.SOURCE_DIRECTORY,
            STRICT.PRETRUTH_PLAN_PATH,
            STRICT.NADIG_DIRECTORY,
            artifact_manifest_path=STRICT.ARTIFACT_MANIFEST_PATH,
            receipt_path=STRICT.OUTCOME_ACCESS_RECEIPT_PATH,
        )
        _require(
            plan == authenticated["plan"],
            "authorized plan changed before Nadig outcome load",
        )
        galleries, full_integrations = _load_nadig_galleries(plan, outcome)

        provisional_grid = STRICT.pds_fraction_grid(
            tuple(item.incumbent for item in galleries),
            tuple(item.source for item in galleries),
            tuple(item.truth for item in galleries),
            tuple(item.transferable for item in galleries),
            fractions=STRICT.FRACTIONS,
        )
        provisional_nested = STRICT.nested_blend_selection(provisional_grid)
        fallback = _fallback_evidence(
            full_integrations,
            provisional_nested.fractions_by_held_gallery,
        )
        final = STRICT.final_fraction(provisional_grid)
        held = {
            surface: evaluate_held_surface(
                surface,
                final,
                harness_path=authenticated["held_paths"][surface]["harness"],
                reference_de_path=authenticated["held_paths"][surface][
                    "reference_de"
                ],
            )
            for surface in STRICT.AUTHORIZING_HELD_SURFACES
        }
        held_metrics = {
            surface: values["metrics"] for surface, values in held.items()
        }
        evaluation = evaluate_nadig_statistics(
            galleries,
            fallback_evidence=fallback,
            held_context_metrics=held_metrics,
        )
        _require(
            evaluation["selected_fraction"] == final,
            "final fraction changed between evaluation passes",
        )
        created_at = datetime.now(UTC).isoformat()
        summary = {
            "schema_version": 1,
            "experiment": EXPERIMENT,
            "status": "strict_evaluation_complete",
            "created_at_utc": created_at,
            "benchmark_role": plan["benchmark_role"],
            "execution_amendment": STRICT.file_binding(
                authenticated["amendment_path"], root=ROOT
            ),
            "artifact_binding_sha256": authenticated[
                "artifact_binding_sha256"
            ],
            "final_fraction": evaluation["selected_fraction"],
            "nested_selection": _nested_payload(evaluation["nested"]),
            "pds_fraction_gallery_means": {
                str(fraction): [
                    float(np.mean(values))
                    for values in evaluation["pds_grid"][fraction]
                ]
                for fraction in STRICT.FRACTIONS
            },
            "pure_source_gallery_pds": evaluation[
                "pure_source_gallery_pds"
            ].tolist(),
            "pure_source_gallery_cosine": evaluation[
                "pure_source_gallery_cosine"
            ].tolist(),
            "permutation_null": {
                "replicates": STRICT.PERMUTATION_REPLICATES,
                "maximum_gain": float(np.max(evaluation["permutation_gains"])),
            },
            "cluster_bootstrap": {
                key: value
                for key, value in evaluation["bootstrap"].items()
                if key != "draws"
            },
            "gain_concentration": evaluation["concentration"],
            "held_context_metrics": held_metrics,
            "promotion": evaluation["gate_report"],
            "outcome_access_receipt": STRICT.file_binding(
                STRICT.OUTCOME_ACCESS_RECEIPT_PATH, root=ROOT
            ),
            "submission_performed": False,
        }
        arrays = _result_arrays(evaluation, held)
        bound_inputs = [
            STRICT.file_binding(authenticated["amendment_path"], root=ROOT),
            STRICT.file_binding(authenticated["plan_path"], root=ROOT),
            STRICT.file_binding(authenticated["validator_path"], root=ROOT),
            STRICT.file_binding(authenticated["runner_path"], root=ROOT),
            STRICT.file_binding(STRICT.ARTIFACT_MANIFEST_PATH, root=ROOT),
            STRICT.file_binding(STRICT.OUTCOME_ACCESS_RECEIPT_PATH, root=ROOT),
        ]
        for surface in STRICT.AUTHORIZING_HELD_SURFACES:
            bound_inputs.extend(
                STRICT.file_binding(
                    authenticated["held_paths"][surface][name], root=ROOT
                )
                for name in ("harness", "reference_de")
            )
        manifest = write_immutable_results(
            RESULT_DIRECTORY,
            summary=summary,
            arrays=arrays,
            input_bindings=bound_inputs,
        )
        return {
            "status": summary["status"],
            "passes": bool(evaluation["gate_report"]["passes"]),
            "final_fraction": evaluation["selected_fraction"],
            "result_directory": RESULT_DIRECTORY.relative_to(ROOT).as_posix(),
            "result_binding_sha256": manifest["binding_sha256"],
        }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    commands.add_parser(
        "authenticate",
        help="verify execution bindings, including opaque held-file hashes",
    )
    run = commands.add_parser("run", help="perform the one-shot real evaluation")
    run.add_argument(
        "--confirm-open-held-outcomes",
        action="store_true",
        help="required acknowledgement that the one-shot outcome phase will begin",
    )
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.command == "authenticate":
        report = authenticate_execution_amendment()
        print(
            json.dumps(
                {
                    "status": "execution_amendment_authenticated",
                    "artifact_binding_sha256": report[
                        "artifact_binding_sha256"
                    ],
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    if not args.confirm_open_held_outcomes:
        raise SystemExit("run requires --confirm-open-held-outcomes")
    print(json.dumps(run_real_evaluation(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
