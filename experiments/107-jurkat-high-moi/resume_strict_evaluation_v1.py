#!/usr/bin/env python3
"""One-shot recovery for the failed Experiment 107 strict evaluation.

The original runner, amendment, pretruth plan, outcome-access receipt, and
failed-run lock remain immutable. This entrypoint authenticates that complete
lineage, verifies the axis-only failure diagnosis, and permits exactly one
replay under a distinct exclusive lock.

Importing this module and running ``authenticate`` do not load any expression
or delta values. The real outcomes are loaded only by ``run`` after the
existing access receipt has authenticated.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import os
import stat
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np

import run_strict_evaluation as ORIGINAL
import strict_nadig_validation as STRICT


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXPERIMENT = "107-jurkat-high-moi"
RECOVERY_AMENDMENT_PATH = (
    HERE / "STRICT_EVALUATION_FAILURE_RECOVERY.v1.json"
)
EXPECTED_RECOVERY_AMENDMENT_SHA256 = (
    "73d0f30b5e9e9fa941624e99fa067d58c90f30777a14688abadd38fd06198c66"
)
RESUME_RUNNER_PATH = Path(__file__).resolve()
RECOVERY_LOCK_PATH = HERE / "work/.strict-evaluation-recovery-v1.lock"
RESULT_DIRECTORY = ORIGINAL.RESULT_DIRECTORY
EXPECTED_HELD_ONLY = ("TBCE-1", "HSPA14-1", "TMSB15B-1")


class RecoveryError(RuntimeError):
    """Raised when the recovery lineage or execution contract is violated."""


@dataclass(frozen=True)
class RecoveryLease:
    """Capability proving that this process owns the one-shot recovery lock."""

    path: Path
    descriptor: int
    device: int
    inode: int


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RecoveryError(message)


def _load_json(path: Path) -> dict[str, Any]:
    try:
        return ORIGINAL._load_json(path)  # noqa: SLF001
    except ORIGINAL.EvaluationDriverError as error:
        raise RecoveryError(str(error)) from error


def _verify_binding(
    root: Path,
    binding: Mapping[str, Any],
    label: str,
) -> Path:
    try:
        return ORIGINAL._verify_binding(root, binding, label)  # noqa: SLF001
    except ORIGINAL.EvaluationDriverError as error:
        raise RecoveryError(str(error)) from error


def _sha256_strings(values: Sequence[str]) -> str:
    return STRICT.sha256_strings(tuple(map(str, values)))


def _decode_strings(values: np.ndarray) -> tuple[str, ...]:
    return tuple(
        item.decode("utf-8") if isinstance(item, bytes) else str(item)
        for item in np.asarray(values).tolist()
    )


def _read_official_axis(path: Path) -> tuple[str, ...]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        _require(
            reader.fieldnames == ["gene_name"],
            "official gene-axis CSV schema changed",
        )
        values = tuple(str(row["gene_name"]) for row in reader)
    _require(values and len(values) == len(set(values)), "official axis is invalid")
    return values


def _read_h5ad_var_axis_metadata_only(path: Path) -> tuple[str, ...]:
    """Read only the H5AD var index; never access X, layers, or raw."""

    import h5py

    try:
        with h5py.File(path, "r") as handle:
            _require("var" in handle, f"{path.name} has no var metadata")
            group = handle["var"]
            key = group.attrs.get("_index", "_index")
            if isinstance(key, bytes):
                key = key.decode("utf-8")
            _require(isinstance(key, str) and key in group, "H5AD var index changed")
            node = group[key]
            if isinstance(node, h5py.Group):
                _require(
                    set(node.keys()) == {"mask", "values"},
                    "nullable H5AD var index schema changed",
                )
                mask = np.asarray(node["mask"][...], dtype=bool)
                _require(not np.any(mask), "H5AD var index contains missing symbols")
                raw = node["values"][...]
            else:
                raw = node[...]
    except (OSError, KeyError, TypeError, ValueError) as error:
        raise RecoveryError(f"cannot read H5AD var metadata from {path}: {error}") from error
    values = _decode_strings(raw)
    _require(values and len(values) == len(set(values)), "held axis is invalid")
    return values


def _read_npz_axes_only(path: Path, names: Sequence[str]) -> dict[str, tuple[str, ...]]:
    """Read only named axis members from an NPZ archive."""

    try:
        with np.load(path, allow_pickle=False) as archive:
            _require(
                all(name in archive.files for name in names),
                f"{path.name} lacks a required axis",
            )
            return {name: _decode_strings(archive[name]) for name in names}
    except (OSError, ValueError, KeyError) as error:
        raise RecoveryError(f"cannot read axes from {path}: {error}") from error


def build_axis_evidence(
    *,
    official_genes: Sequence[str],
    held_genes: Sequence[str],
    reference_genes: Sequence[str],
    integration_genes: Sequence[str],
    score_genes: Sequence[str],
) -> dict[str, Any]:
    """Build deterministic exact-symbol intersection evidence from axes only."""

    official = tuple(map(str, official_genes))
    held = tuple(map(str, held_genes))
    reference = tuple(map(str, reference_genes))
    integration = tuple(map(str, integration_genes))
    score = tuple(map(str, score_genes))
    for label, values in (
        ("official", official),
        ("held", held),
        ("reference", reference),
        ("integration", integration),
        ("score", score),
    ):
        _require(values and len(values) == len(set(values)), f"{label} axis is invalid")
    _require(reference == held, "reference-DE axis differs from held metadata")
    _require(integration == official, "integration axis differs from official axis")

    official_set = set(official)
    held_set = set(held)
    representable = tuple(gene for gene in official if gene in held_set)
    held_only = tuple(gene for gene in held if gene not in official_set)
    official_only = tuple(gene for gene in official if gene not in held_set)
    _require(
        all(gene in held_set and gene in official_set for gene in score),
        "PDS score axis is not fully representable",
    )
    _require("HSPA14" in representable, "exact HSPA14 symbol is not representable")

    return {
        "held_axis_sha256": _sha256_strings(held),
        "held_gene_count": len(held),
        "held_only_gene_count": len(held_only),
        "held_only_genes": list(held_only),
        "held_only_genes_sha256": _sha256_strings(held_only),
        "hspa14_exact_symbol": {
            "held_alias_index": (
                held.index("HSPA14-1") if "HSPA14-1" in held_set else None
            ),
            "held_index": held.index("HSPA14"),
            "official_index": official.index("HSPA14"),
            "representable_index": representable.index("HSPA14"),
        },
        "integration_axis_sha256": _sha256_strings(integration),
        "integration_gene_count": len(integration),
        "official_axis_sha256": _sha256_strings(official),
        "official_gene_count": len(official),
        "official_only_gene_count": len(official_only),
        "official_only_genes_sha256": _sha256_strings(official_only),
        "reference_axis_sha256": _sha256_strings(reference),
        "reference_gene_count": len(reference),
        "representable_axis_sha256": _sha256_strings(representable),
        "representable_gene_count": len(representable),
        "score_axis_all_representable": True,
        "score_axis_sha256": _sha256_strings(score),
        "score_gene_count": len(score),
    }


def _reconstruct_axis_evidence(
    amendment: Mapping[str, Any],
    *,
    root: Path,
) -> dict[str, dict[str, Any]]:
    sources = amendment["axis_sources"]
    official_path = _verify_binding(
        root,
        sources["official_gene_axis"],
        "recovery official gene axis",
    )
    official = _read_official_axis(official_path)
    observed: dict[str, dict[str, Any]] = {}
    for surface in STRICT.AUTHORIZING_HELD_SURFACES:
        bindings = sources[surface]
        _require(
            isinstance(bindings, Mapping)
            and set(bindings) == {"held_harness", "reference_de", "integration"},
            f"{surface} recovery axis-source schema changed",
        )
        held_path = _verify_binding(
            root, bindings["held_harness"], f"{surface} held harness"
        )
        reference_path = _verify_binding(
            root, bindings["reference_de"], f"{surface} reference DE"
        )
        integration_path = _verify_binding(
            root, bindings["integration"], f"{surface} integration"
        )
        held = _read_h5ad_var_axis_metadata_only(held_path)
        reference = _read_npz_axes_only(reference_path, ("genes",))["genes"]
        integration_axes = _read_npz_axes_only(
            integration_path, ("genes", "score_genes")
        )
        observed[surface] = build_axis_evidence(
            official_genes=official,
            held_genes=held,
            reference_genes=reference,
            integration_genes=integration_axes["genes"],
            score_genes=integration_axes["score_genes"],
        )
    return observed


def _authenticate_access_receipt(
    receipt: Mapping[str, Any],
    *,
    receipt_path: Path,
    plan: Mapping[str, Any],
    predecessors: Mapping[str, Any],
) -> None:
    _require(
        stat.S_IMODE(receipt_path.stat().st_mode) & 0o222 == 0,
        "outcome-access receipt is no longer immutable",
    )
    _require(
        receipt.get("schema_version") == 1
        and receipt.get("experiment") == EXPERIMENT
        and receipt.get("status") == "nadig_outcome_access_committed",
        "outcome-access receipt identity changed",
    )
    _require(
        receipt.get("benchmark_role")
        == "reused_development_benchmark_not_independent_confirmation",
        "outcome-access receipt benchmark role changed",
    )
    _require(
        receipt.get("artifact_bindings_sha256")
        == predecessors["artifact_bindings_sha256"]
        == plan["artifact_bindings"]["binding_sha256"],
        "outcome-access receipt artifact binding changed",
    )
    _require(
        receipt.get("pretruth_plan") == dict(predecessors["pretruth_plan"]),
        "outcome-access receipt pretruth-plan binding changed",
    )
    _require(
        receipt.get("validation_implementation")
        == dict(predecessors["strict_validator"]),
        "outcome-access receipt validator binding changed",
    )
    _require(
        receipt.get("nadig_artifacts") == plan["nadig_artifacts"],
        "outcome-access receipt Nadig binding changed",
    )
    _require(
        receipt.get("nadig_artifacts", {}).get("array_values_opened") is False,
        "outcome-access receipt no longer records pre-load commitment",
    )


def authenticate_recovery_amendment(
    amendment_path: Path = RECOVERY_AMENDMENT_PATH,
    *,
    root: Path = ROOT,
    expected_amendment_sha256: str = EXPECTED_RECOVERY_AMENDMENT_SHA256,
    require_replay_available: bool = True,
    recovery_lock_path: Path = RECOVERY_LOCK_PATH,
    result_directory: Path = RESULT_DIRECTORY,
) -> dict[str, Any]:
    """Authenticate the frozen failed-run lineage using metadata and axes only."""

    amendment_path = Path(amendment_path).resolve()
    _require(
        STRICT.sha256_file(amendment_path) == expected_amendment_sha256,
        "recovery amendment SHA-256 changed",
    )
    amendment = _load_json(amendment_path)
    _require(
        set(amendment)
        == {
            "schema_version",
            "experiment",
            "status",
            "failure",
            "predecessor_bindings",
            "axis_sources",
            "axis_evidence",
            "recovery_contract",
        },
        "recovery amendment top-level schema changed",
    )
    _require(
        amendment["schema_version"] == 1
        and amendment["experiment"] == EXPERIMENT
        and amendment["status"]
        == "frozen_recovery_candidate_for_independent_review",
        "recovery amendment identity changed",
    )
    failure = amendment["failure"]
    _require(
        failure
        == {
            "diagnosis_scope": "metadata_and_axes_only",
            "error": "h1 held genes are absent from official axis",
            "failed_run_lock_must_remain": (
                "experiments/107-jurkat-high-moi/work/"
                ".strict-evaluation-run.lock"
            ),
            "original_result_directory_absent_at_freeze": True,
        },
        "failure record changed",
    )
    predecessors = amendment["predecessor_bindings"]
    _require(
        isinstance(predecessors, Mapping)
        and set(predecessors)
        == {
            "artifact_bindings_sha256",
            "failed_run_lock",
            "original_execution_amendment",
            "original_runner",
            "outcome_access_receipt",
            "pretruth_plan",
            "strict_validator",
        },
        "predecessor binding schema changed",
    )
    paths = {
        name: _verify_binding(root, predecessors[name], name.replace("_", " "))
        for name in (
            "failed_run_lock",
            "original_execution_amendment",
            "original_runner",
            "outcome_access_receipt",
            "pretruth_plan",
            "strict_validator",
        )
    }
    _require(
        paths["original_runner"] == ORIGINAL.RUNNER_PATH,
        "recovery binds a different original runner",
    )
    _require(
        paths["strict_validator"] == Path(STRICT.__file__).resolve(),
        "recovery binds a different strict validator",
    )
    _require(
        paths["pretruth_plan"] == STRICT.PRETRUTH_PLAN_PATH.resolve(),
        "recovery binds a different pretruth plan",
    )
    _require(
        paths["outcome_access_receipt"]
        == STRICT.OUTCOME_ACCESS_RECEIPT_PATH.resolve(),
        "recovery binds a different outcome-access receipt",
    )
    _require(
        paths["failed_run_lock"] == ORIGINAL.LOCK_PATH.resolve(),
        "recovery binds a different failed-run lock",
    )

    try:
        original_auth = ORIGINAL.authenticate_execution_amendment(
            paths["original_execution_amendment"],
            root=root,
            runner_path=paths["original_runner"],
            plan_path=paths["pretruth_plan"],
            validator_path=paths["strict_validator"],
            expected_artifact_binding=str(
                predecessors["artifact_bindings_sha256"]
            ),
        )
    except ORIGINAL.EvaluationDriverError as error:
        raise RecoveryError(f"original execution lineage failed: {error}") from error

    receipt = _load_json(paths["outcome_access_receipt"])
    _authenticate_access_receipt(
        receipt,
        receipt_path=paths["outcome_access_receipt"],
        plan=original_auth["plan"],
        predecessors=predecessors,
    )
    observed_axes = _reconstruct_axis_evidence(amendment, root=root)
    _require(
        observed_axes == amendment["axis_evidence"],
        "axis-only failure evidence changed",
    )
    _require(
        tuple(observed_axes["h1"]["held_only_genes"]) == EXPECTED_HELD_ONLY,
        "H1 held-only mismatch is not the frozen three-label set",
    )
    _require(
        observed_axes["orion_hek293t"]["held_only_gene_count"] == 0
        and observed_axes["orion_hek293t"]["representable_gene_count"] == 18_533,
        "HEK293T is no longer an exact official-axis match",
    )

    contract = amendment["recovery_contract"]
    _require(
        contract
        == {
            "alias_merging": False,
            "directional_credit_axis": (
                "official_order_exact_symbol_intersection_with_held_axis"
            ),
            "existing_access_receipt_must_authenticate_before_outcome_reload": True,
            "frozen_components_unchanged": [
                "source_model",
                "integration_candidates",
                "fraction_grid",
                "nested_selection",
                "permutation_null",
                "cluster_bootstrap",
                "promotion_gates",
                "pds_score_axes",
            ],
            "new_exclusive_lock": (
                "experiments/107-jurkat-high-moi/work/"
                ".strict-evaluation-recovery-v1.lock"
            ),
            "original_result_directory": (
                "experiments/107-jurkat-high-moi/work/"
                "strict-evaluation-results"
            ),
            "pds_policy": (
                "byte_identical_frozen_score_axes_and_existing_pds_implementation"
            ),
            "replay_limit": 1,
            "result_directory_must_be_absent_before_replay": True,
        },
        "recovery execution contract changed",
    )
    expected_result = (Path(root) / contract["original_result_directory"]).resolve()
    expected_lock = (Path(root) / contract["new_exclusive_lock"]).resolve()
    _require(
        Path(result_directory).resolve() == expected_result,
        "recovery targets a different result directory",
    )
    _require(
        Path(recovery_lock_path).resolve() == expected_lock,
        "recovery uses a different exclusive lock",
    )
    _require(
        not expected_result.exists() and not expected_result.is_symlink(),
        "original result directory is no longer absent",
    )
    if require_replay_available:
        _require(
            not expected_lock.exists() and not expected_lock.is_symlink(),
            "the one allowed recovery replay has already been claimed",
        )
    else:
        _require(
            expected_lock.is_file()
            and not expected_lock.is_symlink()
            and expected_lock.stat().st_size == 0,
            "active recovery lock is missing or malformed",
        )
    return {
        "amendment": amendment,
        "amendment_path": amendment_path,
        "axis_evidence": observed_axes,
        "original_auth": original_auth,
        "paths": paths,
        "receipt": receipt,
        "result_directory": expected_result,
        "recovery_lock_path": expected_lock,
    }


def _reload_nadig_after_existing_receipt(
    authenticated: Mapping[str, Any],
    *,
    lease: RecoveryLease,
    outcome_loader: Callable[..., np.ndarray] = np.load,
) -> tuple[Any, dict[str, Any], np.ndarray]:
    """Reload Nadig only after the immutable access receipt has authenticated."""

    _require(
        lease.path.resolve() == authenticated["recovery_lock_path"],
        "outcome reload lacks the authorized recovery lock",
    )
    try:
        descriptor_stat = os.fstat(lease.descriptor)
        path_stat = lease.path.stat()
    except OSError as error:
        raise RecoveryError("recovery-lock lease is no longer valid") from error
    _require(
        descriptor_stat.st_dev == lease.device == path_stat.st_dev
        and descriptor_stat.st_ino == lease.inode == path_stat.st_ino,
        "recovery-lock lease does not match the one-shot lock",
    )
    source, plan = STRICT.authorize_outcome_phase(
        STRICT.SOURCE_DIRECTORY,
        STRICT.PRETRUTH_PLAN_PATH,
        STRICT.ARTIFACT_MANIFEST_PATH,
    )
    _require(
        plan == authenticated["original_auth"]["plan"],
        "pretruth plan changed before outcome reload",
    )
    receipt = authenticated["receipt"]
    _require(
        source.manifest_sha256 == receipt["source_manifest_sha256"]
        and source.model_sha256 == receipt["source_model_sha256"],
        "source model differs from the existing access receipt",
    )
    observed_nadig = STRICT.verify_nadig_artifacts(
        plan["nadig_artifacts"],
        STRICT.NADIG_DIRECTORY,
    )
    _require(
        observed_nadig == receipt["nadig_artifacts"],
        "Nadig artifacts differ from the existing access receipt",
    )
    outcome_path = STRICT.NADIG_DIRECTORY / "delta.npy"
    try:
        outcome = outcome_loader(outcome_path, allow_pickle=False)
    except (OSError, ValueError) as error:
        raise RecoveryError(f"cannot reload Nadig outcome: {error}") from error
    outcome = np.asarray(outcome)
    _require(
        outcome.shape
        == (int(plan["nadig_target_count"]), int(plan["nadig_gene_count"]))
        and np.all(np.isfinite(outcome)),
        "reloaded Nadig outcome matrix is invalid",
    )
    return source, plan, outcome


def authorize_then_reload_nadig(
    *,
    lease: RecoveryLease,
    outcome_loader: Callable[..., np.ndarray] = np.load,
    authentication_kwargs: Mapping[str, Any] | None = None,
) -> tuple[dict[str, Any], Any, dict[str, Any], np.ndarray]:
    """Authenticate recovery and its existing receipt before outcome reload."""

    kwargs = dict(authentication_kwargs or {})
    kwargs.setdefault("require_replay_available", False)
    kwargs.setdefault("recovery_lock_path", lease.path)
    authenticated = authenticate_recovery_amendment(
        **kwargs
    )
    source, plan, outcome = _reload_nadig_after_existing_receipt(
        authenticated,
        lease=lease,
        outcome_loader=outcome_loader,
    )
    return authenticated, source, plan, outcome


def exact_symbol_projection(
    official_genes: Sequence[str],
    held_genes: Sequence[str],
) -> tuple[np.ndarray, np.ndarray, tuple[str, ...], tuple[str, ...]]:
    """Return official/held indices for their ordered exact-symbol intersection."""

    official = tuple(map(str, official_genes))
    held = tuple(map(str, held_genes))
    _require(
        len(official) == len(set(official))
        and len(held) == len(set(held)),
        "projection axes contain duplicates",
    )
    held_lookup = {gene: index for index, gene in enumerate(held)}
    official_set = set(official)
    representable = tuple(gene for gene in official if gene in held_lookup)
    held_only = tuple(gene for gene in held if gene not in official_set)
    official_indices = np.asarray(
        [index for index, gene in enumerate(official) if gene in held_lookup],
        dtype=np.int64,
    )
    held_indices = np.asarray(
        [held_lookup[gene] for gene in representable], dtype=np.int64
    )
    _require(
        tuple(np.asarray(official)[official_indices]) == representable
        and tuple(np.asarray(held)[held_indices]) == representable,
        "exact-symbol projection is inconsistent",
    )
    return official_indices, held_indices, representable, held_only


def project_reference_de(
    reference: ORIGINAL.ReferenceDE,
    held_indices: np.ndarray,
) -> ORIGINAL.ReferenceDE:
    """Project reference DE onto selected held coordinates without alias merging."""

    selected = np.asarray(held_indices, dtype=np.int64)
    _require(
        selected.ndim == 1
        and np.all((0 <= selected) & (selected < len(reference.genes)))
        and len(selected) == len(set(map(int, selected))),
        "reference projection indices are invalid",
    )
    _require(
        np.array_equal(reference.gate_indices, np.flatnonzero(reference.gate)),
        "reference gate indices changed",
    )
    old_gate_lookup = {
        int(full_index): gate_index
        for gate_index, full_index in enumerate(reference.gate_indices)
    }
    old_gate_columns = np.asarray(
        [old_gate_lookup[int(index)] for index in selected if reference.gate[index]],
        dtype=np.int64,
    )
    remap = {
        int(old_gate): new_gate
        for new_gate, old_gate in enumerate(old_gate_columns)
    }
    significant = tuple(
        np.asarray(
            [remap[int(index)] for index in values if int(index) in remap],
            dtype=np.int64,
        )
        for values in reference.significant
    )
    gate = np.asarray(reference.gate[selected], dtype=bool)
    projected = ORIGINAL.ReferenceDE(
        genes=np.asarray(reference.genes[selected]).astype(str),
        gate=gate,
        gate_indices=np.flatnonzero(gate),
        control_mean_cpm=np.asarray(
            reference.control_mean_cpm[old_gate_columns], dtype=np.float64
        ),
        targets=np.asarray(reference.targets).astype(str),
        lfc=np.asarray(reference.lfc[:, old_gate_columns], dtype=np.float64),
        significant=significant,
    )
    _require(
        len(projected.gate_indices) == projected.lfc.shape[1]
        == len(projected.control_mean_cpm),
        "projected reference DE is malformed",
    )
    return projected


def pds_on_frozen_score_axis(
    *,
    incumbent_full: np.ndarray,
    candidate_full: np.ndarray,
    held_truth: np.ndarray,
    official_genes: Sequence[str],
    held_genes: Sequence[str],
    score_genes: Sequence[str],
    score_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply the original PDS implementation on the untouched frozen score axis."""

    official = np.asarray(official_genes).astype(str)
    held = np.asarray(held_genes).astype(str)
    score = np.asarray(score_genes).astype(str)
    indices = np.asarray(score_indices, dtype=np.int64)
    _require(
        score.tolist() == official[indices].tolist(),
        "integration PDS score symbols changed",
    )
    held_lookup = {gene: index for index, gene in enumerate(held)}
    _require(
        len(held_lookup) == len(held)
        and all(gene in held_lookup for gene in score),
        "held outcome lacks a frozen PDS score gene",
    )
    truth_indices = np.asarray(
        [held_lookup[gene] for gene in score], dtype=np.int64
    )
    truth = np.asarray(held_truth[:, truth_indices], dtype=np.float64)
    return (
        STRICT.pds_per_target(
            np.asarray(incumbent_full[:, indices], dtype=np.float64),
            truth,
        ),
        STRICT.pds_per_target(
            np.asarray(candidate_full[:, indices], dtype=np.float64),
            truth,
        ),
    )


def evaluate_held_surface_v1(
    surface: str,
    fraction: float,
    *,
    harness_path: Path,
    reference_de_path: Path,
) -> dict[str, Any]:
    """Evaluate one held surface with only the directional axis repaired."""

    arrays = ORIGINAL._load_integration(surface)  # noqa: SLF001
    ORIGINAL._validate_integration_candidates(arrays, surface=surface)  # noqa: SLF001
    targets = tuple(np.asarray(arrays["targets"]).astype(str))
    official_genes = np.asarray(arrays["genes"]).astype(str)
    score_genes = np.asarray(arrays["score_genes"]).astype(str)
    score_indices = np.asarray(arrays["score_indices"], dtype=np.int64)
    held_genes, control, held_truth = ORIGINAL._load_held_effects(  # noqa: SLF001
        harness_path, targets
    )
    reference = ORIGINAL._load_reference_de(reference_de_path)  # noqa: SLF001
    _require(
        np.array_equal(reference.genes, held_genes),
        f"{surface} reference-DE gene axis differs from held outcome",
    )
    suffix = ORIGINAL._fraction_suffix(fraction)  # noqa: SLF001
    incumbent_full = np.asarray(arrays["incumbent_effects"], dtype=np.float64)
    candidate_full = np.asarray(arrays[f"candidate_{suffix}"], dtype=np.float64)
    incumbent_pds, candidate_pds = pds_on_frozen_score_axis(
        incumbent_full=incumbent_full,
        candidate_full=candidate_full,
        held_truth=held_truth,
        official_genes=official_genes,
        held_genes=held_genes,
        score_genes=score_genes,
        score_indices=score_indices,
    )

    official_indices, held_indices, representable, held_only = (
        exact_symbol_projection(official_genes, held_genes)
    )
    projected_reference = project_reference_de(reference, held_indices)
    _require(
        tuple(projected_reference.genes.astype(str)) == representable,
        f"{surface} projected directional axis changed",
    )
    projected_control = np.asarray(control[held_indices], dtype=np.float64)
    incumbent_credit = ORIGINAL.directional_credit_per_target(
        incumbent_full[:, official_indices],
        projected_control,
        targets,
        projected_reference,
    )
    candidate_credit = ORIGINAL.directional_credit_per_target(
        candidate_full[:, official_indices],
        projected_control,
        targets,
        projected_reference,
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
        "directional_genes": np.asarray(representable, dtype="U"),
        "held_only_excluded_genes": np.asarray(held_only, dtype="U"),
    }


@contextmanager
def _one_shot_recovery_lock(
    path: Path = RECOVERY_LOCK_PATH,
) -> Iterable[RecoveryLease]:
    """Create one permanent attempt marker and hold its exclusive advisory lock."""

    lock = Path(path)
    lock.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(
            lock,
            os.O_CREAT | os.O_EXCL | os.O_WRONLY,
            0o444,
        )
    except FileExistsError as error:
        raise RecoveryError("the one allowed recovery replay was already claimed") from error
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        os.fsync(descriptor)
        parent_descriptor = os.open(lock.parent, os.O_RDONLY)
        try:
            os.fsync(parent_descriptor)
        finally:
            os.close(parent_descriptor)
        descriptor_stat = os.fstat(descriptor)
        yield RecoveryLease(
            path=lock.resolve(),
            descriptor=descriptor,
            device=descriptor_stat.st_dev,
            inode=descriptor_stat.st_ino,
        )
    finally:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            os.close(descriptor)


def _result_arrays(
    evaluation: Mapping[str, Any],
    held: Mapping[str, Mapping[str, Any]],
) -> dict[str, np.ndarray]:
    arrays = ORIGINAL._result_arrays(evaluation, held)  # noqa: SLF001
    for surface, values in held.items():
        arrays[f"{surface}_directional_genes"] = np.asarray(
            values["directional_genes"]
        ).astype("U")
        arrays[f"{surface}_held_only_excluded_genes"] = np.asarray(
            values["held_only_excluded_genes"]
        ).astype("U")
    return arrays


def run_recovery_evaluation() -> dict[str, Any]:
    """Perform the sole authorized replay and write the original result path."""

    with _one_shot_recovery_lock() as lease:
        authenticated = authenticate_recovery_amendment(
            require_replay_available=False
        )
        ORIGINAL.preflight_frozen_integrations()
        _source, plan, outcome = _reload_nadig_after_existing_receipt(
            authenticated,
            lease=lease,
        )
        galleries, full_integrations = ORIGINAL._load_nadig_galleries(  # noqa: SLF001
            plan, outcome
        )

        provisional_grid = STRICT.pds_fraction_grid(
            tuple(item.incumbent for item in galleries),
            tuple(item.source for item in galleries),
            tuple(item.truth for item in galleries),
            tuple(item.transferable for item in galleries),
            fractions=STRICT.FRACTIONS,
        )
        provisional_nested = STRICT.nested_blend_selection(provisional_grid)
        fallback = ORIGINAL._fallback_evidence(  # noqa: SLF001
            full_integrations,
            provisional_nested.fractions_by_held_gallery,
        )
        final = STRICT.final_fraction(provisional_grid)
        held = {
            surface: evaluate_held_surface_v1(
                surface,
                final,
                harness_path=authenticated["original_auth"]["held_paths"][surface][
                    "harness"
                ],
                reference_de_path=authenticated["original_auth"]["held_paths"][
                    surface
                ]["reference_de"],
            )
            for surface in STRICT.AUTHORIZING_HELD_SURFACES
        }
        held_metrics = {
            surface: values["metrics"] for surface, values in held.items()
        }
        evaluation = ORIGINAL.evaluate_nadig_statistics(
            galleries,
            fallback_evidence=fallback,
            held_context_metrics=held_metrics,
        )
        _require(
            evaluation["selected_fraction"] == final,
            "final fraction changed between frozen evaluation passes",
        )
        summary = {
            "schema_version": 1,
            "experiment": EXPERIMENT,
            "status": "strict_evaluation_complete_after_axis_recovery_v1",
            "created_at_utc": datetime.now(UTC).isoformat(),
            "benchmark_role": plan["benchmark_role"],
            "recovery_amendment": STRICT.file_binding(
                authenticated["amendment_path"], root=ROOT
            ),
            "resume_runner": STRICT.file_binding(RESUME_RUNNER_PATH, root=ROOT),
            "original_execution_amendment": STRICT.file_binding(
                authenticated["paths"]["original_execution_amendment"], root=ROOT
            ),
            "artifact_binding_sha256": authenticated["original_auth"][
                "artifact_binding_sha256"
            ],
            "axis_recovery": authenticated["axis_evidence"],
            "final_fraction": evaluation["selected_fraction"],
            "nested_selection": ORIGINAL._nested_payload(  # noqa: SLF001
                evaluation["nested"]
            ),
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
                authenticated["paths"]["outcome_access_receipt"], root=ROOT
            ),
            "failed_run_lock_preserved": STRICT.file_binding(
                authenticated["paths"]["failed_run_lock"], root=ROOT
            ),
            "submission_performed": False,
        }
        arrays = _result_arrays(evaluation, held)
        input_paths = [
            authenticated["amendment_path"],
            RESUME_RUNNER_PATH,
            authenticated["recovery_lock_path"],
            *authenticated["paths"].values(),
            STRICT.ARTIFACT_MANIFEST_PATH,
        ]
        axis_sources = authenticated["amendment"]["axis_sources"]
        input_paths.append(ROOT / axis_sources["official_gene_axis"]["path"])
        for surface in STRICT.AUTHORIZING_HELD_SURFACES:
            input_paths.extend(
                ROOT / axis_sources[surface][name]["path"]
                for name in ("held_harness", "reference_de", "integration")
            )
        for surface in STRICT.AUTHORIZING_HELD_SURFACES:
            input_paths.extend(
                authenticated["original_auth"]["held_paths"][surface][name]
                for name in ("harness", "reference_de")
            )
        unique_input_paths: list[Path] = []
        seen_inputs: set[Path] = set()
        for path in map(Path, input_paths):
            resolved = path.resolve()
            if resolved not in seen_inputs:
                seen_inputs.add(resolved)
                unique_input_paths.append(resolved)
        bound_inputs = [
            STRICT.file_binding(
                path,
                root=ROOT,
            )
            for path in unique_input_paths
        ]
        manifest = ORIGINAL.write_immutable_results(
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
        help="verify the failed-run lineage and axis-only recovery",
    )
    run = commands.add_parser("run", help="perform the sole recovery replay")
    run.add_argument(
        "--confirm-reload-held-outcomes",
        action="store_true",
        help="required acknowledgement that frozen outcomes will be reloaded",
    )
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.command == "authenticate":
        report = authenticate_recovery_amendment()
        print(
            json.dumps(
                {
                    "status": "recovery_amendment_authenticated",
                    "artifact_binding_sha256": report["original_auth"][
                        "artifact_binding_sha256"
                    ],
                    "h1_representable_gene_count": report["axis_evidence"]["h1"][
                        "representable_gene_count"
                    ],
                    "h1_held_only_excluded_genes": report["axis_evidence"]["h1"][
                        "held_only_genes"
                    ],
                    "result_directory_absent": True,
                    "recovery_replay_available": True,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    if not args.confirm_reload_held_outcomes:
        raise SystemExit("run requires --confirm-reload-held-outcomes")
    print(json.dumps(run_recovery_evaluation(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
