#!/usr/bin/env python3
"""Second one-shot recovery for Experiment 107 directional-credit empty sets.

The v1 recovery failed after outcome reload because PCCB had no eligible
reference directional-credit genes after the frozen filters. This entrypoint
preserves the complete v1 lineage and changes only the mathematically defined
value of an empty directional-credit sum: zero absolute credit for both
incumbent and candidate, while retaining the target and its equal weight.

Importing this module and running ``authenticate`` do not load any real
expression or delta arrays. Real outcomes are reloaded only by ``run`` after
the original access receipt and both failed-run lineages authenticate.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import os
import stat
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np

import resume_strict_evaluation_v1 as V1
import run_strict_evaluation as ORIGINAL
import strict_nadig_validation as STRICT


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXPERIMENT = "107-jurkat-high-moi"
V2_AMENDMENT_PATH = HERE / "STRICT_EVALUATION_FAILURE_RECOVERY.v2.json"
EXPECTED_V2_AMENDMENT_SHA256 = (
    "9ad8d1ed3fe3ce19cef7da6dd6b100d472a6b4311a73281aabfcbcc3824d6730"
)
V2_RUNNER_PATH = Path(__file__).resolve()
V2_LOCK_PATH = HERE / "work/.strict-evaluation-recovery-v2.lock"
V2_RESULT_DIRECTORY = HERE / "work/strict-evaluation-results-v2"
V1_RESULT_DIRECTORY = HERE / "work/strict-evaluation-results-v1"
EXPECTED_ERROR = "PCCB has no directional-credit genes"


class RecoveryV2Error(RuntimeError):
    """Raised when the v2 recovery lineage or execution contract is violated."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RecoveryV2Error(message)


def _load_json(path: Path) -> dict[str, Any]:
    try:
        return V1._load_json(path)  # noqa: SLF001
    except V1.RecoveryError as error:
        raise RecoveryV2Error(str(error)) from error


def _verify_binding(
    root: Path,
    binding: Mapping[str, Any],
    label: str,
) -> Path:
    try:
        return V1._verify_binding(root, binding, label)  # noqa: SLF001
    except V1.RecoveryError as error:
        raise RecoveryV2Error(str(error)) from error


def _expected_empty_set_policy() -> dict[str, Any]:
    return {
        "applies_uniformly_to_all_targets_and_authorizing_surfaces": True,
        "candidate_and_incumbent_absolute_credit": 0.0,
        "candidate_minus_incumbent_delta": 0.0,
        "denominator_definition": (
            "max(sum(abs(reference_lfc_over_eligible_set)),1e-300)"
        ),
        "eligible_set": (
            "reference_significant_genes_after_exact_axis_projection_"
            "excluding_the_perturbed_gene_nonfinite_lfc_and_zero_lfc"
        ),
        "empty_sum_definition": 0.0,
        "mathematical_extension": (
            "sum(empty)=0_and_0_divided_by_max(0,1e-300)=0"
        ),
        "nonempty_target_semantics": (
            "bit_identical_to_the_original_directional_credit_implementation"
        ),
        "pds_semantics": (
            "unchanged_v1_frozen_score_axes_and_original_pds_implementation"
        ),
        "target_handling": (
            "retain_every_target_in_original_order_and_equal_weighting_"
            "no_exclusion_no_reweighting_no_alias_merge"
        ),
    }


def _expected_recovery_contract() -> dict[str, Any]:
    return {
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
            "exact_symbol_axis_projection",
            "target_order",
            "target_weighting",
        ],
        "new_exclusive_lock": (
            "experiments/107-jurkat-high-moi/work/"
            ".strict-evaluation-recovery-v2.lock"
        ),
        "output_binding": {
            "arrays": ORIGINAL.RESULT_ARRAYS_NAME,
            "directory": (
                "experiments/107-jurkat-high-moi/work/"
                "strict-evaluation-results-v2"
            ),
            "manifest": ORIGINAL.RESULT_MANIFEST_NAME,
            "summary": ORIGINAL.RESULT_SUMMARY_NAME,
        },
        "predecessor_result_directories_must_be_absent": [
            "experiments/107-jurkat-high-moi/work/strict-evaluation-results",
            "experiments/107-jurkat-high-moi/work/strict-evaluation-results-v1",
        ],
        "replay_limit": 1,
        "v2_result_directory_must_be_absent_before_replay": True,
    }


def authenticate_v2_amendment(
    amendment_path: Path = V2_AMENDMENT_PATH,
    *,
    root: Path = ROOT,
    expected_amendment_sha256: str = EXPECTED_V2_AMENDMENT_SHA256,
    require_replay_available: bool = True,
    v2_lock_path: Path = V2_LOCK_PATH,
    v2_result_directory: Path = V2_RESULT_DIRECTORY,
) -> dict[str, Any]:
    """Authenticate v1's failed replay and the source-only v2 amendment."""

    amendment_path = Path(amendment_path).resolve()
    _require(
        STRICT.sha256_file(amendment_path) == expected_amendment_sha256,
        "v2 recovery amendment SHA-256 changed",
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
            "empty_set_policy",
            "recovery_contract",
        },
        "v2 recovery amendment top-level schema changed",
    )
    _require(
        amendment["schema_version"] == 2
        and amendment["experiment"] == EXPERIMENT
        and amendment["status"]
        == "frozen_v2_recovery_candidate_for_independent_review",
        "v2 recovery amendment identity changed",
    )
    _require(
        amendment["failure"]
        == {
            "exception_type": "EvaluationDriverError",
            "exact_error": EXPECTED_ERROR,
            "failing_callable": "ORIGINAL.directional_credit_per_target",
            "observation_source": "authorized_operator_report",
            "reported_stage": "after_outcome_reload_before_result_publication",
            "v1_result_published": False,
        },
        "v1 failure record changed",
    )
    _require(
        amendment["empty_set_policy"] == _expected_empty_set_policy(),
        "directional-credit empty-set policy changed",
    )
    contract = amendment["recovery_contract"]
    _require(
        contract == _expected_recovery_contract(),
        "v2 recovery execution contract changed",
    )

    predecessors = amendment["predecessor_bindings"]
    _require(
        isinstance(predecessors, Mapping)
        and set(predecessors)
        == {
            "artifact_bindings_sha256",
            "original_outcome_access_receipt",
            "strict_validator",
            "v1_failure_recovery_amendment",
            "v1_failed_recovery_lock",
            "v1_resume_runner",
        },
        "v2 predecessor binding schema changed",
    )
    paths = {
        name: _verify_binding(root, predecessors[name], name.replace("_", " "))
        for name in (
            "original_outcome_access_receipt",
            "strict_validator",
            "v1_failure_recovery_amendment",
            "v1_failed_recovery_lock",
            "v1_resume_runner",
        )
    }
    _require(
        paths["original_outcome_access_receipt"]
        == STRICT.OUTCOME_ACCESS_RECEIPT_PATH.resolve(),
        "v2 binds a different original outcome-access receipt",
    )
    _require(
        paths["strict_validator"] == Path(STRICT.__file__).resolve(),
        "v2 binds a different strict validator",
    )
    _require(
        paths["v1_failure_recovery_amendment"]
        == V1.RECOVERY_AMENDMENT_PATH.resolve(),
        "v2 binds a different v1 amendment",
    )
    _require(
        paths["v1_failed_recovery_lock"] == V1.RECOVERY_LOCK_PATH.resolve(),
        "v2 binds a different v1 failed-run lock",
    )
    _require(
        paths["v1_resume_runner"] == Path(V1.__file__).resolve(),
        "v2 binds a different v1 resume runner",
    )
    _require(
        stat.S_IMODE(paths["v1_failed_recovery_lock"].stat().st_mode) & 0o222 == 0,
        "v1 failed-run lock is no longer immutable",
    )
    _require(
        stat.S_IMODE(paths["original_outcome_access_receipt"].stat().st_mode)
        & 0o222
        == 0,
        "original outcome-access receipt is no longer immutable",
    )

    try:
        v1_auth = V1.authenticate_recovery_amendment(
            paths["v1_failure_recovery_amendment"],
            root=root,
            expected_amendment_sha256=str(
                predecessors["v1_failure_recovery_amendment"]["sha256"]
            ),
            require_replay_available=False,
            recovery_lock_path=paths["v1_failed_recovery_lock"],
            result_directory=ORIGINAL.RESULT_DIRECTORY,
        )
    except V1.RecoveryError as error:
        raise RecoveryV2Error(f"v1 recovery lineage failed: {error}") from error
    _require(
        v1_auth["original_auth"]["artifact_binding_sha256"]
        == predecessors["artifact_bindings_sha256"],
        "v1 artifact binding differs from v2 amendment",
    )

    expected_predecessor_results = tuple(
        (Path(root) / relative).resolve()
        for relative in contract["predecessor_result_directories_must_be_absent"]
    )
    _require(
        expected_predecessor_results
        == (ORIGINAL.RESULT_DIRECTORY.resolve(), V1_RESULT_DIRECTORY.resolve()),
        "v2 predecessor result paths changed",
    )
    for path in expected_predecessor_results:
        _require(
            not path.exists() and not path.is_symlink(),
            f"predecessor result directory is no longer absent: {path}",
        )
    expected_result = (
        Path(root) / contract["output_binding"]["directory"]
    ).resolve()
    expected_lock = (Path(root) / contract["new_exclusive_lock"]).resolve()
    _require(
        Path(v2_result_directory).resolve() == expected_result,
        "v2 targets a different result directory",
    )
    _require(
        Path(v2_lock_path).resolve() == expected_lock,
        "v2 uses a different exclusive lock",
    )
    _require(
        not expected_result.exists() and not expected_result.is_symlink(),
        "v2 result directory is no longer absent",
    )
    if require_replay_available:
        _require(
            not expected_lock.exists() and not expected_lock.is_symlink(),
            "the one allowed v2 replay has already been claimed",
        )
    else:
        _require(
            expected_lock.is_file()
            and not expected_lock.is_symlink()
            and expected_lock.stat().st_size == 0,
            "active v2 recovery lock is missing or malformed",
        )
    return {
        "amendment": amendment,
        "amendment_path": amendment_path,
        "contract": contract,
        "paths": paths,
        "v1_auth": v1_auth,
        "v2_lock_path": expected_lock,
        "v2_result_directory": expected_result,
    }


def _verify_v2_lease(
    authenticated: Mapping[str, Any],
    lease: V1.RecoveryLease,
) -> None:
    _require(
        lease.path.resolve() == authenticated["v2_lock_path"],
        "outcome reload lacks the authorized v2 recovery lock",
    )
    try:
        descriptor_stat = os.fstat(lease.descriptor)
        path_stat = lease.path.stat()
    except OSError as error:
        raise RecoveryV2Error("v2 recovery-lock lease is no longer valid") from error
    _require(
        descriptor_stat.st_dev == lease.device == path_stat.st_dev
        and descriptor_stat.st_ino == lease.inode == path_stat.st_ino,
        "v2 recovery-lock lease does not match the one-shot lock",
    )


def _reload_nadig_after_v2_authentication(
    authenticated: Mapping[str, Any],
    *,
    lease: V1.RecoveryLease,
    outcome_loader: Callable[..., np.ndarray] = np.load,
) -> tuple[Any, dict[str, Any], np.ndarray]:
    """Reload Nadig only after v2 and the original receipt authenticate."""

    _verify_v2_lease(authenticated, lease)
    v1_auth = authenticated["v1_auth"]
    source, plan = STRICT.authorize_outcome_phase(
        STRICT.SOURCE_DIRECTORY,
        STRICT.PRETRUTH_PLAN_PATH,
        STRICT.ARTIFACT_MANIFEST_PATH,
    )
    _require(
        plan == v1_auth["original_auth"]["plan"],
        "pretruth plan changed before v2 outcome reload",
    )
    receipt = v1_auth["receipt"]
    _require(
        source.manifest_sha256 == receipt["source_manifest_sha256"]
        and source.model_sha256 == receipt["source_model_sha256"],
        "source model differs from the original access receipt",
    )
    observed_nadig = STRICT.verify_nadig_artifacts(
        plan["nadig_artifacts"],
        STRICT.NADIG_DIRECTORY,
    )
    _require(
        observed_nadig == receipt["nadig_artifacts"],
        "Nadig artifacts differ from the original access receipt",
    )
    outcome_path = STRICT.NADIG_DIRECTORY / "delta.npy"
    try:
        outcome = outcome_loader(outcome_path, allow_pickle=False)
    except (OSError, ValueError) as error:
        raise RecoveryV2Error(f"cannot reload Nadig outcome: {error}") from error
    outcome = np.asarray(outcome)
    _require(
        outcome.shape
        == (int(plan["nadig_target_count"]), int(plan["nadig_gene_count"]))
        and np.all(np.isfinite(outcome)),
        "reloaded Nadig outcome matrix is invalid",
    )
    return source, plan, outcome


def authorize_then_reload_nadig_v2(
    *,
    lease: V1.RecoveryLease,
    outcome_loader: Callable[..., np.ndarray] = np.load,
    authentication_kwargs: Mapping[str, Any] | None = None,
) -> tuple[dict[str, Any], Any, dict[str, Any], np.ndarray]:
    """Authenticate the complete v2 lineage before invoking the outcome loader."""

    kwargs = dict(authentication_kwargs or {})
    kwargs.setdefault("require_replay_available", False)
    kwargs.setdefault("v2_lock_path", lease.path)
    authenticated = authenticate_v2_amendment(**kwargs)
    source, plan, outcome = _reload_nadig_after_v2_authentication(
        authenticated,
        lease=lease,
        outcome_loader=outcome_loader,
    )
    return authenticated, source, plan, outcome


def _directional_credit_with_empty_mask(
    effect: np.ndarray,
    control_log_profile: np.ndarray,
    targets: Sequence[str],
    reference: ORIGINAL.ReferenceDE,
) -> tuple[np.ndarray, np.ndarray]:
    """Original directional credit with the empty-sum value defined as zero."""

    values = np.asarray(effect, dtype=np.float64)
    control = np.asarray(control_log_profile, dtype=np.float64)
    _require(
        values.ndim == 2 and values.shape[1] == len(control),
        "directional effect and control axes do not align",
    )
    profile = np.expm1(control[None, :] + values) / ORIGINAL.TS
    np.clip(profile, 0.0, None, out=profile)
    predicted_lfc = np.log2(
        (1e6 * profile[:, reference.gate] + 1e-9)
        / (reference.control_mean_cpm[None, :] + 1e-9)
    )
    reference_rows = {
        str(target): index for index, target in enumerate(reference.targets)
    }
    result = np.zeros(len(targets), dtype=np.float64)
    empty = np.zeros(len(targets), dtype=bool)
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
        if significant.size == 0:
            empty[row] = True
            continue
        result[row] = float(
            np.sign(truth) @ predicted_lfc[row, significant]
        ) / max(float(np.abs(truth).sum()), 1e-300)
    _require(np.all(np.isfinite(result)), "directional credit is non-finite")
    return result, empty


def directional_credit_per_target_v2(
    effect: np.ndarray,
    control_log_profile: np.ndarray,
    targets: Sequence[str],
    reference: ORIGINAL.ReferenceDE,
) -> np.ndarray:
    """Return directional credit, assigning the empty eligible set credit zero."""

    values, _empty = _directional_credit_with_empty_mask(
        effect,
        control_log_profile,
        targets,
        reference,
    )
    return values


def evaluate_held_surface_v2(
    surface: str,
    fraction: float,
    *,
    harness_path: Path,
    reference_de_path: Path,
) -> dict[str, Any]:
    """Evaluate one held surface with v1 axes and v2 empty-set semantics."""

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
    incumbent_pds, candidate_pds = V1.pds_on_frozen_score_axis(
        incumbent_full=incumbent_full,
        candidate_full=candidate_full,
        held_truth=held_truth,
        official_genes=official_genes,
        held_genes=held_genes,
        score_genes=score_genes,
        score_indices=score_indices,
    )

    official_indices, held_indices, representable, held_only = (
        V1.exact_symbol_projection(official_genes, held_genes)
    )
    projected_reference = V1.project_reference_de(reference, held_indices)
    _require(
        tuple(projected_reference.genes.astype(str)) == representable,
        f"{surface} projected directional axis changed",
    )
    projected_control = np.asarray(control[held_indices], dtype=np.float64)
    incumbent_credit, incumbent_empty = _directional_credit_with_empty_mask(
        incumbent_full[:, official_indices],
        projected_control,
        targets,
        projected_reference,
    )
    candidate_credit, candidate_empty = _directional_credit_with_empty_mask(
        candidate_full[:, official_indices],
        projected_control,
        targets,
        projected_reference,
    )
    _require(
        np.array_equal(incumbent_empty, candidate_empty),
        f"{surface} empty directional set depends on prediction values",
    )
    empty_targets = np.asarray(targets, dtype="U")[incumbent_empty]
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
        "empty_directional_credit_targets": empty_targets,
    }


@contextmanager
def _one_shot_v2_lock(path: Path = V2_LOCK_PATH) -> Iterable[V1.RecoveryLease]:
    """Create the permanent v2 attempt marker and hold its exclusive lock."""

    lock = Path(path)
    lock.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(
            lock,
            os.O_CREAT | os.O_EXCL | os.O_WRONLY,
            0o444,
        )
    except FileExistsError as error:
        raise RecoveryV2Error("the one allowed v2 replay was already claimed") from error
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        os.fsync(descriptor)
        parent_descriptor = os.open(lock.parent, os.O_RDONLY)
        try:
            os.fsync(parent_descriptor)
        finally:
            os.close(parent_descriptor)
        descriptor_stat = os.fstat(descriptor)
        yield V1.RecoveryLease(
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


def _result_arrays_v2(
    evaluation: Mapping[str, Any],
    held: Mapping[str, Mapping[str, Any]],
) -> dict[str, np.ndarray]:
    arrays = V1._result_arrays(evaluation, held)  # noqa: SLF001
    for surface, values in held.items():
        arrays[f"{surface}_empty_directional_credit_targets"] = np.asarray(
            values["empty_directional_credit_targets"]
        ).astype("U")
    return arrays


def _unique_file_bindings(paths: Sequence[Path]) -> list[dict[str, Any]]:
    unique: list[Path] = []
    seen: set[Path] = set()
    for path in map(Path, paths):
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique.append(resolved)
    return [STRICT.file_binding(path, root=ROOT) for path in unique]


def run_v2_recovery_evaluation() -> dict[str, Any]:
    """Perform the sole v2 replay and publish only to the v2 result directory."""

    with _one_shot_v2_lock() as lease:
        authenticated = authenticate_v2_amendment(require_replay_available=False)
        ORIGINAL.preflight_frozen_integrations()
        _source, plan, outcome = _reload_nadig_after_v2_authentication(
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
        v1_auth = authenticated["v1_auth"]
        held = {
            surface: evaluate_held_surface_v2(
                surface,
                final,
                harness_path=v1_auth["original_auth"]["held_paths"][surface][
                    "harness"
                ],
                reference_de_path=v1_auth["original_auth"]["held_paths"][surface][
                    "reference_de"
                ],
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
            "schema_version": 2,
            "experiment": EXPERIMENT,
            "status": "strict_evaluation_complete_after_empty_set_recovery_v2",
            "created_at_utc": datetime.now(UTC).isoformat(),
            "benchmark_role": plan["benchmark_role"],
            "v2_recovery_amendment": STRICT.file_binding(
                authenticated["amendment_path"], root=ROOT
            ),
            "v2_resume_runner": STRICT.file_binding(V2_RUNNER_PATH, root=ROOT),
            "v1_failure": authenticated["amendment"]["failure"],
            "empty_set_policy": authenticated["amendment"]["empty_set_policy"],
            "artifact_binding_sha256": v1_auth["original_auth"][
                "artifact_binding_sha256"
            ],
            "axis_recovery": v1_auth["axis_evidence"],
            "empty_directional_credit_targets": {
                surface: values["empty_directional_credit_targets"].tolist()
                for surface, values in held.items()
            },
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
            "original_outcome_access_receipt": STRICT.file_binding(
                authenticated["paths"]["original_outcome_access_receipt"],
                root=ROOT,
            ),
            "v1_failed_recovery_lock": STRICT.file_binding(
                authenticated["paths"]["v1_failed_recovery_lock"],
                root=ROOT,
            ),
            "v2_recovery_lock": STRICT.file_binding(
                authenticated["v2_lock_path"],
                root=ROOT,
            ),
            "submission_performed": False,
        }

        input_paths = [
            authenticated["amendment_path"],
            V2_RUNNER_PATH,
            authenticated["v2_lock_path"],
            *authenticated["paths"].values(),
            STRICT.ARTIFACT_MANIFEST_PATH,
            *v1_auth["paths"].values(),
        ]
        axis_sources = v1_auth["amendment"]["axis_sources"]
        input_paths.append(ROOT / axis_sources["official_gene_axis"]["path"])
        for surface in STRICT.AUTHORIZING_HELD_SURFACES:
            input_paths.extend(
                ROOT / axis_sources[surface][name]["path"]
                for name in ("held_harness", "reference_de", "integration")
            )
            input_paths.extend(
                v1_auth["original_auth"]["held_paths"][surface][name]
                for name in ("harness", "reference_de")
            )
        manifest = ORIGINAL.write_immutable_results(
            V2_RESULT_DIRECTORY,
            summary=summary,
            arrays=_result_arrays_v2(evaluation, held),
            input_bindings=_unique_file_bindings(input_paths),
        )
        return {
            "status": summary["status"],
            "passes": bool(evaluation["gate_report"]["passes"]),
            "final_fraction": evaluation["selected_fraction"],
            "result_directory": V2_RESULT_DIRECTORY.relative_to(ROOT).as_posix(),
            "result_binding_sha256": manifest["binding_sha256"],
        }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    commands.add_parser(
        "authenticate",
        help="verify the v1 failure lineage and source-only v2 recovery",
    )
    run = commands.add_parser("run", help="perform the sole v2 recovery replay")
    run.add_argument(
        "--confirm-reload-held-outcomes",
        action="store_true",
        help="required acknowledgement that frozen outcomes will be reloaded",
    )
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.command == "authenticate":
        report = authenticate_v2_amendment()
        print(
            json.dumps(
                {
                    "status": "v2_recovery_amendment_authenticated",
                    "artifact_binding_sha256": report["v1_auth"]["original_auth"][
                        "artifact_binding_sha256"
                    ],
                    "exact_v1_error": report["amendment"]["failure"]["exact_error"],
                    "empty_set_absolute_credit": report["amendment"][
                        "empty_set_policy"
                    ]["candidate_and_incumbent_absolute_credit"],
                    "v2_result_directory_absent": True,
                    "v2_replay_available": True,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    if not args.confirm_reload_held_outcomes:
        raise SystemExit("run requires --confirm-reload-held-outcomes")
    print(json.dumps(run_v2_recovery_evaluation(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
