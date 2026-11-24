#!/usr/bin/env python3
"""Audit or execute Experiment 114's frozen HEK293T PDS-only diagnostic.

The audit command validates the complete source-only freeze and metadata-fixed
row/coordinate mappings without hashing or opening either held NPY array.  The
evaluate command requires explicit confirmation, creates and fsyncs an
exclusive access receipt, and only then authenticates and opens held truth.
No caller-supplied input or output paths are accepted.
"""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import importlib.util
import json
import os
import sys
from pathlib import Path
from types import ModuleType
from typing import Any, Mapping, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUN_PATH = HERE / "run.py"
FREEZER_PATH = HERE / "freeze_source_candidates_v1.py"
APPROVAL_PATH = HERE / "HELD_SURFACE_EVALUATION_APPROVAL.v1.json"
EXPECTED_APPROVAL_SHA256 = (
    "9fd85f9df8cce2d917dcd5844eacccf9ba40fff4d70046a34d68bd018547e72c"
)
EXPECTED_AUTHORIZATION = {
    "metadata_only_audit": True,
    "held_truth_access": True,
    "evaluation": True,
    "source_candidate_freeze": False,
    "model_fitting": False,
    "count_generation": False,
    "challenge_build": False,
    "submission": False,
    "commit": False,
    "push": False,
}
EXPECTED_SURFACE_ID = "orion_hek293t_full_outer_holdout"
EXPECTED_ROLE = "orion_hek293t"
EXPECTED_TARGET_COUNT = 35
EXPECTED_TARGET_ORDERED_SHA256 = (
    "e81c113e76a627b4d8d44a6998d33f3fc533ee11e539b413016552a0c35b3587"
)
EXPECTED_OUTPUT_DIRECTORY = (
    "experiments/114-depmap-coessentiality/work/"
    "held-surface-evaluation-v1/orion_hek293t_full"
)
EXPECTED_ACCESS_RECEIPT = f"{EXPECTED_OUTPUT_DIRECTORY}/OUTCOME_ACCESS_RECEIPT.v1.json"
EXPECTED_RESULT_PATH = f"{EXPECTED_OUTPUT_DIRECTORY}/result.json"


def _load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


RUN = _load_module("experiment_114_held_eval_v1_run", RUN_PATH)
FREEZER = _load_module("experiment_114_held_eval_v1_freezer", FREEZER_PATH)
BUILDER = FREEZER.BUILDER
ContractError = RUN.ContractError


def load_approval() -> dict[str, Any]:
    if RUN.sha256_file(APPROVAL_PATH) != EXPECTED_APPROVAL_SHA256:
        raise ContractError("held-evaluation approval amendment changed")
    approval = RUN.load_json(APPROVAL_PATH)
    validate_approval(approval)
    return approval


def validate_approval(approval: Mapping[str, Any]) -> None:
    expected_keys = {
        "schema_version",
        "experiment",
        "amendment_id",
        "date_approved",
        "stage",
        "status",
        "entrypoint_path",
        "predecessors",
        "frozen_bundles",
        "authorization",
        "surface",
        "invariants",
    }
    if set(approval) != expected_keys:
        raise ContractError("held-evaluation approval schema drifted")
    expected_identity = {
        "schema_version": 1,
        "experiment": "114-depmap-coessentiality",
        "amendment_id": "experiment-114-held-surface-evaluation-approval-v1",
        "date_approved": "2026-08-30",
        "stage": "held_surface_evaluation_approval",
        "status": (
            "approved_for_full_orion_hek293t_outer_holdout_pds_only_diagnostic_"
            "pending_independent_review_and_execution"
        ),
        "entrypoint_path": (
            "experiments/114-depmap-coessentiality/"
            "evaluate_frozen_candidates_v1.py"
        ),
    }
    for key, expected in expected_identity.items():
        if approval.get(key) != expected:
            raise ContractError(f"held-evaluation approval identity drifted: {key}")
    if approval.get("authorization") != EXPECTED_AUTHORIZATION:
        raise ContractError("held-evaluation authorization drifted")
    if tuple(approval.get("frozen_bundles", {})) != ("h1", "orion_hek293t"):
        raise ContractError("held-evaluation frozen bundle set drifted")
    expected_invariants = {
        "all_frozen_bundles_must_validate_before_truth_path_acceptance": True,
        "caller_supplied_paths_forbidden": True,
        "evaluation_is_diagnostic_only": True,
        "exclusive_fsynced_access_receipt_before_truth_hash_or_open": True,
        "frozen_predecessors_must_not_be_mutated": True,
        "minimum_changed_target_rows_for_pds": 2,
        "no_overwrite": True,
        "null_mappings_must_regenerate_before_truth_path_acceptance": True,
        "only_predeclared_target_and_gene_mappings_allowed": True,
        "pds_only": True,
        "promotion_authorized": False,
    }
    if approval.get("invariants") != expected_invariants:
        raise ContractError("held-evaluation invariants drifted")
    surface = approval.get("surface")
    if not isinstance(surface, dict):
        raise ContractError("held-evaluation surface is missing")
    expected_surface = {
        "id": EXPECTED_SURFACE_ID,
        "frozen_role": EXPECTED_ROLE,
        "diagnostic_only": True,
        "direction_valid": False,
    }
    for key, expected in expected_surface.items():
        if surface.get(key) != expected:
            raise ContractError(f"held-evaluation surface drifted: {key}")
    metric = surface.get("metric")
    if not isinstance(metric, dict) or metric != {
        "name": "raw_uncentered_pds",
        "centering": "none",
        "gallery": "exact_metadata_fixed_35_target_overlap",
        "ties": "midrank",
        "response_coordinates": (
            "truth_representable_official_genes_excluding_all_300_vcc_panel_targets"
        ),
        "promotion_authorized": False,
    }:
        raise ContractError("held-evaluation metric contract drifted")
    output = surface.get("output")
    if not isinstance(output, dict) or output != {
        "directory": EXPECTED_OUTPUT_DIRECTORY,
        "access_receipt": EXPECTED_ACCESS_RECEIPT,
        "result": EXPECTED_RESULT_PATH,
    }:
        raise ContractError("held-evaluation output paths drifted")
    truth = surface.get("truth")
    if not isinstance(truth, dict) or tuple(truth) != (
        "b_ctrl",
        "delta",
        "genes",
        "meta",
        "provenance_manifest",
        "targets",
    ):
        raise ContractError("held-evaluation truth bundle schema drifted")
    expected_truth_paths = {
        "b_ctrl": "data/donors/orion_hek293t/b_ctrl.npy",
        "delta": "data/donors/orion_hek293t/delta.npy",
        "genes": "data/donors/orion_hek293t/genes.txt",
        "meta": "data/donors/orion_hek293t/meta.json",
        "provenance_manifest": (
            "experiments/107-jurkat-high-moi/work/strict-validation/"
            "orion_hek293t/donor.json"
        ),
        "targets": "data/donors/orion_hek293t/perts.txt",
    }
    if {
        name: str(record.get("path", "")) for name, record in truth.items()
    } != expected_truth_paths:
        raise ContractError("held-evaluation truth paths drifted")
    target_mapping = surface.get("target_mapping")
    if not isinstance(target_mapping, dict):
        raise ContractError("held-evaluation target mapping is missing")
    targets = target_mapping.get("targets")
    missing = target_mapping.get("eligible_missing_targets")
    if (
        not isinstance(targets, list)
        or len(targets) != EXPECTED_TARGET_COUNT
        or len(set(targets)) != len(targets)
        or RUN.sha256_lines(targets) != EXPECTED_TARGET_ORDERED_SHA256
        or target_mapping.get("changed_target_count") != EXPECTED_TARGET_COUNT
        or target_mapping.get("eligible_overlap_count") != EXPECTED_TARGET_COUNT
        or target_mapping.get("eligible_target_count") != 41
        or target_mapping.get("ordered_sha256") != EXPECTED_TARGET_ORDERED_SHA256
    ):
        raise ContractError("held-evaluation target axis drifted")
    if (
        missing != ["C5orf22", "CDK8", "GLYR1", "IL10RB", "RFK", "UNK"]
        or RUN.sha256_lines(missing)
        != target_mapping.get("eligible_missing_ordered_sha256")
    ):
        raise ContractError("held-evaluation missing-target axis drifted")
    if EXPECTED_TARGET_COUNT < int(
        approval["invariants"]["minimum_changed_target_rows_for_pds"]
    ):
        raise ContractError("held-evaluation gallery is too small for PDS")


def _verify_file_record(
    record: Mapping[str, Any],
    label: str,
    *,
    hash_contents: bool,
) -> Path:
    required = {"path", "sha256"}
    if not required.issubset(record):
        raise ContractError(f"{label} lacks path or sha256")
    path = RUN.resolve_repository_path(str(record["path"]), label)
    RUN.reject_symlink_components(path, label)
    if not path.is_file():
        raise ContractError(f"{label} is missing: {path}")
    if "bytes" in record and path.stat().st_size != int(record["bytes"]):
        raise ContractError(f"{label} byte count drifted")
    if hash_contents and RUN.sha256_file(path) != record["sha256"]:
        raise ContractError(f"{label} hash drifted")
    return path


def validate_frozen_output_without_transitive_truth(
    directory: Path,
    *,
    role: str,
    bundle: Mapping[str, Any],
    protocol: Mapping[str, Any],
) -> None:
    """Validate one frozen payload without replaying role-source input hashes."""

    manifest = RUN.load_json(directory / "manifest.json")
    expected_identity = {
        "schema_version": 1,
        "experiment": "114-depmap-coessentiality",
        "stage": "source_only_freeze",
        "contains_held_truth": False,
        "sealed_before_held_truth": True,
        "role": role,
    }
    for key, expected in expected_identity.items():
        if manifest.get(key) != expected:
            raise ContractError(f"{role} frozen output identity drifted: {key}")
    if manifest.get("input_manifest_sha256") != bundle["role_manifest"]["sha256"]:
        raise ContractError(f"{role} frozen role-manifest binding drifted")
    if manifest.get("outputs") != {
        "arms.npz": bundle["arms"]["sha256"],
        "permutation_nulls.npz": bundle["permutation_nulls"]["sha256"],
    }:
        raise ContractError(f"{role} frozen output bindings drifted")

    with np.load(directory / "arms.npz", allow_pickle=False) as archive:
        required = {
            "targets",
            "genes",
            *(f"effect__{name}" for name in RUN.ARMS),
            *(f"confidence__{name}" for name in RUN.ARMS),
        }
        if set(archive.files) != required:
            raise ContractError(f"{role} frozen arm keys drifted")
        arrays = {name: archive[name] for name in required}
    targets = arrays["targets"].astype(str)
    genes = arrays["genes"].astype(str)
    canonical_targets = np.asarray(
        RUN.read_csv_column(RUN.PANEL_PATH, "target_gene"), dtype=str
    )
    canonical_genes = np.asarray(
        RUN.read_csv_column(RUN.GENE_NAMES_PATH, "gene_name"), dtype=str
    )
    if not np.array_equal(targets, canonical_targets):
        raise ContractError(f"{role} frozen target axis drifted")
    if not np.array_equal(genes, canonical_genes):
        raise ContractError(f"{role} frozen gene axis drifted")
    if (
        RUN.sha256_lines(targets) != manifest["panel_targets_ordered_sha256"]
        or RUN.sha256_lines(genes) != manifest["gene_axis_ordered_sha256"]
    ):
        raise ContractError(f"{role} frozen axis hash drifted")
    shape = (len(targets), len(genes))
    for name in RUN.ARMS:
        effect = arrays[f"effect__{name}"]
        confidence = arrays[f"confidence__{name}"]
        if effect.shape != shape or confidence.shape != shape:
            raise ContractError(f"{role} frozen arm shape drifted: {name}")
        if not np.isfinite(effect).all() or not np.isfinite(confidence).all():
            raise ContractError(f"{role} frozen arm is nonfinite: {name}")
        if RUN.sha256_array(effect) != manifest["arms"][name]["effect_sha256"]:
            raise ContractError(f"{role} frozen effect hash drifted: {name}")
        if (
            RUN.sha256_array(confidence)
            != manifest["arms"][name]["confidence_sha256"]
        ):
            raise ContractError(f"{role} frozen confidence hash drifted: {name}")
        if not np.array_equal(
            confidence, np.abs(np.asarray(effect, dtype=np.float64))
        ):
            raise ContractError(f"{role} confidence is not abs(effect): {name}")
    if not np.array_equal(
        arrays["effect__identity"], arrays["effect__current_fallback"]
    ):
        raise ContractError(f"{role} identity/current_fallback differ")

    with np.load(directory / "permutation_nulls.npz", allow_pickle=False) as archive:
        if set(archive.files) != {
            "eligible_rows",
            "effects",
            "confidence",
            "query_feature_targets",
        }:
            raise ContractError(f"{role} frozen null keys drifted")
        eligible_rows = np.asarray(archive["eligible_rows"], dtype=np.int64)
        null_effects = archive["effects"]
        null_confidence = archive["confidence"]
        query_targets = archive["query_feature_targets"].astype(str)
    expected_eligible = np.flatnonzero(
        np.isin(targets, np.asarray(protocol["eligibility"]["eligible_targets"]))
    )
    if not np.array_equal(eligible_rows, expected_eligible):
        raise ContractError(f"{role} frozen null eligible rows drifted")
    expected_shape = (RUN.NULL_DRAWS, len(eligible_rows), len(genes))
    if (
        null_effects.shape != expected_shape
        or null_confidence.shape != expected_shape
        or query_targets.shape != (RUN.NULL_DRAWS, len(eligible_rows))
    ):
        raise ContractError(f"{role} frozen null shapes drifted")
    if not np.isfinite(null_effects).all() or not np.isfinite(null_confidence).all():
        raise ContractError(f"{role} frozen null arrays are nonfinite")
    if (
        RUN.sha256_array(null_effects) != manifest["permutation_effects_sha256"]
        or RUN.sha256_array(null_confidence)
        != manifest["permutation_confidence_sha256"]
        or RUN.sha256_array(query_targets)
        != manifest["permutation_query_feature_targets_sha256"]
    ):
        raise ContractError(f"{role} frozen null logical hash drifted")
    if not np.array_equal(
        null_confidence, np.abs(np.asarray(null_effects, dtype=np.float64))
    ):
        raise ContractError(f"{role} frozen null confidence is not abs(effect)")
    regenerated = [
        RUN.permuted_query_labels(
            protocol["eligibility"]["incumbent_uncovered_targets"],
            protocol["eligibility"]["low_confidence_targets"],
            draw,
        )
        for draw in range(RUN.NULL_DRAWS)
    ]
    expected_queries = np.asarray(
        [
            [mapping[target] for target in targets[eligible_rows]]
            for mapping in regenerated
        ]
    )
    if not np.array_equal(query_targets, expected_queries):
        raise ContractError(f"{role} frozen null labels did not regenerate")


def verify_frozen_chain(
    approval: Mapping[str, Any],
) -> tuple[dict[str, tuple[dict[str, Any], dict[str, np.ndarray]]], dict[str, Any]]:
    for name, record in approval["predecessors"].items():
        RUN.verify_fixed_record(record, f"held-evaluation predecessor {name}")

    freeze_approval = FREEZER.load_approval()
    FREEZER.verify_predecessors(freeze_approval)
    effective_protocol, effective_boundary = FREEZER.effective_authorization(
        freeze_approval
    )
    loaded: dict[str, tuple[dict[str, Any], dict[str, np.ndarray]]] = {}
    original_loader = RUN.load_freeze_boundary
    try:
        RUN.load_freeze_boundary = lambda: copy.deepcopy(effective_boundary)
        for role in ("h1", "orion_hek293t"):
            bundle = approval["frozen_bundles"][role]
            expected_directory = RUN.resolve_repository_path(
                str(bundle["directory"]), f"{role} frozen directory"
            )
            expected_role_directory = RUN.resolve_repository_path(
                freeze_approval["approved_roles"][role]["role_directory"],
                f"{role} role directory",
            )
            for name in ("manifest", "arms", "permutation_nulls"):
                path = RUN.verify_fixed_record(
                    bundle[name],
                    f"{role} frozen {name}",
                    required_under=expected_directory,
                )
                if path.parent != expected_directory:
                    raise ContractError(f"{role} frozen {name} path drifted")
            for name, record in bundle["source_role_artifacts"].items():
                RUN.verify_fixed_record(
                    record,
                    f"{role} source role artifact {name}",
                    required_under=expected_role_directory,
                )
            if (
                RUN.verify_fixed_record(
                    bundle["role_manifest"],
                    f"{role} role manifest",
                    required_under=expected_role_directory,
                )
                != expected_role_directory / "manifest.json"
            ):
                raise ContractError(f"{role} role manifest path drifted")
            if role == "h1":
                validate_frozen_output_without_transitive_truth(
                    expected_directory,
                    role=role,
                    bundle=bundle,
                    protocol=effective_protocol,
                )
            else:
                loaded[role] = RUN.load_frozen_directory(expected_directory)
    finally:
        RUN.load_freeze_boundary = original_loader
    return loaded, effective_protocol


def _read_axis_metadata(
    record: Mapping[str, Any],
    label: str,
    *,
    unique: bool,
) -> tuple[Path, np.ndarray]:
    path = _verify_file_record(record, label, hash_contents=True)
    values = np.asarray(
        BUILDER.read_axis(path, unique=unique, label=label),
        dtype=str,
    )
    if len(values) != int(record["count"]):
        raise ContractError(f"{label} count drifted")
    if RUN.sha256_lines(values) != record["ordered_sha256"]:
        raise ContractError(f"{label} order drifted")
    return path, values


def _int64(values: Sequence[int]) -> np.ndarray:
    return np.asarray(values, dtype="<i8")


def validate_metadata_mapping(
    approval: Mapping[str, Any],
    frozen: Mapping[str, Any],
    arrays: Mapping[str, np.ndarray],
    protocol: Mapping[str, Any],
) -> dict[str, Any]:
    surface = approval["surface"]
    target_record = surface["truth"]["targets"]
    gene_record = surface["truth"]["genes"]
    _, truth_targets = _read_axis_metadata(
        target_record, "HEK293T truth target metadata", unique=True
    )
    _, truth_genes = _read_axis_metadata(
        gene_record, "HEK293T truth gene metadata", unique=False
    )
    for name in ("meta", "provenance_manifest"):
        _verify_file_record(
            surface["truth"][name],
            f"HEK293T truth metadata {name}",
            hash_contents=True,
        )

    frozen_targets = arrays["targets"].astype(str)
    frozen_genes = arrays["genes"].astype(str)
    if frozen.get("role") != EXPECTED_ROLE:
        raise ContractError("evaluation loaded the wrong frozen role")
    primary = np.asarray(arrays[f"effect__{RUN.PRIMARY_ARM}"])
    fallback = np.asarray(arrays["effect__current_fallback"])
    changed = np.any(primary != fallback, axis=1)
    truth_target_index = {target: row for row, target in enumerate(truth_targets)}
    observed_targets = np.asarray(
        [
            target
            for row, target in enumerate(frozen_targets)
            if changed[row] and target in truth_target_index
        ],
        dtype=str,
    )
    mapping = surface["target_mapping"]
    expected_targets = np.asarray(mapping["targets"], dtype=str)
    if not np.array_equal(observed_targets, expected_targets):
        raise ContractError("metadata-fixed changed-target overlap drifted")
    eligible_targets = tuple(protocol["eligibility"]["eligible_targets"])
    observed_missing = [
        target for target in eligible_targets if target not in truth_target_index
    ]
    if observed_missing != mapping["eligible_missing_targets"]:
        raise ContractError("metadata-fixed missing-target set drifted")
    if len(observed_targets) < int(
        approval["invariants"]["minimum_changed_target_rows_for_pds"]
    ):
        raise ContractError("fewer than two changed targets are available for PDS")
    frozen_rows = np.asarray(
        [int(np.flatnonzero(frozen_targets == target)[0]) for target in observed_targets],
        dtype="<i8",
    )
    truth_rows = np.asarray(
        [truth_target_index[target] for target in observed_targets], dtype="<i8"
    )
    eligible_rows = np.asarray(arrays["null_eligible_rows"], dtype=np.int64)
    eligible_slot = {int(row): slot for slot, row in enumerate(eligible_rows)}
    if any(int(row) not in eligible_slot for row in frozen_rows):
        raise ContractError("changed target is absent from the frozen null axis")
    null_rows = np.asarray(
        [eligible_slot[int(row)] for row in frozen_rows], dtype="<i8"
    )
    for name, observed in (
        ("frozen_row_indices", frozen_rows),
        ("truth_row_indices", truth_rows),
        ("null_row_indices", null_rows),
    ):
        expected = _int64(mapping[name])
        if not np.array_equal(observed, expected):
            raise ContractError(f"{name} drifted")
        if (
            RUN.sha256_array(observed)
            != mapping[f"{name}_int64_sha256"]
        ):
            raise ContractError(f"{name} logical hash drifted")

    canonical_columns, truth_columns = BUILDER.last_occurrence_map(
        truth_genes, frozen_genes
    )
    coordinate = surface["response_coordinate_mapping"]
    if len(canonical_columns) != int(coordinate["mapped_count"]):
        raise ContractError("truth/canonical mapped gene count drifted")
    if (
        RUN.sha256_array(np.asarray(canonical_columns, dtype="<i8"))
        != coordinate["canonical_indices_int64_sha256"]
        or RUN.sha256_array(np.asarray(truth_columns, dtype="<i8"))
        != coordinate["truth_gene_indices_int64_sha256"]
    ):
        raise ContractError("truth/canonical gene mapping drifted")
    mapped = np.zeros(len(frozen_genes), dtype=bool)
    mapped[canonical_columns] = True
    panel = np.asarray(
        RUN.read_csv_column(RUN.PANEL_PATH, "target_gene"), dtype=str
    )
    if RUN.sha256_lines(panel) != frozen["panel_targets_ordered_sha256"]:
        raise ContractError("frozen panel axis drifted")
    panel_mask = np.isin(frozen_genes, panel)
    if int(np.sum(mapped & panel_mask)) != len(panel):
        raise ContractError("truth gene axis no longer represents every panel target")
    allowed = mapped & ~panel_mask
    if (
        int(mapped.sum()) != coordinate["mapped_count"]
        or int((~mapped).sum()) != coordinate["unmapped_canonical_gene_count"]
        or RUN.sha256_array(mapped.astype(np.uint8))
        != coordinate["mapped_mask_uint8_sha256"]
        or int(allowed.sum()) != coordinate["allowed_count"]
        or RUN.sha256_array(allowed.astype(np.uint8))
        != coordinate["allowed_mask_uint8_sha256"]
        or RUN.sha256_array(np.flatnonzero(allowed).astype("<i8"))
        != coordinate["allowed_indices_int64_sha256"]
    ):
        raise ContractError("response-coordinate mapping drifted")
    for name in RUN.ARMS:
        effect = np.asarray(arrays[f"effect__{name}"])
        confidence = np.asarray(arrays[f"confidence__{name}"])
        if not np.isfinite(effect).all() or not np.isfinite(confidence).all():
            raise ContractError(f"evaluation frozen arm is nonfinite: {name}")
        if np.any(effect[:, panel_mask] != 0) or np.any(
            confidence[:, panel_mask] != 0
        ):
            raise ContractError(f"evaluation frozen arm retains panel coordinates: {name}")
    if not np.isfinite(arrays["null_effects"]).all() or not np.isfinite(
        arrays["null_confidence"]
    ).all():
        raise ContractError("evaluation frozen null arrays are nonfinite")
    if np.any(arrays["null_effects"][:, :, panel_mask] != 0) or np.any(
        arrays["null_confidence"][:, :, panel_mask] != 0
    ):
        raise ContractError("evaluation frozen nulls retain panel coordinates")

    regenerated = [
        RUN.permuted_query_labels(
            protocol["eligibility"]["incumbent_uncovered_targets"],
            protocol["eligibility"]["low_confidence_targets"],
            draw,
        )
        for draw in range(RUN.NULL_DRAWS)
    ]
    expected_null_labels = np.asarray(
        [
            [draw[target] for target in frozen_targets[eligible_rows]]
            for draw in regenerated
        ]
    )
    if not np.array_equal(
        arrays["null_query_feature_targets"].astype(str), expected_null_labels
    ):
        raise ContractError("permutation null labels did not regenerate")

    return {
        "targets": observed_targets,
        "frozen_rows": frozen_rows,
        "truth_rows": truth_rows,
        "null_rows": null_rows,
        "canonical_columns": np.asarray(canonical_columns, dtype=np.int64),
        "truth_columns": np.asarray(truth_columns, dtype=np.int64),
        "allowed": allowed,
        "truth_target_count": len(truth_targets),
        "truth_gene_count": len(truth_genes),
    }


def resolve_output_paths(approval: Mapping[str, Any]) -> tuple[Path, Path, Path]:
    output = approval["surface"]["output"]
    directory = RUN.resolve_repository_path(output["directory"], "evaluation output")
    receipt = RUN.resolve_repository_path(output["access_receipt"], "access receipt")
    result = RUN.resolve_repository_path(output["result"], "evaluation result")
    if (
        directory.relative_to(ROOT.resolve()).as_posix() != EXPECTED_OUTPUT_DIRECTORY
        or receipt != directory / "OUTCOME_ACCESS_RECEIPT.v1.json"
        or result != directory / "result.json"
    ):
        raise ContractError("evaluation output layout drifted")
    for path, label in (
        (directory, "evaluation output"),
        (receipt, "access receipt"),
        (result, "evaluation result"),
    ):
        RUN.require_under_experiment(path)
        RUN.reject_symlink_components(path, label)
    return directory, receipt, result


def audit_pretruth(*, require_output_absent: bool = True) -> dict[str, Any]:
    approval = load_approval()
    loaded, protocol = verify_frozen_chain(approval)
    frozen, arrays = loaded[EXPECTED_ROLE]
    mapping = validate_metadata_mapping(approval, frozen, arrays, protocol)
    output_directory, access_receipt, result_path = resolve_output_paths(approval)
    output_exists = (
        output_directory.exists()
        or output_directory.is_symlink()
        or access_receipt.exists()
        or access_receipt.is_symlink()
        or result_path.exists()
        or result_path.is_symlink()
    )
    if require_output_absent and output_exists:
        raise ContractError("held-evaluation output already exists")
    delta_path = _verify_file_record(
        approval["surface"]["truth"]["delta"],
        "HEK293T held truth delta",
        hash_contents=False,
    )
    control_path = _verify_file_record(
        approval["surface"]["truth"]["b_ctrl"],
        "HEK293T held truth control",
        hash_contents=False,
    )
    return {
        "schema_version": 1,
        "experiment": "114-depmap-coessentiality",
        "stage": "held_surface_evaluation_pretruth_audit",
        "approval_sha256": RUN.sha256_file(APPROVAL_PATH),
        "surface": EXPECTED_SURFACE_ID,
        "role": EXPECTED_ROLE,
        "diagnostic_only": True,
        "promotion_authorized": False,
        "metrics": ["raw_uncentered_pds"],
        "direction_valid": False,
        "changed_target_count": len(mapping["targets"]),
        "changed_target_ordered_sha256": RUN.sha256_lines(mapping["targets"]),
        "allowed_coordinate_count": int(mapping["allowed"].sum()),
        "all_frozen_bundles_validated": True,
        "null_mappings_regenerated": True,
        "truth_paths_fixed": True,
        "truth_files_present_with_expected_sizes": True,
        "truth_delta_path": delta_path.relative_to(ROOT.resolve()).as_posix(),
        "truth_control_path": control_path.relative_to(ROOT.resolve()).as_posix(),
        "truth_matrix_hashed": False,
        "truth_matrix_opened": False,
        "output_directory": output_directory.relative_to(ROOT.resolve()).as_posix(),
        "output_absent": not output_exists,
        "ready_for_independent_review": True,
        "_loaded": loaded,
        "_protocol": protocol,
        "_mapping": mapping,
        "_approval": approval,
    }


def _public_audit(report: Mapping[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in report.items() if not key.startswith("_")}


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def create_access_receipt(
    approval: Mapping[str, Any],
    output_directory: Path,
    access_receipt: Path,
    result_path: Path,
) -> None:
    parent = output_directory.parent
    RUN.require_under_experiment(parent)
    RUN.reject_symlink_components(parent, "evaluation output parent")
    if output_directory.exists() or output_directory.is_symlink():
        raise ContractError("held-evaluation output already exists")
    parent.mkdir(parents=True, exist_ok=True)
    RUN.reject_symlink_components(parent, "evaluation output parent")
    try:
        output_directory.mkdir(mode=0o700)
    except FileExistsError as error:
        raise ContractError("held-evaluation output already exists") from error
    _fsync_directory(parent)
    payload = {
        "schema_version": 1,
        "experiment": "114-depmap-coessentiality",
        "stage": "held_truth_access_started",
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "approval_sha256": RUN.sha256_file(APPROVAL_PATH),
        "entrypoint_sha256": RUN.sha256_file(Path(__file__)),
        "surface": EXPECTED_SURFACE_ID,
        "diagnostic_only": True,
        "promotion_authorized": False,
        "truth_delta": approval["surface"]["truth"]["delta"],
        "result_path": result_path.relative_to(ROOT.resolve()).as_posix(),
        "truth_hashed_before_receipt": False,
        "truth_opened_before_receipt": False,
    }
    encoded = (
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False)
        + "\n"
    ).encode()
    descriptor = os.open(
        access_receipt,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL,
        0o600,
    )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        access_receipt.unlink(missing_ok=True)
        raise
    _fsync_directory(output_directory)


def _verify_truth_bundle_after_receipt(
    approval: Mapping[str, Any], access_receipt: Path
) -> dict[str, Path]:
    if not access_receipt.is_file():
        raise ContractError("held truth cannot be accessed before the access receipt")
    paths: dict[str, Path] = {}
    for name, record in approval["surface"]["truth"].items():
        paths[name] = _verify_file_record(
            record, f"HEK293T held truth {name}", hash_contents=True
        )
    return paths


def _load_truth_effect(
    approval: Mapping[str, Any],
    paths: Mapping[str, Path],
    mapping: Mapping[str, Any],
    canonical_gene_count: int,
) -> np.ndarray:
    record = approval["surface"]["truth"]["delta"]
    delta = np.load(paths["delta"], mmap_mode="r", allow_pickle=False)
    try:
        if list(delta.shape) != record["shape"] or delta.dtype.str != record["dtype"]:
            raise ContractError("HEK293T held truth array shape or dtype drifted")
        selected = np.asarray(delta[mapping["truth_rows"]], dtype=np.float64)
        truth = np.zeros(
            (len(mapping["truth_rows"]), canonical_gene_count), dtype=np.float64
        )
        truth[:, mapping["canonical_columns"]] = selected[:, mapping["truth_columns"]]
    finally:
        mmap = getattr(delta, "_mmap", None)
        if mmap is not None:
            mmap.close()
    if not np.isfinite(truth[:, mapping["allowed"]]).all():
        raise ContractError("HEK293T held truth contains nonfinite scored values")
    return truth


def evaluate(*, confirmed: bool) -> dict[str, Any]:
    if not confirmed:
        raise ContractError("real evaluation requires --confirm-held-evaluation")
    audit = audit_pretruth(require_output_absent=True)
    approval = audit["_approval"]
    loaded = audit["_loaded"]
    mapping = audit["_mapping"]
    frozen, arrays = loaded[EXPECTED_ROLE]
    output_directory, access_receipt, result_path = resolve_output_paths(approval)
    create_access_receipt(
        approval, output_directory, access_receipt, result_path
    )

    paths = _verify_truth_bundle_after_receipt(approval, access_receipt)
    truth = _load_truth_effect(
        approval,
        paths,
        mapping,
        len(arrays["genes"]),
    )
    rows = mapping["frozen_rows"]
    allowed = mapping["allowed"]
    arms = {
        name: np.asarray(arrays[f"effect__{name}"])[rows] for name in RUN.ARMS
    }
    scores: dict[str, Any] = {}
    for name in RUN.ARMS:
        raw_pds, per_target = RUN.raw_uncentered_pds(
            arms[name], truth, allowed
        )
        scores[name] = {
            "raw_pds": raw_pds,
            "per_target_raw_pds": per_target.tolist(),
        }
    fallback_pds = float(scores["current_fallback"]["raw_pds"])
    null_deltas = np.empty(RUN.NULL_DRAWS, dtype=np.float64)
    null_rows = mapping["null_rows"]
    for draw in range(RUN.NULL_DRAWS):
        null_effect = np.asarray(arrays["null_effects"][draw, null_rows])
        null_pds, _ = RUN.raw_uncentered_pds(null_effect, truth, allowed)
        null_deltas[draw] = null_pds - fallback_pds

    result = {
        "schema_version": 1,
        "experiment": "114-depmap-coessentiality",
        "stage": "held_surface_pds_only_diagnostic",
        "surface": EXPECTED_SURFACE_ID,
        "role": EXPECTED_ROLE,
        "diagnostic_only": True,
        "promotion_authorized": False,
        "directional_evaluation_performed": False,
        "approval_sha256": RUN.sha256_file(APPROVAL_PATH),
        "entrypoint_sha256": RUN.sha256_file(Path(__file__)),
        "access_receipt_sha256": RUN.sha256_file(access_receipt),
        "frozen_manifest_sha256": RUN.sha256_file(
            RUN.resolve_repository_path(
                approval["frozen_bundles"][EXPECTED_ROLE]["manifest"]["path"],
                "HEK293T frozen manifest",
            )
        ),
        "truth_delta_sha256": approval["surface"]["truth"]["delta"]["sha256"],
        "target_count": len(mapping["targets"]),
        "target_ordered_sha256": RUN.sha256_lines(mapping["targets"]),
        "allowed_coordinate_count": int(allowed.sum()),
        "allowed_mask_sha256": RUN.sha256_array(allowed.astype(np.uint8)),
        "arms": scores,
        "primary_delta_raw_pds": (
            float(scores[RUN.PRIMARY_ARM]["raw_pds"]) - fallback_pds
        ),
        "permutation_raw_pds_deltas": null_deltas.tolist(),
        "permutation_p95": float(
            np.quantile(null_deltas, 0.95, method="higher")
        ),
        "no_post_truth_fit": True,
    }
    RUN.atomic_json(result_path, result)
    _fsync_directory(output_directory)
    return result


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    commands.add_parser(
        "audit",
        help="validate freezes and metadata mappings without opening held truth",
    )
    execute = commands.add_parser(
        "evaluate",
        help="run the fixed HEK293T PDS-only diagnostic",
    )
    execute.add_argument("--confirm-held-evaluation", action="store_true")
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.command == "audit":
        print(json.dumps(_public_audit(audit_pretruth()), indent=2, sort_keys=True))
        return 0
    print(
        json.dumps(
            evaluate(confirmed=args.confirm_held_evaluation),
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
