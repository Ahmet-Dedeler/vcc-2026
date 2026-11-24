#!/usr/bin/env python3
"""Audit and execute Experiment 114's versioned source-only freeze approval.

This entrypoint leaves the sealed predecessor protocol, boundary, implementation,
role-source contract, and role bundles byte-identical.  It has no held-truth or
evaluation interface.
"""

from __future__ import annotations

import argparse
import copy
import gc
import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType
from typing import Any, Mapping

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUN_PATH = HERE / "run.py"
BUILDER_PATH = HERE / "build_role_inputs.py"
APPROVAL_PATH = HERE / "SOURCE_CANDIDATE_FREEZE_APPROVAL.v1.json"
EXPECTED_APPROVAL_SHA256 = (
    "fa360e7b04c61d5effb137de3d990925b899f9f15cf0e2551459fb28752113a4"
)
EXPECTED_PREDECESSORS = {
    "preregistered_protocol": {
        "path": "experiments/114-depmap-coessentiality/preregistered_protocol.json",
        "sha256": (
            "cfd495175d87cba7e7cad1ff171b7a817664eabc3e836aa311ec7417f0c59b91"
        ),
    },
    "source_candidate_freeze_boundary": {
        "path": (
            "experiments/114-depmap-coessentiality/"
            "SOURCE_CANDIDATE_FREEZE_BOUNDARY.json"
        ),
        "sha256": (
            "a22ca21f51e7d6729c9c7f09c2901ba0fa9cf1ecb067fc75ab91d9324ef18025"
        ),
    },
    "run_implementation": {
        "path": "experiments/114-depmap-coessentiality/run.py",
        "sha256": (
            "737d0bf225c5ee5e5108a54a07e9fdfc3ee6e8728f1ad559def5ff2c03569fd9"
        ),
    },
    "role_source_manifest": {
        "path": (
            "experiments/114-depmap-coessentiality/ROLE_SOURCE_MANIFEST.json"
        ),
        "sha256": (
            "e05d877b5ea0a7aeae06f4c53ad363e984379678b5b61fd852e3e877d17d0248"
        ),
    },
}
EXPECTED_ROLE_MANIFEST_SHA256 = {
    "h1": "357244f5515db1e51686facc1d923da6209fe9ac95c3ce686c28b7474bece017",
    "orion_hek293t": (
        "939c43916c3addb0a2fd2ff7753a8c69dd866d39e33a0481d901f0dc87247ad3"
    ),
}
EXPECTED_OUTPUTS = {
    "h1": (
        "experiments/114-depmap-coessentiality/"
        "work/source-candidate-freeze-v1/h1"
    ),
    "orion_hek293t": (
        "experiments/114-depmap-coessentiality/"
        "work/source-candidate-freeze-v1/orion_hek293t"
    ),
}
EXPECTED_AUTHORIZATION = {
    "source_candidate_freeze": True,
    "held_truth_access": False,
    "model_fitting": False,
    "evaluation": False,
    "count_generation": False,
    "challenge_build": False,
    "submission": False,
    "commit": False,
    "push": False,
}


def _load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


RUN = _load_module("experiment_114_freeze_v1_run", RUN_PATH)
BUILDER = _load_module("experiment_114_freeze_v1_builder", BUILDER_PATH)
ContractError = RUN.ContractError


def load_approval() -> dict[str, Any]:
    if RUN.sha256_file(APPROVAL_PATH) != EXPECTED_APPROVAL_SHA256:
        raise ContractError("source-candidate approval amendment changed")
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
        "authorization",
        "approved_roles",
        "invariants",
    }
    if set(approval) != expected_keys:
        raise ContractError("approval amendment schema drifted")
    expected_identity = {
        "schema_version": 1,
        "experiment": "114-depmap-coessentiality",
        "amendment_id": "experiment-114-source-candidate-freeze-approval-v1",
        "date_approved": "2026-08-30",
        "stage": "source_candidate_freeze_approval",
        "status": (
            "approved_for_source_only_freeze_pending_independent_review_and_execution"
        ),
        "entrypoint_path": (
            "experiments/114-depmap-coessentiality/freeze_source_candidates_v1.py"
        ),
    }
    for key, expected in expected_identity.items():
        if approval.get(key) != expected:
            raise ContractError(f"approval amendment identity drifted: {key}")
    if approval.get("predecessors") != EXPECTED_PREDECESSORS:
        raise ContractError("approval predecessor chain drifted")
    if approval.get("authorization") != EXPECTED_AUTHORIZATION:
        raise ContractError("approval authorizes more than source_candidate_freeze")
    invariants = approval.get("invariants")
    required_invariants = {
        "allowed_action": "source_candidate_freeze",
        "canonical_axes_required": True,
        "exact_a14_parameters_required": True,
        "exact_eligibility_required": True,
        "full_artifact_receipt_and_transitive_hash_verification_required": True,
        "held_surface_coordinates_zero_and_uncovered_required": True,
        "held_truth_must_not_be_opened": True,
        "operator_mask_must_equal_source_bank_coverage_required": True,
        "outputs_must_be_absent_and_must_not_be_overwritten": True,
        "panel_coordinates_zero_and_uncovered_required": True,
        "real_freeze_executed_by_this_amendment": False,
    }
    if invariants != required_invariants:
        raise ContractError("approval invariant set drifted")
    roles = approval.get("approved_roles")
    if not isinstance(roles, dict) or tuple(roles) != tuple(
        EXPECTED_ROLE_MANIFEST_SHA256
    ):
        raise ContractError("approval role order or names drifted")
    for role, record in roles.items():
        if record.get("manifest_sha256") != EXPECTED_ROLE_MANIFEST_SHA256[role]:
            raise ContractError(f"{role}: approved role-manifest hash drifted")
        if record.get("output_directory") != EXPECTED_OUTPUTS[role]:
            raise ContractError(f"{role}: approved freeze output path drifted")
        RUN.approved_role_for_manifest(
            RUN.resolve_repository_path(record["manifest_path"], f"{role} manifest"),
            {"approved_roles": {role: record}},
        )


def verify_predecessors(approval: Mapping[str, Any]) -> None:
    for name, record in approval["predecessors"].items():
        RUN.verify_fixed_record(record, f"approval predecessor {name}")
    protocol = RUN.load_json(RUN.PROTOCOL_PATH)
    source_manifest = RUN.load_json(RUN.SOURCE_MANIFEST_PATH)
    boundary = RUN.load_json(RUN.SOURCE_CANDIDATE_FREEZE_BOUNDARY_PATH)
    RUN.validate_protocol(protocol)
    RUN.validate_source_manifest(source_manifest)
    RUN.validate_source_candidate_freeze_boundary(boundary)
    BUILDER.validate_checked_in_locks()
    BUILDER.load_role_source_manifest()


def effective_authorization(
    approval: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    protocol = copy.deepcopy(RUN.load_json(RUN.PROTOCOL_PATH))
    protocol["current_authorization"] = {
        key: key
        in {
            "small_public_code_and_metadata_inspection",
            "local_synthetic_tests",
            "source_candidate_freeze",
        }
        for key in RUN.EXPECTED_AUTHORIZATION
    }
    boundary = copy.deepcopy(RUN.load_json(RUN.SOURCE_CANDIDATE_FREEZE_BOUNDARY_PATH))
    boundary["source_candidate_freeze_authorized"] = True
    boundary["approved_roles"] = copy.deepcopy(approval["approved_roles"])
    return protocol, boundary


def _load_role_receipts(
    role_directory: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    incumbent = RUN.load_json(role_directory / "incumbent-receipt.json")
    source = RUN.load_json(role_directory / "source-bank-receipt.json")
    return incumbent, source


def _verify_exact_transitive_inputs(
    *,
    role: str,
    role_spec: Mapping[str, Any],
    operator: Any,
    role_source_manifest: Mapping[str, Any],
    incumbent_receipt: Mapping[str, Any],
    source_receipt: Mapping[str, Any],
) -> None:
    hash_cache = {
        operator.artifact_path.resolve(): operator.artifact_sha256,
        operator.manifest_path.resolve(): operator.manifest_sha256,
    }
    expected_incumbent = BUILDER._source_records_for_role(
        role_source_manifest,
        role_spec,
        include_incumbent_sources=True,
        operator=operator,
        hash_cache=hash_cache,
    )
    expected_source = BUILDER._source_records_for_role(
        role_source_manifest,
        role_spec,
        include_incumbent_sources=False,
        operator=operator,
        hash_cache=hash_cache,
    )
    if incumbent_receipt.get("source_inputs") != expected_incumbent:
        raise ContractError(f"{role}: incumbent transitive source set drifted")
    if source_receipt.get("source_inputs") != expected_source:
        raise ContractError(f"{role}: source-bank transitive source set drifted")


def validate_operator_mask_and_zeroing(
    *,
    role: str,
    incumbent: Mapping[str, np.ndarray],
    source_bank: Any,
    canonical_targets: np.ndarray,
    canonical_genes: np.ndarray,
    operator: Any,
    source_receipt: Mapping[str, Any],
    chunk_rows: int = 128,
) -> dict[str, int]:
    columns = np.asarray(operator.columns, dtype=np.int64)
    if len(columns) == 0 or not np.array_equal(columns, np.unique(columns)):
        raise ContractError(f"{role}: operator columns are invalid")
    panel_positions = np.flatnonzero(np.isin(canonical_genes, canonical_targets))
    held_positions = np.flatnonzero(
        np.isin(canonical_genes, np.asarray(operator.held_targets, dtype=str))
    )
    excluded = np.union1d(panel_positions, held_positions)
    if np.intersect1d(columns, excluded).size:
        raise ContractError(f"{role}: operator retains panel or held coordinates")
    if np.any(np.asarray(incumbent["effects"])[:, panel_positions] != 0):
        raise ContractError(f"{role}: incumbent panel coordinates are nonzero")
    if np.any(np.asarray(incumbent["effects"])[:, held_positions] != 0):
        raise ContractError(f"{role}: incumbent held coordinates are nonzero")
    for start in range(0, len(source_bank.targets), chunk_rows):
        stop = min(start + chunk_rows, len(source_bank.targets))
        covered = source_bank.covered[start:stop]
        effects = source_bank.effects[start:stop]
        if not np.all(covered[:, columns]):
            raise ContractError(f"{role}: operator coordinates are not fully covered")
        if np.count_nonzero(covered) != len(covered) * len(columns):
            raise ContractError(f"{role}: coverage differs from the operator mask")
        if np.any(effects[:, panel_positions] != 0):
            raise ContractError(f"{role}: source-bank panel coordinates are nonzero")
        if np.any(covered[:, panel_positions]):
            raise ContractError(f"{role}: source-bank panel coordinates are covered")
        if np.any(effects[:, held_positions] != 0):
            raise ContractError(f"{role}: source-bank held coordinates are nonzero")
        if np.any(covered[:, held_positions]):
            raise ContractError(f"{role}: source-bank held coordinates are covered")
    artifact_audit = source_receipt.get("parameters", {}).get("artifact", {})
    if artifact_audit.get("covered_coordinates_per_row") != len(columns):
        raise ContractError(f"{role}: receipt operator coverage count drifted")
    if artifact_audit.get("panel_coordinates_zero_and_uncovered") is not True:
        raise ContractError(f"{role}: receipt panel-zero assertion drifted")
    if artifact_audit.get("held_surface_coordinates_zero_and_uncovered") is not True:
        raise ContractError(f"{role}: receipt held-zero assertion drifted")
    return {
        "operator_column_count": len(columns),
        "panel_coordinate_count": len(panel_positions),
        "held_coordinate_count_on_gene_axis": len(held_positions),
    }


def audit_role(
    role: str,
    approval: Mapping[str, Any],
    *,
    require_output_absent: bool = True,
) -> dict[str, Any]:
    effective_protocol, effective_boundary = effective_authorization(approval)
    source_manifest = RUN.load_json(RUN.SOURCE_MANIFEST_PATH)
    role_source_manifest = BUILDER.load_role_source_manifest()
    approved = approval["approved_roles"][role]
    role_spec = role_source_manifest["roles"][role]
    manifest_path = RUN.resolve_repository_path(
        approved["manifest_path"], f"{role} role manifest"
    )
    role_name, selected = RUN.approved_role_for_manifest(
        manifest_path, effective_boundary
    )
    if role_name != role or selected != approved:
        raise ContractError(f"{role}: effective approval selection drifted")
    role_manifest = RUN.load_json(manifest_path)
    RUN.validate_input_manifest(role_manifest, role=role, approved=approved)
    role_directory = RUN.validate_role_manifest(
        role_manifest,
        role=role,
        approved=approved,
    )
    expected_files = {
        "incumbent.npz",
        "incumbent-receipt.json",
        "source-bank.npz",
        "source-bank-receipt.json",
        "manifest.json",
    }
    if {path.name for path in role_directory.iterdir()} != expected_files:
        raise ContractError(f"{role}: sealed role-directory contents drifted")
    output = RUN.resolve_repository_path(
        approved["output_directory"], f"{role} freeze output"
    )
    RUN.require_under_experiment(output)
    RUN.reject_symlink_components(output, f"{role} freeze output")
    if require_output_absent and (output.exists() or output.is_symlink()):
        raise ContractError(f"{role}: freeze output already exists")

    canonical_targets, canonical_genes = RUN.validate_canonical_axes(
        effective_boundary, source_manifest
    )
    sealed_protocol = RUN.load_json(RUN.PROTOCOL_PATH)
    RUN.validate_depmap_feature_artifact(
        effective_boundary, sealed_protocol, source_manifest
    )
    RUN.validate_string_feature_artifact(effective_boundary, source_manifest)
    incumbent, source_bank, role_audit = RUN.load_validated_role_artifacts(
        role_manifest,
        role=role,
        approved=approved,
        role_directory=role_directory,
        canonical_targets=canonical_targets,
        canonical_genes=canonical_genes,
        protocol=effective_protocol,
        source_manifest=source_manifest,
    )
    incumbent_receipt, source_receipt = _load_role_receipts(role_directory)
    if incumbent_receipt.get("parameters", {}).get("a14_parameters") != (
        role_source_manifest["a14"]
    ):
        raise ContractError(f"{role}: incumbent exact a14 parameters drifted")
    if source_receipt.get("parameters", {}).get("a14") != role_source_manifest["a14"]:
        raise ContractError(f"{role}: source-bank exact a14 parameters drifted")
    if source_receipt.get("parameters", {}).get("source_bank") != (
        role_source_manifest["source_bank"]
    ):
        raise ContractError(f"{role}: source-bank construction contract drifted")
    if incumbent_receipt.get("parameters", {}).get("role_donor_order") != (
        role_spec["allowed_source_ids"]
    ):
        raise ContractError(f"{role}: incumbent donor order drifted")

    operator = BUILDER.load_role_operator(role, role_spec)
    operator_receipt = source_receipt.get("parameters", {}).get("operator")
    expected_operator_receipt = {
        "artifact_path": operator.artifact_path.relative_to(ROOT.resolve()).as_posix(),
        "artifact_sha256": operator.artifact_sha256,
        "manifest_path": operator.manifest_path.relative_to(ROOT.resolve()).as_posix(),
        "manifest_sha256": operator.manifest_sha256,
        "columns_int64_sha256": BUILDER.sha256_int64(operator.columns),
        "held_target_ordered_sha256": RUN.sha256_lines(operator.held_targets),
    }
    if operator_receipt != expected_operator_receipt:
        raise ContractError(f"{role}: Experiment 107 operator receipt drifted")
    _verify_exact_transitive_inputs(
        role=role,
        role_spec=role_spec,
        operator=operator,
        role_source_manifest=role_source_manifest,
        incumbent_receipt=incumbent_receipt,
        source_receipt=source_receipt,
    )
    mask_audit = validate_operator_mask_and_zeroing(
        role=role,
        incumbent=incumbent,
        source_bank=source_bank,
        canonical_targets=canonical_targets,
        canonical_genes=canonical_genes,
        operator=operator,
        source_receipt=source_receipt,
    )
    result = {
        "role": role,
        "manifest_sha256": RUN.sha256_file(manifest_path),
        "role_directory": role_directory.relative_to(ROOT.resolve()).as_posix(),
        "output_directory": output.relative_to(ROOT.resolve()).as_posix(),
        "output_absent": not output.exists() and not output.is_symlink(),
        "incumbent_artifact_sha256": role_manifest["incumbent"]["artifact"]["sha256"],
        "incumbent_receipt_sha256": role_audit["incumbent_receipt_sha256"],
        "source_bank_artifact_sha256": role_manifest["source_bank"]["artifact"][
            "sha256"
        ],
        "source_bank_receipt_sha256": role_audit["source_bank_receipt_sha256"],
        "allowed_source_ids": approved["allowed_source_ids"],
        "forbidden_source_ids": approved["forbidden_source_ids"],
        "held_surface_source_id": approved["held_surface_source_id"],
        "eligibility": role_audit["eligibility"],
        "mask_audit": mask_audit,
        "held_truth_opened": False,
    }
    del incumbent, source_bank, operator
    gc.collect()
    return result


def audit_all() -> dict[str, Any]:
    approval = load_approval()
    verify_predecessors(approval)
    roles = {
        role: audit_role(role, approval)
        for role in EXPECTED_ROLE_MANIFEST_SHA256
    }
    return {
        "experiment": "114-depmap-coessentiality",
        "stage": "source_candidate_freeze_preflight",
        "approval_sha256": RUN.sha256_file(APPROVAL_PATH),
        "predecessor_hashes_verified": True,
        "authorized_action": "source_candidate_freeze",
        "held_truth_access_authorized": False,
        "evaluation_authorized": False,
        "count_generation_authorized": False,
        "submission_authorized": False,
        "roles": roles,
        "ready_for_independent_review": True,
        "real_candidates_frozen": False,
    }


def freeze_role(role: str, *, confirmed: bool) -> dict[str, Any]:
    if not confirmed:
        raise ContractError("real freeze requires --confirm-source-only-freeze")
    approval = load_approval()
    verify_predecessors(approval)
    audit = audit_role(role, approval)
    effective_protocol, effective_boundary = effective_authorization(approval)
    approved = approval["approved_roles"][role]
    manifest_path = RUN.resolve_repository_path(
        approved["manifest_path"], f"{role} role manifest"
    )
    output = RUN.resolve_repository_path(
        approved["output_directory"], f"{role} freeze output"
    )
    original_loader = RUN.load_freeze_boundary
    original_depmap_validator = RUN.validate_depmap_feature_artifact
    sealed_protocol = RUN.load_json(RUN.PROTOCOL_PATH)

    def validate_prebuilt_depmap(
        boundary: Mapping[str, Any],
        _effective_protocol: Mapping[str, Any],
        source_manifest: Mapping[str, Any],
    ) -> Any:
        return original_depmap_validator(
            boundary,
            sealed_protocol,
            source_manifest,
        )

    try:
        RUN.load_freeze_boundary = lambda: copy.deepcopy(effective_boundary)
        RUN.validate_depmap_feature_artifact = validate_prebuilt_depmap
        frozen = RUN.freeze_from_manifest(
            manifest_path,
            output,
            protocol=effective_protocol,
        )
    finally:
        RUN.load_freeze_boundary = original_loader
        RUN.validate_depmap_feature_artifact = original_depmap_validator
    return {
        "approval_sha256": RUN.sha256_file(APPROVAL_PATH),
        "preflight": audit,
        "frozen_manifest": frozen,
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    commands.add_parser(
        "audit",
        help="read-only audit of both approved source-only role bundles",
    )
    freeze = commands.add_parser(
        "freeze",
        help="freeze one source-only role; no held-truth or evaluation interface exists",
    )
    freeze.add_argument("--role", choices=tuple(EXPECTED_ROLE_MANIFEST_SHA256), required=True)
    freeze.add_argument("--confirm-source-only-freeze", action="store_true")
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.command == "audit":
        print(json.dumps(audit_all(), indent=2, sort_keys=True))
        return 0
    result = freeze_role(args.role, confirmed=args.confirm_source_only_freeze)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
