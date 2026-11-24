#!/usr/bin/env python3
"""Source-only DepMap co-essentiality target-kernel scaffold.

The checked-in protocol authorizes only the pinned source download and
deterministic feature construction. Candidate freezing and held-out scoring
remain closed behind separate protocol amendments.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
SOURCE_MANIFEST_PATH = HERE / "SOURCE_MANIFEST.json"
SOURCE_SCHEMA_AMENDMENT_PATH = HERE / "SOURCE_SCHEMA_AMENDMENT.json"
SOURCE_CANDIDATE_FREEZE_BOUNDARY_PATH = HERE / "SOURCE_CANDIDATE_FREEZE_BOUNDARY.json"
ROLE_SOURCE_MANIFEST_PATH = HERE / "ROLE_SOURCE_MANIFEST.json"
SOURCE_FEATURE_DIRECTORY = HERE / "work" / "source-features"
DEPMAP_CSV_PATH = SOURCE_FEATURE_DIRECTORY / "CRISPRGeneEffect.csv"
DEPMAP_SOURCE_RECEIPT_PATH = SOURCE_FEATURE_DIRECTORY / "source-receipt.json"
DEPMAP_FEATURES_PATH = SOURCE_FEATURE_DIRECTORY / "depmap-features.npz"
DEPMAP_FEATURES_RECEIPT_PATH = SOURCE_FEATURE_DIRECTORY / "feature-receipt.json"
STRING_FEATURES_PATH = (
    ROOT / "data" / "external" / "slim" / "gene_string_embeddings.v0.3.h5"
)
PANEL_PATH = ROOT / "data" / "challenge" / "pert_counts.csv"
GENE_NAMES_PATH = ROOT / "data" / "challenge" / "gene_names.csv"
EXPECTED_PROTOCOL_SHA256 = (
    "cfd495175d87cba7e7cad1ff171b7a817664eabc3e836aa311ec7417f0c59b91"
)
FEATURE_BUILD_PROTOCOL_SHA256 = (
    "90c3c458f182ce4fa39a62757fdffac69ff105c26fa57bd6e28cec2e5804edc0"
)
EXPECTED_SOURCE_MANIFEST_SHA256 = (
    "5834bcd29a7e31b23aaafdcd02a5f56f509e06c87ad1fb0c2568df40510a5497"
)
EXPECTED_SOURCE_SCHEMA_AMENDMENT_SHA256 = (
    "266f805bc80d011b53322bbf81c54bd69507e39c3752fa9a47c8b95cae886075"
)
EXPECTED_SOURCE_CANDIDATE_FREEZE_BOUNDARY_SHA256 = (
    "a22ca21f51e7d6729c9c7f09c2901ba0fa9cf1ecb067fc75ab91d9324ef18025"
)
EXPECTED_ROLE_SOURCE_MANIFEST_SHA256 = (
    "e05d877b5ea0a7aeae06f4c53ad363e984379678b5b61fd852e3e877d17d0248"
)
ACQUISITION_PROTOCOL_SHA256 = (
    "c1b7241a1423affe334c4bd52b2c0c730ade6925c0f7249b5e86eb0da2feac27"
)
EXPECTED_SOURCE_RECEIPT_SHA256 = (
    "ff33216dd9997026e0cfb756b2139f80724ed4d03fbaa75b283048e488e8330c"
)
EXPECTED_FEATURE_RECEIPT_SHA256 = (
    "6e01f6daa6bb39203d8178ea693e39caad6246fbe96fa38dbc18d25deadd7b7b"
)
FEATURE_BUILD_IMPLEMENTATION_SHA256 = (
    "e32cfd9f6ceecb140bb02214f74c40f6e1fde796e9be2bfe812423931081ce67"
)
EXPECTED_DEPMAP_SHAPE = (1_178, 17_916)
ARMS = (
    "identity",
    "current_fallback",
    "string_neighbor",
    "depmap_neighbor",
    "string_depmap_blend_w1",
)
PRIMARY_ARM = "string_depmap_blend_w1"
NULL_DRAWS = 100
NEIGHBOURS = 25
NEIGHBOUR_POWER = 3.0
DIRECTIONAL_FRACTION = 0.20
EXPECTED_AUTHORIZATION = {
    "small_public_code_and_metadata_inspection": True,
    "local_synthetic_tests": True,
    "depmap_large_download": True,
    "source_feature_build": True,
    "source_candidate_freeze": False,
    "held_truth_access": False,
    "model_fitting": False,
    "evaluation": False,
    "count_generation": False,
    "challenge_build": False,
    "submission": False,
    "commit": False,
    "push": False,
}


class ContractError(RuntimeError):
    """A frozen source, eligibility, leakage, or output contract failed."""


@dataclass(frozen=True)
class GeneFeatures:
    genes: np.ndarray
    matrix: np.ndarray

    def __post_init__(self) -> None:
        genes = np.asarray(self.genes).astype(str)
        matrix = np.asarray(self.matrix, dtype=np.float32)
        if genes.ndim != 1 or len(genes) == 0:
            raise ContractError("feature gene axis must be a nonempty vector")
        if len(set(genes.tolist())) != len(genes):
            raise ContractError("feature gene axis contains duplicates")
        if matrix.ndim != 2 or matrix.shape[0] != len(genes):
            raise ContractError("feature matrix does not align with its gene axis")
        if not np.isfinite(matrix).all():
            raise ContractError("feature matrix contains nonfinite values")
        object.__setattr__(self, "genes", genes)
        object.__setattr__(self, "matrix", matrix)


@dataclass(frozen=True)
class SourceBank:
    targets: np.ndarray
    genes: np.ndarray
    effects: np.ndarray
    covered: np.ndarray
    source_ids: np.ndarray


@dataclass(frozen=True)
class ConsensusBank:
    targets: np.ndarray
    genes: np.ndarray
    effects: np.ndarray
    covered: np.ndarray
    carrier_count: np.ndarray


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def md5_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return sha256_bytes(("\n".join(map(str, values)) + "\n").encode())


def sha256_array(values: np.ndarray) -> str:
    array = np.ascontiguousarray(values)
    header = json.dumps(
        {"dtype": array.dtype.str, "shape": list(array.shape)},
        separators=(",", ":"),
        sort_keys=True,
    ).encode()
    return sha256_bytes(header + b"\0" + array.view(np.uint8).tobytes())


def stable_seed(*parts: object) -> int:
    return int.from_bytes(
        hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:8],
        "little",
    )


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _reject_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON constant: {value}")


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
            parse_constant=_reject_constant,
        )
    except (OSError, json.JSONDecodeError, ValueError) as error:
        raise ContractError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise ContractError(f"{path} must contain a JSON object")
    return value


def require_under_experiment(path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_relative_to(HERE.resolve()):
        raise ContractError(f"output escapes experiment directory: {resolved}")
    return resolved


def reject_symlink_components(path: Path, label: str) -> None:
    absolute = path.expanduser().absolute()
    current = Path(absolute.anchor)
    for part in absolute.parts[1:]:
        current /= part
        if current.is_symlink():
            raise ContractError(f"{label} contains a symlink component: {current}")


def resolve_repository_path(path_value: str, label: str) -> Path:
    if not isinstance(path_value, str) or not path_value:
        raise ContractError(f"{label} path must be a nonempty string")
    relative = Path(path_value)
    if relative.is_absolute() or ".." in relative.parts:
        raise ContractError(f"{label} path must be repository-relative")
    candidate = ROOT / relative
    reject_symlink_components(candidate, label)
    resolved = candidate.resolve()
    if not resolved.is_relative_to(ROOT.resolve()):
        raise ContractError(f"{label} escapes the repository: {resolved}")
    return resolved


def verify_fixed_record(
    record: Mapping[str, Any],
    label: str,
    *,
    required_under: Path | None = None,
) -> Path:
    path_value = record.get("path")
    expected_hash = record.get("sha256")
    if not isinstance(path_value, str) or not isinstance(expected_hash, str):
        raise ContractError(f"{label} requires repository-relative path and sha256")
    path = resolve_repository_path(path_value, label)
    if required_under is not None and not path.is_relative_to(required_under.resolve()):
        raise ContractError(f"{label} is outside its frozen role directory")
    if not path.is_file():
        raise ContractError(f"{label} is missing: {path}")
    if sha256_file(path) != expected_hash:
        raise ContractError(f"{label} hash drifted")
    return path


def atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    path = require_under_experiment(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.is_symlink():
        raise ContractError(f"refusing to overwrite JSON artifact: {path}")
    descriptor, temporary = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, path)
        Path(temporary).unlink()
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def atomic_npz(path: Path, **values: object) -> None:
    path = require_under_experiment(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.is_symlink():
        raise ContractError(f"refusing to overwrite NPZ artifact: {path}")
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("wb") as handle:
            np.savez_compressed(handle, **values)
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, path)
        temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ContractError(f"{path}: missing CSV column {column!r}")
        values = tuple(
            str(row[column]).strip()
            for row in reader
            if row.get(column) and str(row[column]).strip()
        )
    if not values or len(values) != len(set(values)):
        raise ContractError(f"{path}: {column!r} is empty or contains duplicates")
    return values


def read_lines(path: Path) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if not values or len(values) != len(set(values)):
        raise ContractError(f"{path}: expected nonempty unique lines")
    return values


def row_normalize(matrix: np.ndarray) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float32)
    if values.ndim != 2:
        raise ContractError("row_normalize expects a matrix")
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return values / norms


def validate_protocol(protocol: Mapping[str, Any]) -> None:
    if sha256_file(PROTOCOL_PATH) != EXPECTED_PROTOCOL_SHA256:
        raise ContractError("preregistered protocol changed")
    if protocol.get("experiment") != "114-depmap-coessentiality":
        raise ContractError("wrong experiment protocol")
    if (
        protocol.get("status")
        != "role_source_builders_implemented_real_bundles_not_built_or_authorized"
    ):
        raise ContractError("protocol status drifted")
    if protocol.get("primary_arm") != PRIMARY_ARM:
        raise ContractError("primary arm drifted")
    incumbent = protocol.get("incumbent_contract")
    if not isinstance(incumbent, dict):
        raise ContractError("incumbent support contract is missing")
    expected_incumbent_provenance = {
        "name": "a14_support_contract",
        "reference_metadata_role": (
            "audit_only_cross_check_of_pre_cd4_perturbmatch_support_fields_shared_by_a14_and_a15"
        ),
        "reference_metadata_fields_used": [
            "effects.perturbmatch.neither_supported_preserved"
        ],
        "eligibility_primary_source": (
            "hash_pinned_panel_k562_perturbmatch_and_orion_target_lists"
        ),
        "alpha": 0.982,
        "donor_weights": "equal_over_role_donors_that_carry_the_target",
        "gene_norm": False,
        "empirical_bayes": False,
        "beta": 0.0,
        "perturbmatch_blend": 0.5,
        "whitener_id": "union_pruned",
        "whiten_gamma": 0.65,
        "whiten_lambda": 0.0,
        "basis_prune": "common_essential_after_sampling",
        "basis_row_cap_per_donor": 12_000,
        "basis_rng": "numpy_default_rng_seed_0_one_stream_in_role_donor_order",
        "cd4_fraction": 0.05,
        "cd4_target_solve": None,
    }
    for key, value in expected_incumbent_provenance.items():
        if incumbent.get(key) != value:
            raise ContractError(f"incumbent support provenance drifted: {key}")
    arms = protocol.get("arms")
    if not isinstance(arms, dict) or tuple(arms) != ARMS:
        raise ContractError(f"arm order or names drifted: {tuple(arms or {})}")
    eligibility = protocol.get("eligibility")
    if not isinstance(eligibility, dict):
        raise ContractError("missing eligibility contract")
    for prefix, count in (
        ("incumbent_uncovered", 28),
        ("low_confidence", 13),
        ("eligible", 41),
    ):
        values = eligibility.get(f"{prefix}_targets")
        if not isinstance(values, list) or len(values) != count:
            raise ContractError(f"{prefix} target count drifted")
        if sha256_lines(values) != eligibility.get(f"{prefix}_ordered_sha256"):
            raise ContractError(f"{prefix} ordered hash drifted")
    fallback = set(eligibility["incumbent_uncovered_targets"])
    low = set(eligibility["low_confidence_targets"])
    eligible = set(eligibility["eligible_targets"])
    if fallback & low or fallback | low != eligible:
        raise ContractError("eligibility strata are not disjoint and exhaustive")
    construction = protocol.get("candidate_construction", {})
    if construction.get("neighbours") != NEIGHBOURS:
        raise ContractError("neighbor count drifted")
    if construction.get("neighbour_tie_policy") != (
        "exactly 25 positive neighbours maximum; descending similarity with ties "
        "resolved by frozen sorted source-target order"
    ):
        raise ContractError("neighbor tie policy drifted")
    if construction.get("neighbour_weight") != "max(cosine,0)^3":
        raise ContractError("neighbor power drifted")
    if construction.get("confidence") != (
        "absolute value of the final stored effect array after dtype conversion"
    ):
        raise ContractError("confidence contract drifted")
    nulls = protocol.get("evaluation", {}).get("target_label_permutation_nulls", {})
    if nulls.get("draws") != NULL_DRAWS:
        raise ContractError("permutation draw count drifted")
    auth = protocol.get("current_authorization")
    if auth != EXPECTED_AUTHORIZATION:
        raise ContractError("authorization block drifted from the exact reviewed map")
    amendment_link = protocol.get("source_schema_amendment")
    expected_amendment_link = {
        "path": SOURCE_SCHEMA_AMENDMENT_PATH.relative_to(ROOT).as_posix(),
        "sha256": EXPECTED_SOURCE_SCHEMA_AMENDMENT_SHA256,
        "timing": (
            "after_verified_source_download_before_feature_construction_or_held_truth_access"
        ),
        "amended_gene_columns": EXPECTED_DEPMAP_SHAPE[1],
    }
    if amendment_link != expected_amendment_link:
        raise ContractError("source schema amendment link drifted")
    validate_source_schema_amendment(load_json(SOURCE_SCHEMA_AMENDMENT_PATH))
    freeze_boundary_link = protocol.get("source_candidate_freeze_boundary")
    expected_freeze_boundary_link = {
        "path": SOURCE_CANDIDATE_FREEZE_BOUNDARY_PATH.relative_to(ROOT).as_posix(),
        "sha256": EXPECTED_SOURCE_CANDIDATE_FREEZE_BOUNDARY_SHA256,
        "feature_build_protocol_sha256": FEATURE_BUILD_PROTOCOL_SHA256,
        "authorization": (
            "false_until_each_role_manifest_and_artifact_bundle_is_hash_pinned_and_"
            "independently_reviewed"
        ),
    }
    if freeze_boundary_link != expected_freeze_boundary_link:
        raise ContractError("source-candidate freeze boundary link drifted")
    role_source_link = protocol.get("role_source_manifest")
    expected_role_source_link = {
        "path": ROLE_SOURCE_MANIFEST_PATH.relative_to(ROOT).as_posix(),
        "sha256": EXPECTED_ROLE_SOURCE_MANIFEST_SHA256,
        "status": "builders_defined_real_role_bundles_not_built_or_approved",
        "roles": ["h1", "orion_hek293t"],
    }
    if role_source_link != expected_role_source_link:
        raise ContractError("role source manifest link drifted")
    validate_role_source_manifest(load_json(ROLE_SOURCE_MANIFEST_PATH))
    freeze_boundary = load_json(SOURCE_CANDIDATE_FREEZE_BOUNDARY_PATH)
    validate_source_candidate_freeze_boundary(freeze_boundary)
    for prefix in ("incumbent_uncovered", "low_confidence", "eligible"):
        boundary_record = freeze_boundary["eligibility"]
        if boundary_record[f"{prefix}_count"] != len(eligibility[f"{prefix}_targets"]):
            raise ContractError(f"freeze-boundary {prefix} count drifted")
        if (
            boundary_record[f"{prefix}_ordered_sha256"]
            != eligibility[f"{prefix}_ordered_sha256"]
        ):
            raise ContractError(f"freeze-boundary {prefix} hash drifted")


def validate_source_schema_amendment(amendment: Mapping[str, Any]) -> None:
    if (
        sha256_file(SOURCE_SCHEMA_AMENDMENT_PATH)
        != EXPECTED_SOURCE_SCHEMA_AMENDMENT_SHA256
    ):
        raise ContractError("source schema amendment changed")
    expected = {
        "schema_version": 1,
        "experiment": "114-depmap-coessentiality",
        "date_amended": "2026-08-30",
        "status": (
            "frozen_after_source_download_before_feature_construction_or_held_truth_access"
        ),
        "acquisition_protocol_sha256": ACQUISITION_PROTOCOL_SHA256,
        "rows_cell_lines": EXPECTED_DEPMAP_SHAPE[0],
        "held_truth_opened": False,
        "candidate_construction_started": False,
        "feature_construction_started": False,
    }
    for key, value in expected.items():
        if amendment.get(key) != value:
            raise ContractError(f"source schema amendment drifted: {key}")
    manifest = amendment.get("source_manifest")
    if not isinstance(manifest, dict) or manifest != {
        "path": SOURCE_MANIFEST_PATH.relative_to(ROOT).as_posix(),
        "sha256": EXPECTED_SOURCE_MANIFEST_SHA256,
        "superseded_field": "depmap.expected_matrix.gene_columns",
        "superseded_value": 17_289,
    }:
        raise ContractError("source schema amendment manifest binding drifted")
    header = amendment.get("header")
    if not isinstance(header, dict):
        raise ContractError("source schema amendment header audit is missing")
    expected_header = {
        "bytes": 254_613,
        "sha256": "1a72fa781f9514231b0631a8bc1ee7199659224c553f8084bbfb848e9e92fa1e",
        "index_column_label": "",
        "gene_columns": EXPECTED_DEPMAP_SHAPE[1],
        "unique_symbol_prefixes": EXPECTED_DEPMAP_SHAPE[1],
        "all_gene_labels_match_symbol_entrez_shape": True,
    }
    if header != expected_header:
        raise ContractError("source schema amendment header audit drifted")


def validate_source_candidate_freeze_boundary(
    boundary: Mapping[str, Any],
) -> None:
    if (
        sha256_file(SOURCE_CANDIDATE_FREEZE_BOUNDARY_PATH)
        != EXPECTED_SOURCE_CANDIDATE_FREEZE_BOUNDARY_SHA256
    ):
        raise ContractError("source-candidate freeze boundary changed")
    expected = {
        "schema_version": 1,
        "experiment": "114-depmap-coessentiality",
        "date_frozen": "2026-08-30",
        "status": (
            "implemented_not_authorized_pending_role_artifacts_and_independent_review"
        ),
        "source_candidate_freeze_authorized": False,
    }
    for key, value in expected.items():
        if boundary.get(key) != value:
            raise ContractError(f"source-candidate freeze boundary drifted: {key}")
    lineage = boundary.get("lineage")
    if not isinstance(lineage, dict) or lineage != {
        "feature_build_protocol_sha256": FEATURE_BUILD_PROTOCOL_SHA256,
        "acquisition_protocol_sha256": ACQUISITION_PROTOCOL_SHA256,
        "source_manifest_sha256": EXPECTED_SOURCE_MANIFEST_SHA256,
        "source_schema_amendment_sha256": EXPECTED_SOURCE_SCHEMA_AMENDMENT_SHA256,
        "source_receipt_sha256": EXPECTED_SOURCE_RECEIPT_SHA256,
        "feature_receipt_sha256": EXPECTED_FEATURE_RECEIPT_SHA256,
        "feature_build_implementation_sha256": FEATURE_BUILD_IMPLEMENTATION_SHA256,
        "role_source_manifest_path": (
            "experiments/114-depmap-coessentiality/ROLE_SOURCE_MANIFEST.json"
        ),
        "role_source_manifest_sha256": EXPECTED_ROLE_SOURCE_MANIFEST_SHA256,
    }:
        raise ContractError("source-candidate freeze lineage drifted")
    canonical_axes = boundary.get("canonical_axes")
    expected_axes = {
        "panel": {
            "path": "data/challenge/pert_counts.csv",
            "sha256": "f57edd7b912ebd718efc7ee9d0f334772513e7cc418d133ce525470e373b3276",
            "column": "target_gene",
            "count": 300,
            "ordered_sha256": (
                "cd3d00f1cd5da660609ad41d85b53e2d5557f7f5edb3e1f64559f78f7fb3cbc4"
            ),
        },
        "genes": {
            "path": "data/challenge/gene_names.csv",
            "sha256": "25bfa66715e186bebabce7ac788bbcea47e2bf59ca70be1f8f3a06f2f0e47201",
            "column": "gene_name",
            "count": 18_533,
            "ordered_sha256": (
                "c2b8dec58c02810d64018164155b655c4e7c7be451a3ae579965db193af33723"
            ),
        },
    }
    if canonical_axes != expected_axes:
        raise ContractError("source-candidate canonical-axis contract drifted")
    features = boundary.get("features")
    expected_features = {
        "depmap": {
            "artifact_path": (
                "experiments/114-depmap-coessentiality/work/source-features/"
                "depmap-features.npz"
            ),
            "artifact_sha256": (
                "7f42b585f788db95cdc345d44b44926bd76ce0fc5ae94c544fd2f6929b9ddc51"
            ),
            "receipt_path": (
                "experiments/114-depmap-coessentiality/work/source-features/"
                "feature-receipt.json"
            ),
            "receipt_sha256": EXPECTED_FEATURE_RECEIPT_SHA256,
            "genes": 17_916,
            "dimensions": 1_178,
            "gene_order_sha256": (
                "0f595cfb921ec95012e07ecd1928f00d29512cdadbbbc1b52733c3923112bb9c"
            ),
            "matrix_sha256": (
                "2444fadb3feec310b38823dd962694aea487ea9556e776c7fa545de7ebfe648e"
            ),
        },
        "string": {
            "artifact_path": "data/external/slim/gene_string_embeddings.v0.3.h5",
            "artifact_sha256": (
                "789416877b8701ef6f800106d26bf7bb97ea8e72744e6ab93e24933a717f247d"
            ),
            "genes": 21_688,
            "dimensions": 64,
        },
    }
    if features != expected_features:
        raise ContractError("source-candidate feature contract drifted")
    role_contract = boundary.get("role_manifest_contract")
    expected_role_contract = {
        "schema_version": 1,
        "required_stage": "source_only_role_inputs",
        "required_artifacts": ["incumbent", "source_bank"],
        "artifact_directory": (
            "experiments/114-depmap-coessentiality/work/roles/<role>"
        ),
        "incumbent_stage": "role_incumbent_source_only",
        "source_bank_stage": "role_source_bank_source_only",
        "all_artifacts_must_be_sealed_before_held_truth": True,
        "all_paths_and_receipts_must_be_hash_pinned_by_the_approved_role_manifest": True,
        "held_surface_source_id_must_be_forbidden": True,
        "source_ids_must_equal_the_role_allowlist": True,
        "canonical_target_axis_required": True,
        "canonical_gene_axis_required": True,
        "exact_preregistered_eligibility_required": True,
        "role_source_manifest_binding_required": True,
        "role_specific_experiment_107_operator_required": True,
        "held_surface_response_coordinates_must_be_excluded": True,
        "source_bank_target_universe": (
            "exact_top25_neighbour_closure_for_all_41_eligible_queries_and_three_"
            "noncontrol_arms"
        ),
    }
    if role_contract != expected_role_contract:
        raise ContractError("source-candidate role-manifest contract drifted")
    approved_roles = boundary.get("approved_roles")
    if not isinstance(approved_roles, dict):
        raise ContractError("source-candidate approved_roles must be an object")
    if approved_roles:
        raise ContractError(
            "checked-in boundary unexpectedly contains approved roles before review"
        )


def validate_role_source_manifest(value: Mapping[str, Any]) -> None:
    if sha256_file(ROLE_SOURCE_MANIFEST_PATH) != EXPECTED_ROLE_SOURCE_MANIFEST_SHA256:
        raise ContractError("role source manifest changed")
    expected = {
        "schema_version": 1,
        "experiment": "114-depmap-coessentiality",
        "status": "builders_defined_real_role_bundles_not_built_or_approved",
    }
    for key, item in expected.items():
        if value.get(key) != item:
            raise ContractError(f"role source manifest drifted: {key}")
    if tuple(value.get("roles", {})) != ("h1", "orion_hek293t"):
        raise ContractError("role source names or order drifted")
    authorization = value.get("authorization")
    if not isinstance(authorization, dict):
        raise ContractError("role source authorization is missing")
    if authorization.get("source_candidate_freeze") is not False:
        raise ContractError("role source manifest unexpectedly authorizes freeze")
    if authorization.get("held_truth_access") is not False:
        raise ContractError("role source manifest unexpectedly authorizes held truth")
    if authorization.get("approved_roles") != {}:
        raise ContractError("role source manifest contains premature role approvals")
    axes = value.get("canonical_axes", {})
    if axes.get("targets") != {
        "path": "data/challenge/pert_counts.csv",
        "sha256": "f57edd7b912ebd718efc7ee9d0f334772513e7cc418d133ce525470e373b3276",
        "column": "target_gene",
        "count": 300,
        "ordered_sha256": (
            "cd3d00f1cd5da660609ad41d85b53e2d5557f7f5edb3e1f64559f78f7fb3cbc4"
        ),
    }:
        raise ContractError("role source target-axis contract drifted")
    if axes.get("genes") != {
        "path": "data/challenge/gene_names.csv",
        "sha256": "25bfa66715e186bebabce7ac788bbcea47e2bf59ca70be1f8f3a06f2f0e47201",
        "column": "gene_name",
        "count": 18_533,
        "ordered_sha256": (
            "c2b8dec58c02810d64018164155b655c4e7c7be451a3ae579965db193af33723"
        ),
    }:
        raise ContractError("role source gene-axis contract drifted")


def validate_source_manifest(source: Mapping[str, Any]) -> None:
    if sha256_file(SOURCE_MANIFEST_PATH) != EXPECTED_SOURCE_MANIFEST_SHA256:
        raise ContractError("source manifest changed")
    depmap = source.get("depmap")
    if not isinstance(depmap, dict):
        raise ContractError("missing DepMap provenance")
    expected = {
        "release": "DepMap Public 24Q4",
        "article_id": 27993248,
        "doi": "10.25452/figshare.plus.27993248.v1",
    }
    for key, value in expected.items():
        if depmap.get(key) != value:
            raise ContractError(f"DepMap {key} drifted")
    file_record = depmap.get("file")
    if not isinstance(file_record, dict):
        raise ContractError("missing DepMap file record")
    if file_record.get("id") != 51064667:
        raise ContractError("wrong DepMap file id")
    if file_record.get("name") != "CRISPRGeneEffect.csv":
        raise ContractError("wrong DepMap filename")
    if file_record.get("size_bytes") != 428_678_699:
        raise ContractError("DepMap byte count drifted")
    if file_record.get("md5") != "6edf7ade09b9b34199210b559d4745d3":
        raise ContractError("DepMap MD5 drifted")
    if file_record.get("download_url") != (
        "https://ndownloader.figshare.com/files/51064667"
    ):
        raise ContractError("DepMap download URL drifted")


def download_depmap(
    *,
    protocol: Mapping[str, Any],
    source_manifest: Mapping[str, Any],
    confirmed: bool,
) -> dict[str, Any]:
    assert_authorized(protocol, "depmap_large_download")
    if not confirmed:
        raise ContractError("large download requires --confirm-large-download")
    csv_path = require_under_experiment(DEPMAP_CSV_PATH)
    receipt_path = require_under_experiment(DEPMAP_SOURCE_RECEIPT_PATH)
    if (
        csv_path.exists()
        or csv_path.is_symlink()
        or receipt_path.exists()
        or receipt_path.is_symlink()
    ):
        raise ContractError("refusing to overwrite existing DepMap source artifacts")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    record = source_manifest["depmap"]["file"]
    source_url = str(record["download_url"])
    if source_url != "https://ndownloader.figshare.com/files/51064667":
        raise ContractError("DepMap source URL is not the frozen allowlisted URL")
    temporary = csv_path.with_name(f".{csv_path.name}.{os.getpid()}.part")
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    total = 0
    try:
        request = urllib.request.Request(
            source_url,
            headers={"User-Agent": "vcc-exp114-source-acquisition/1"},
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            final_url = str(response.geturl())
            if not final_url.startswith("https://"):
                raise ContractError("DepMap download redirected outside HTTPS")
            length = response.headers.get("Content-Length")
            if length is not None and int(length) != int(record["size_bytes"]):
                raise ContractError(
                    "DepMap Content-Length differs from source manifest"
                )
            with temporary.open("xb") as handle:
                while block := response.read(8 << 20):
                    total += len(block)
                    if total > int(record["size_bytes"]):
                        raise ContractError(
                            "DepMap download exceeded expected byte count"
                        )
                    md5.update(block)
                    sha256.update(block)
                    handle.write(block)
                handle.flush()
                os.fsync(handle.fileno())
        if total != int(record["size_bytes"]):
            raise ContractError(
                "DepMap download byte count differs from source manifest"
            )
        if md5.hexdigest() != str(record["md5"]):
            raise ContractError("DepMap download MD5 differs from source manifest")
        os.link(temporary, csv_path)
        temporary.unlink()
        receipt = {
            "schema_version": 1,
            "experiment": "114-depmap-coessentiality",
            "stage": "depmap_source_acquisition",
            "source_url": source_url,
            "resolved_https_url": final_url,
            "figshare_file_id": int(record["id"]),
            "source_bytes": total,
            "source_md5": md5.hexdigest(),
            "source_sha256": sha256.hexdigest(),
            "source_csv": csv_path.relative_to(ROOT).as_posix(),
            "protocol_sha256": sha256_file(PROTOCOL_PATH),
            "source_manifest_sha256": sha256_file(SOURCE_MANIFEST_PATH),
            "implementation_sha256": sha256_file(Path(__file__)),
            "held_truth_opened": False,
        }
        atomic_json(receipt_path, receipt)
        return receipt
    except BaseException:
        temporary.unlink(missing_ok=True)
        if not receipt_path.exists():
            csv_path.unlink(missing_ok=True)
        raise


def verify_depmap_source(
    csv_path: Path,
    *,
    protocol: Mapping[str, Any],
    source_manifest: Mapping[str, Any],
) -> dict[str, Any]:
    assert_authorized(protocol, "depmap_large_download")
    if csv_path.resolve() != DEPMAP_CSV_PATH.resolve():
        raise ContractError("DepMap CSV is not the frozen experiment source path")
    if not csv_path.is_file() or csv_path.is_symlink():
        raise ContractError("DepMap CSV is missing or unsafe")
    if (
        not DEPMAP_SOURCE_RECEIPT_PATH.is_file()
        or DEPMAP_SOURCE_RECEIPT_PATH.is_symlink()
    ):
        raise ContractError("DepMap source receipt is missing or unsafe")
    receipt = load_json(DEPMAP_SOURCE_RECEIPT_PATH)
    amendment = load_json(SOURCE_SCHEMA_AMENDMENT_PATH)
    validate_source_schema_amendment(amendment)
    record = source_manifest["depmap"]["file"]
    expected = {
        "schema_version": 1,
        "experiment": "114-depmap-coessentiality",
        "stage": "depmap_source_acquisition",
        "source_url": record["download_url"],
        "figshare_file_id": record["id"],
        "source_bytes": record["size_bytes"],
        "source_md5": record["md5"],
        "source_csv": csv_path.relative_to(ROOT).as_posix(),
        "protocol_sha256": ACQUISITION_PROTOCOL_SHA256,
        "source_manifest_sha256": sha256_file(SOURCE_MANIFEST_PATH),
        "held_truth_opened": False,
    }
    for key, value in expected.items():
        if receipt.get(key) != value:
            raise ContractError(f"DepMap source receipt drifted: {key}")
    if csv_path.stat().st_size != int(record["size_bytes"]):
        raise ContractError("DepMap file byte count differs from source manifest")
    observed_md5 = md5_file(csv_path)
    observed_sha256 = sha256_file(csv_path)
    if observed_md5 != record["md5"] or receipt.get("source_md5") != observed_md5:
        raise ContractError("DepMap file MD5 differs from its source receipt")
    if receipt.get("source_sha256") != observed_sha256:
        raise ContractError("DepMap file SHA-256 differs from its source receipt")
    amendment_receipt = amendment.get("source_receipt")
    if not isinstance(amendment_receipt, dict) or amendment_receipt != {
        "path": DEPMAP_SOURCE_RECEIPT_PATH.relative_to(ROOT).as_posix(),
        "sha256": sha256_file(DEPMAP_SOURCE_RECEIPT_PATH),
    }:
        raise ContractError("source schema amendment receipt binding drifted")
    amendment_source = amendment.get("source_file")
    if not isinstance(amendment_source, dict) or amendment_source != {
        "path": csv_path.relative_to(ROOT).as_posix(),
        "bytes": csv_path.stat().st_size,
        "md5": observed_md5,
        "sha256": observed_sha256,
    }:
        raise ContractError("source schema amendment file binding drifted")
    return receipt


def require_local_input_hash(source: Mapping[str, Any], path: Path) -> str:
    try:
        relative = path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError as error:
        raise ContractError(f"local input escapes repository: {path}") from error
    expected = source.get("local_inputs", {}).get(relative)
    if not isinstance(expected, str):
        raise ContractError(f"source manifest lacks local input hash: {relative}")
    observed = sha256_file(path)
    if observed != expected:
        raise ContractError(f"local input hash drifted: {relative}")
    return observed


def derive_production_eligibility(
    protocol: Mapping[str, Any],
    source_manifest: Mapping[str, Any],
    *,
    panel_path: Path,
    incumbent_meta_path: Path,
    donor_pert_paths: Mapping[str, Path],
    perturbmatch_pert_path: Path,
) -> dict[str, Any]:
    require_local_input_hash(source_manifest, panel_path)
    panel = read_csv_column(panel_path, "target_gene")
    if len(panel) != 300:
        raise ContractError("challenge panel must contain exactly 300 targets")
    incumbent = load_json(incumbent_meta_path)
    expected_hash = protocol["incumbent_contract"]["reference_metadata_sha256"]
    if sha256_file(incumbent_meta_path) != expected_hash:
        raise ContractError("incumbent metadata hash drifted")
    frozen_fallback_cross_check = incumbent["effects"]["perturbmatch"][
        "neither_supported_preserved"
    ]
    require_local_input_hash(source_manifest, perturbmatch_pert_path)
    perturbmatch_targets = set(read_lines(perturbmatch_pert_path))
    for path in donor_pert_paths.values():
        require_local_input_hash(source_manifest, path)
    donor_sets = {
        name: set(read_lines(path)) for name, path in donor_pert_paths.items()
    }
    required_names = tuple(protocol["incumbent_contract"]["real_donors"])
    if tuple(donor_pert_paths) != required_names:
        raise ContractError("real-donor order drifted")
    counts = {
        target: sum(target in donor_sets[name] for name in required_names)
        for target in panel
    }
    observed_donor_counts = {
        name: sum(target in donor_sets[name] for target in panel)
        for name in required_names
    }
    if (
        observed_donor_counts
        != protocol["incumbent_contract"]["real_donor_target_counts"]
    ):
        raise ContractError("real-donor target counts drifted")
    fallback = [
        target
        for target in panel
        if target not in donor_sets["k562_gwps"] and target not in perturbmatch_targets
    ]
    if fallback != frozen_fallback_cross_check:
        raise ContractError(
            "direct K562/PerturbMatch support derivation disagrees with metadata cross-check"
        )
    low = [
        target for target in panel if counts[target] < 3 and target not in set(fallback)
    ]
    eligible = [
        target for target in panel if target in set(fallback) or counts[target] < 3
    ]
    expected = protocol["eligibility"]
    observed = {
        "incumbent_uncovered_targets": fallback,
        "low_confidence_targets": low,
        "eligible_targets": eligible,
    }
    for key, values in observed.items():
        if values != expected[key]:
            raise ContractError(f"derived {key} differs from preregistration")
    histogram = {
        str(value): sum(count == value for count in counts.values())
        for value in range(4)
    }
    histogram = {key: value for key, value in histogram.items() if value}
    if histogram != protocol["incumbent_contract"]["source_vote_count_histogram"]:
        raise ContractError("source vote histogram drifted")
    return {
        "panel_count": len(panel),
        "panel_ordered_sha256": sha256_lines(panel),
        "incumbent_uncovered_count": len(fallback),
        "low_confidence_count": len(low),
        "eligible_count": len(eligible),
        "eligible_ordered_sha256": sha256_lines(eligible),
        "source_vote_count_histogram": histogram,
        "real_donor_target_counts": observed_donor_counts,
        "reference_metadata_role": protocol["incumbent_contract"][
            "reference_metadata_role"
        ],
    }


def depmap_features_from_frame(frame: Any) -> GeneFeatures:
    """Apply pphouse's exact 24Q4 row transform to a DataFrame-like object."""
    columns = [str(value) for value in frame.columns]
    symbols = np.asarray([column.split(" (", 1)[0] for column in columns], dtype=str)
    if len(set(symbols.tolist())) != len(symbols):
        raise ContractError("DepMap symbol prefixes are not unique")
    matrix = frame.to_numpy(dtype=np.float32).T
    seen = np.isfinite(matrix).sum(axis=1)
    keep = seen >= 2
    genes = symbols[keep]
    values = matrix[keep]
    centered = np.where(np.isfinite(values), values, np.nan)
    with np.errstate(invalid="ignore"):
        means = np.nanmean(centered, axis=1, keepdims=True)
    centered = np.nan_to_num(centered - means, nan=0.0)
    norms = np.linalg.norm(centered, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return GeneFeatures(genes, centered / norms)


def build_depmap_features(
    csv_path: Path,
    output_path: Path,
    *,
    protocol: Mapping[str, Any],
    source_manifest: Mapping[str, Any],
) -> dict[str, Any]:
    assert_authorized(protocol, "depmap_large_download")
    assert_authorized(protocol, "source_feature_build")
    if csv_path.resolve() != DEPMAP_CSV_PATH.resolve():
        raise ContractError("DepMap CSV is not the frozen experiment source path")
    if output_path.resolve() != DEPMAP_FEATURES_PATH.resolve():
        raise ContractError("DepMap feature output is not the frozen experiment path")
    if (
        output_path.exists()
        or output_path.is_symlink()
        or DEPMAP_FEATURES_RECEIPT_PATH.exists()
        or DEPMAP_FEATURES_RECEIPT_PATH.is_symlink()
    ):
        raise ContractError("refusing to overwrite existing DepMap feature artifacts")
    source_receipt = verify_depmap_source(
        csv_path,
        protocol=protocol,
        source_manifest=source_manifest,
    )
    import pandas as pd

    frame = pd.read_csv(csv_path, index_col=0)
    representation = protocol.get("representations", {}).get("depmap_24q4", {})
    expected_shape = (
        representation.get("rows_cell_lines"),
        representation.get("columns_genes"),
    )
    if expected_shape != EXPECTED_DEPMAP_SHAPE:
        raise ContractError("amended DepMap protocol dimensions drifted")
    if frame.shape != EXPECTED_DEPMAP_SHAPE:
        raise ContractError(f"DepMap matrix shape drifted: {frame.shape}")
    features = depmap_features_from_frame(frame)
    receipt = {
        "schema_version": 1,
        "source_csv": str(csv_path.resolve()),
        "source_bytes": csv_path.stat().st_size,
        "source_md5": md5_file(csv_path),
        "source_sha256": sha256_file(csv_path),
        "source_receipt_sha256": sha256_file(DEPMAP_SOURCE_RECEIPT_PATH),
        "genes": int(len(features.genes)),
        "dimensions": int(features.matrix.shape[1]),
        "gene_order_sha256": sha256_lines(features.genes),
        "matrix_sha256": sha256_array(features.matrix),
        "protocol_sha256": sha256_file(PROTOCOL_PATH),
        "source_manifest_sha256": sha256_file(SOURCE_MANIFEST_PATH),
        "source_schema_amendment_sha256": sha256_file(SOURCE_SCHEMA_AMENDMENT_PATH),
        "acquisition_protocol_sha256": source_receipt["protocol_sha256"],
        "implementation_sha256": sha256_file(Path(__file__)),
        "transformation_version": "pphouse-79cd460-float32-bit-exact-v1",
        "held_truth_opened": False,
    }
    atomic_npz(output_path, genes=features.genes, matrix=features.matrix)
    receipt["output_npz"] = output_path.relative_to(ROOT).as_posix()
    receipt["output_npz_sha256"] = sha256_file(output_path)
    receipt["source_acquisition"] = {
        "source_url": source_receipt["source_url"],
        "figshare_file_id": source_receipt["figshare_file_id"],
    }
    atomic_json(DEPMAP_FEATURES_RECEIPT_PATH, receipt)
    return receipt


def load_npz_features(path: Path) -> GeneFeatures:
    with np.load(path, allow_pickle=False) as archive:
        if set(archive.files) != {"genes", "matrix"}:
            raise ContractError(f"{path}: feature NPZ keys must be genes,matrix")
        return GeneFeatures(archive["genes"], archive["matrix"])


def load_string_features(path: Path) -> GeneFeatures:
    import h5py

    genes: list[str] = []
    rows: list[np.ndarray] = []
    with h5py.File(path, "r") as handle:
        for gene in sorted(handle.keys()):
            value = np.asarray(handle[gene], dtype=np.float32)
            if value.shape != (64,):
                raise ContractError(f"STRING row {gene!r} has shape {value.shape}")
            genes.append(str(gene))
            rows.append(value)
    if len(genes) != 21_688:
        raise ContractError(f"STRING artifact has {len(genes)} rows, expected 21688")
    return GeneFeatures(np.asarray(genes), row_normalize(np.vstack(rows)))


def validate_canonical_axes(
    boundary: Mapping[str, Any],
    source_manifest: Mapping[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    axes = boundary["canonical_axes"]
    panel_record = axes["panel"]
    gene_record = axes["genes"]
    panel_path = verify_fixed_record(panel_record, "canonical panel")
    gene_path = verify_fixed_record(gene_record, "canonical gene axis")
    if panel_path != PANEL_PATH.resolve() or gene_path != GENE_NAMES_PATH.resolve():
        raise ContractError("canonical axis paths drifted from the challenge bundle")
    for path, record in ((panel_path, panel_record), (gene_path, gene_record)):
        relative = path.relative_to(ROOT.resolve()).as_posix()
        if source_manifest.get("local_inputs", {}).get(relative) != record["sha256"]:
            raise ContractError(
                f"canonical axis is not bound by source manifest: {relative}"
            )
    panel = np.asarray(read_csv_column(panel_path, str(panel_record["column"])))
    genes = np.asarray(read_csv_column(gene_path, str(gene_record["column"])))
    for label, values, record in (
        ("panel", panel, panel_record),
        ("gene", genes, gene_record),
    ):
        if len(values) != int(record["count"]):
            raise ContractError(f"canonical {label} count drifted")
        if sha256_lines(values) != record["ordered_sha256"]:
            raise ContractError(f"canonical {label} order drifted")
    return panel, genes


def validate_depmap_feature_artifact(
    boundary: Mapping[str, Any],
    protocol: Mapping[str, Any],
    source_manifest: Mapping[str, Any],
) -> GeneFeatures:
    record = boundary["features"]["depmap"]
    artifact_path = resolve_repository_path(record["artifact_path"], "DepMap feature")
    receipt_path = resolve_repository_path(
        record["receipt_path"], "DepMap feature receipt"
    )
    if artifact_path != DEPMAP_FEATURES_PATH.resolve():
        raise ContractError("DepMap feature artifact path drifted")
    if receipt_path != DEPMAP_FEATURES_RECEIPT_PATH.resolve():
        raise ContractError("DepMap feature receipt path drifted")
    for path, expected_hash, label in (
        (artifact_path, record["artifact_sha256"], "DepMap feature artifact"),
        (receipt_path, record["receipt_sha256"], "DepMap feature receipt"),
        (
            DEPMAP_SOURCE_RECEIPT_PATH,
            boundary["lineage"]["source_receipt_sha256"],
            "DepMap source receipt",
        ),
    ):
        reject_symlink_components(path, label)
        if not path.is_file() or sha256_file(path) != expected_hash:
            raise ContractError(f"{label} is missing or changed")
    source_receipt = verify_depmap_source(
        DEPMAP_CSV_PATH,
        protocol=protocol,
        source_manifest=source_manifest,
    )
    receipt = load_json(receipt_path)
    expected_receipt = {
        "schema_version": 1,
        "source_bytes": DEPMAP_CSV_PATH.stat().st_size,
        "source_md5": md5_file(DEPMAP_CSV_PATH),
        "source_sha256": sha256_file(DEPMAP_CSV_PATH),
        "source_receipt_sha256": sha256_file(DEPMAP_SOURCE_RECEIPT_PATH),
        "genes": int(record["genes"]),
        "dimensions": int(record["dimensions"]),
        "gene_order_sha256": record["gene_order_sha256"],
        "matrix_sha256": record["matrix_sha256"],
        "protocol_sha256": FEATURE_BUILD_PROTOCOL_SHA256,
        "source_manifest_sha256": EXPECTED_SOURCE_MANIFEST_SHA256,
        "source_schema_amendment_sha256": (EXPECTED_SOURCE_SCHEMA_AMENDMENT_SHA256),
        "acquisition_protocol_sha256": ACQUISITION_PROTOCOL_SHA256,
        "implementation_sha256": FEATURE_BUILD_IMPLEMENTATION_SHA256,
        "transformation_version": "pphouse-79cd460-float32-bit-exact-v1",
        "held_truth_opened": False,
        "output_npz": DEPMAP_FEATURES_PATH.relative_to(ROOT).as_posix(),
        "output_npz_sha256": record["artifact_sha256"],
        "source_acquisition": {
            "source_url": source_receipt["source_url"],
            "figshare_file_id": source_receipt["figshare_file_id"],
        },
    }
    for key, value in expected_receipt.items():
        if receipt.get(key) != value:
            raise ContractError(f"DepMap feature receipt drifted: {key}")
    source_csv = receipt.get("source_csv")
    if (
        not isinstance(source_csv, str)
        or Path(source_csv).resolve() != DEPMAP_CSV_PATH.resolve()
    ):
        raise ContractError("DepMap feature receipt source path drifted")
    features = load_npz_features(artifact_path)
    if features.matrix.shape != (int(record["genes"]), int(record["dimensions"])):
        raise ContractError("DepMap feature dimensions drifted")
    if sha256_lines(features.genes) != record["gene_order_sha256"]:
        raise ContractError("DepMap feature gene order drifted")
    if sha256_array(features.matrix) != record["matrix_sha256"]:
        raise ContractError("DepMap logical feature matrix drifted")
    norms = np.linalg.norm(features.matrix, axis=1)
    nonzero = norms > 0
    if np.any(np.abs(norms[nonzero] - 1.0) > 2e-6):
        raise ContractError("DepMap feature rows are not normalized")
    return features


def validate_string_feature_artifact(
    boundary: Mapping[str, Any],
    source_manifest: Mapping[str, Any],
) -> GeneFeatures:
    record = boundary["features"]["string"]
    path = resolve_repository_path(record["artifact_path"], "STRING feature")
    if path != STRING_FEATURES_PATH.resolve():
        raise ContractError("STRING feature path drifted")
    reject_symlink_components(path, "STRING feature")
    if not path.is_file() or sha256_file(path) != record["artifact_sha256"]:
        raise ContractError("STRING feature artifact is missing or changed")
    relative = path.relative_to(ROOT.resolve()).as_posix()
    if (
        source_manifest.get("local_inputs", {}).get(relative)
        != record["artifact_sha256"]
    ):
        raise ContractError("STRING feature is not bound by source manifest")
    features = load_string_features(path)
    if features.matrix.shape != (int(record["genes"]), int(record["dimensions"])):
        raise ContractError("STRING feature dimensions drifted")
    return features


def blend_features(
    left: GeneFeatures,
    right: GeneFeatures,
    weight: float = 1.0,
) -> GeneFeatures:
    if weight != 1.0:
        raise ContractError("experiment 114 permits only published blend weight 1.0")
    genes = np.asarray(sorted(set(left.genes) | set(right.genes)), dtype=str)
    left_rows = np.zeros((len(genes), left.matrix.shape[1]), dtype=np.float32)
    right_rows = np.zeros((len(genes), right.matrix.shape[1]), dtype=np.float32)
    left_index = {gene: index for index, gene in enumerate(left.genes)}
    right_index = {gene: index for index, gene in enumerate(right.genes)}
    for row, gene in enumerate(genes):
        if gene in left_index:
            left_rows[row] = left.matrix[left_index[gene]]
        if gene in right_index:
            right_rows[row] = right.matrix[right_index[gene]]
    return GeneFeatures(genes, row_normalize(np.hstack([left_rows, right_rows])))


def load_source_bank(path: Path) -> SourceBank:
    with np.load(path, allow_pickle=False) as archive:
        required = {"targets", "genes", "effects", "covered", "source_ids"}
        if set(archive.files) != required:
            raise ContractError(f"{path}: source NPZ keys must be {sorted(required)}")
        targets = archive["targets"].astype(str)
        genes = archive["genes"].astype(str)
        effects = np.asarray(archive["effects"], dtype=np.float64)
        covered = np.asarray(archive["covered"], dtype=bool)
        source_ids = archive["source_ids"].astype(str)
    if targets.ndim != 1 or len(targets) == 0 or source_ids.shape != targets.shape:
        raise ContractError("source target/source axes disagree")
    if np.any(targets == "") or np.any(source_ids == ""):
        raise ContractError("source target/source labels must be nonempty")
    if genes.ndim != 1 or len(set(genes.tolist())) != len(genes):
        raise ContractError("source gene axis is invalid")
    if effects.shape != (len(targets), len(genes)) or covered.shape != effects.shape:
        raise ContractError("source effect/coverage matrices do not align")
    if not np.isfinite(effects).all():
        raise ContractError("source effects contain nonfinite values")
    if np.any(effects[~covered] != 0):
        raise ContractError("uncovered source coordinates must be literal zero")
    pairs = list(zip(source_ids.tolist(), targets.tolist(), strict=True))
    if len(pairs) != len(set(pairs)):
        raise ContractError("source bank repeats a source_id,target pair")
    return SourceBank(targets, genes, effects, covered, source_ids)


def build_consensus(bank: SourceBank) -> ConsensusBank:
    targets = np.asarray(sorted(set(bank.targets.tolist())), dtype=str)
    effects = np.zeros((len(targets), len(bank.genes)), dtype=np.float64)
    covered = np.zeros_like(effects, dtype=bool)
    carrier_count = np.zeros(len(targets), dtype=np.int16)
    for out_row, target in enumerate(targets):
        rows = np.flatnonzero(bank.targets == target)
        carrier_count[out_row] = len(set(bank.source_ids[rows].tolist()))
        mask = bank.covered[rows]
        numerator = (bank.effects[rows] * mask).sum(axis=0)
        denominator = mask.sum(axis=0)
        covered[out_row] = denominator > 0
        effects[out_row] = np.divide(
            numerator,
            denominator,
            out=np.zeros(len(bank.genes), dtype=np.float64),
            where=denominator > 0,
        )
    return ConsensusBank(targets, bank.genes.copy(), effects, covered, carrier_count)


def feature_rows(
    features: GeneFeatures, genes: Sequence[str]
) -> tuple[np.ndarray, np.ndarray]:
    index = {gene: row for row, gene in enumerate(features.genes.tolist())}
    matrix = np.zeros((len(genes), features.matrix.shape[1]), dtype=np.float32)
    available = np.zeros(len(genes), dtype=bool)
    for row, gene in enumerate(map(str, genes)):
        source = index.get(gene)
        if source is None:
            continue
        vector = features.matrix[source]
        if np.linalg.norm(vector) <= 0:
            continue
        matrix[row] = vector
        available[row] = True
    return matrix, available


def neighbour_weights(similarity: np.ndarray) -> np.ndarray:
    sim = np.clip(
        np.nan_to_num(similarity, nan=0.0, posinf=0.0, neginf=0.0),
        0.0,
        None,
    )
    if sim.size == 0 or not np.any(sim > 0):
        return np.zeros_like(sim)
    positive = np.flatnonzero(sim > 0)
    order = positive[np.argsort(-sim[positive], kind="stable")]
    selected = order[:NEIGHBOURS]
    weights = np.zeros_like(sim)
    weights[selected] = sim[selected] ** NEIGHBOUR_POWER
    return weights


def coverage_aware_average(
    effects: np.ndarray,
    covered: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    selected = np.flatnonzero(weights > 0)
    if len(selected) == 0:
        return np.zeros(effects.shape[1], dtype=np.float64)
    w = weights[selected, None]
    mask = covered[selected]
    numerator = (effects[selected] * w * mask).sum(axis=0)
    denominator = (w * mask).sum(axis=0)
    return np.divide(
        numerator,
        denominator,
        out=np.zeros(effects.shape[1], dtype=np.float64),
        where=denominator > 0,
    )


def build_arm(
    *,
    incumbent_targets: Sequence[str],
    incumbent_genes: Sequence[str],
    incumbent_effects: np.ndarray,
    consensus: ConsensusBank,
    features: GeneFeatures,
    eligible_targets: set[str],
    panel_targets: set[str],
    query_feature_target: Mapping[str, str] | None = None,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    targets = np.asarray(incumbent_targets).astype(str)
    genes = np.asarray(incumbent_genes).astype(str)
    incumbent = np.asarray(incumbent_effects)
    if incumbent.shape != (len(targets), len(genes)):
        raise ContractError("incumbent effect matrix does not align")
    if not np.array_equal(genes, consensus.genes):
        raise ContractError("incumbent and source gene axes differ")
    panel_positions = np.asarray(
        [index for index, gene in enumerate(genes) if gene in panel_targets],
        dtype=np.int64,
    )
    allowed = np.ones(len(genes), dtype=bool)
    allowed[panel_positions] = False
    if np.any(incumbent[:, panel_positions] != 0):
        raise ContractError("incumbent has nonzero VCC target coordinates")

    source_feature, source_available = feature_rows(features, consensus.targets)
    feature_index = {gene: row for row, gene in enumerate(features.genes.tolist())}
    result = incumbent.copy()
    confidence = np.abs(np.asarray(incumbent, dtype=np.float64))
    changed: list[str] = []
    unsupported: list[str] = []
    neighbor_audit: dict[str, Any] = {}

    for row, target in enumerate(targets):
        if target not in eligible_targets:
            continue
        feature_target = (
            target
            if query_feature_target is None
            else query_feature_target.get(target, target)
        )
        query_row = feature_index.get(feature_target)
        if query_row is None or np.linalg.norm(features.matrix[query_row]) <= 0:
            unsupported.append(target)
            continue
        similarity = source_feature @ features.matrix[query_row]
        similarity[~source_available] = -np.inf
        similarity[consensus.targets == target] = -np.inf
        weights = neighbour_weights(similarity)
        if weights.sum() <= 0:
            unsupported.append(target)
            continue
        candidate = coverage_aware_average(
            consensus.effects, consensus.covered, weights
        )
        candidate[~allowed] = 0.0
        incumbent_norm = float(
            np.linalg.norm(np.asarray(incumbent[row], dtype=np.float64)[allowed])
        )
        candidate_norm = float(np.linalg.norm(candidate[allowed]))
        if (
            not np.isfinite(candidate).all()
            or incumbent_norm <= 0
            or candidate_norm <= 0
        ):
            unsupported.append(target)
            continue
        candidate *= incumbent_norm / candidate_norm
        candidate[~allowed] = 0.0
        result[row] = candidate.astype(result.dtype, copy=False)
        confidence[row] = np.abs(np.asarray(result[row], dtype=np.float64))
        changed.append(target)
        selected = np.flatnonzero(weights > 0)
        ordered = selected[np.argsort(-weights[selected], kind="stable")]
        neighbor_audit[target] = {
            "feature_target": feature_target,
            "targets": consensus.targets[ordered].tolist(),
            "weights": weights[ordered].tolist(),
            "same_target_source_rows": int(
                np.sum(consensus.targets[ordered] == target)
            ),
        }

    ineligible = np.asarray([target not in eligible_targets for target in targets])
    if not np.array_equal(result[ineligible], incumbent[ineligible]):
        raise AssertionError("ineligible rows changed")
    unsupported_rows = np.asarray([target in set(unsupported) for target in targets])
    if not np.array_equal(result[unsupported_rows], incumbent[unsupported_rows]):
        raise AssertionError("unsupported rows changed")
    if np.any(result[:, panel_positions] != 0):
        raise AssertionError("candidate has nonzero VCC target coordinates")
    if any(item["same_target_source_rows"] != 0 for item in neighbor_audit.values()):
        raise AssertionError("same-target source signature survived leave-target-out")
    return (
        result,
        confidence,
        {
            "changed_targets": changed,
            "unsupported_targets": unsupported,
            "ineligible_rows_bit_identical": True,
            "unsupported_rows_bit_identical": True,
            "panel_coordinates_zero": True,
            "neighbor_audit": neighbor_audit,
        },
    )


def permuted_query_labels(
    fallback_targets: Sequence[str],
    low_confidence_targets: Sequence[str],
    draw: int,
) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for name, values in (
        ("incumbent_uncovered", tuple(map(str, fallback_targets))),
        ("low_confidence", tuple(map(str, low_confidence_targets))),
    ):
        rng = np.random.default_rng(
            stable_seed("experiment-114-target-label-null-v1", name, draw)
        )
        permuted = np.asarray(values)[rng.permutation(len(values))]
        mapping.update(zip(values, permuted.tolist(), strict=True))
    return mapping


def raw_uncentered_pds(
    prediction: np.ndarray,
    truth: np.ndarray,
    allowed_mask: np.ndarray | None = None,
) -> tuple[float, np.ndarray]:
    prediction = np.asarray(prediction, dtype=np.float64)
    truth = np.asarray(truth, dtype=np.float64)
    if allowed_mask is not None:
        allowed = np.asarray(allowed_mask, dtype=bool)
        if allowed.shape != (prediction.shape[1],):
            raise ContractError("PDS allowed-coordinate mask does not align")
        prediction = prediction[:, allowed]
        truth = truth[:, allowed]
    pred = row_normalize(prediction)
    real = row_normalize(truth)
    if pred.shape != real.shape or pred.shape[0] < 2:
        raise ContractError("PDS prediction/truth shapes are invalid")
    scores = pred @ real.T
    own = np.diag(scores)[:, None]
    better = (scores > own).sum(axis=1)
    ties = (scores == own).sum(axis=1) - 1
    per_target = 1.0 - (better + 0.5 * ties) / (len(scores) - 1)
    return float(per_target.mean()), per_target


def source_blind_directional_precision(
    prediction: np.ndarray,
    confidence: np.ndarray,
    truth_sign: np.ndarray,
    *,
    fraction: float = DIRECTIONAL_FRACTION,
    allowed_mask: np.ndarray | None = None,
) -> tuple[float, np.ndarray]:
    pred = np.asarray(prediction, dtype=np.float64)
    conf = np.asarray(confidence, dtype=np.float64)
    truth = np.asarray(truth_sign, dtype=np.int8)
    if pred.shape != conf.shape or pred.shape != truth.shape:
        raise ContractError("directional arrays do not align")
    if not (0 < fraction <= 1):
        raise ContractError("directional fraction must lie in (0,1]")
    allowed = (
        np.ones(pred.shape[1], dtype=bool)
        if allowed_mask is None
        else np.asarray(allowed_mask, dtype=bool)
    )
    if allowed.shape != (pred.shape[1],) or not np.any(allowed):
        raise ContractError("directional allowed-coordinate mask is invalid")
    allowed_rows = np.flatnonzero(allowed)
    count = max(1, int(math.ceil(len(allowed_rows) * fraction)))
    per_target = np.zeros(pred.shape[0], dtype=np.float64)
    for row in range(pred.shape[0]):
        ranked = np.argsort(-conf[row, allowed_rows], kind="stable")[:count]
        order = allowed_rows[ranked]
        predicted_sign = np.sign(pred[row, order]).astype(np.int8)
        per_target[row] = np.mean(
            (truth[row, order] != 0) & (predicted_sign == truth[row, order])
        )
    return float(per_target.mean()), per_target


def evaluate_surface(
    arms: Mapping[str, np.ndarray],
    confidence: Mapping[str, np.ndarray],
    truth_effect: np.ndarray,
    truth_sign: np.ndarray,
    allowed_mask: np.ndarray,
) -> dict[str, Any]:
    if set(arms) != set(ARMS) or set(confidence) != set(ARMS):
        raise ContractError("evaluation arms are incomplete")
    if not np.array_equal(arms["identity"], arms["current_fallback"]):
        raise ContractError("identity and current_fallback are not array-identical")
    output: dict[str, Any] = {}
    for name in ARMS:
        pds, per_target_pds = raw_uncentered_pds(arms[name], truth_effect, allowed_mask)
        precision, per_target_precision = source_blind_directional_precision(
            arms[name],
            confidence[name],
            truth_sign,
            allowed_mask=allowed_mask,
        )
        output[name] = {
            "raw_pds": pds,
            "directional_precision": precision,
            "per_target_raw_pds": per_target_pds.tolist(),
            "per_target_directional_precision": per_target_precision.tolist(),
        }
    return output


def promotion_decision(
    surface_results: Mapping[str, Mapping[str, Any]],
    permutation_macro_deltas: Sequence[float],
) -> dict[str, Any]:
    if not surface_results:
        raise ContractError("promotion requires at least one clean held surface")
    deltas = []
    precision_deltas = []
    for surface, result in surface_results.items():
        primary = result[PRIMARY_ARM]
        fallback = result["current_fallback"]
        deltas.append(float(primary["raw_pds"]) - float(fallback["raw_pds"]))
        precision_deltas.append(
            float(primary["directional_precision"])
            - float(fallback["directional_precision"])
        )
        if deltas[-1] < 0:
            return {"promote": False, "reason": f"negative PDS delta on {surface}"}
        if precision_deltas[-1] < 0:
            return {
                "promote": False,
                "reason": f"negative directional precision delta on {surface}",
            }
    macro_pds = float(np.mean(deltas))
    macro_precision = float(np.mean(precision_deltas))
    if macro_pds < 0.003:
        return {"promote": False, "reason": "macro raw PDS delta below 0.003"}
    if macro_precision <= 0:
        return {"promote": False, "reason": "directional precision did not improve"}
    if len(permutation_macro_deltas) != NULL_DRAWS:
        raise ContractError("promotion requires exactly 100 permutation nulls")
    null_p95 = float(np.quantile(permutation_macro_deltas, 0.95, method="higher"))
    if macro_pds <= null_p95:
        return {"promote": False, "reason": "observed PDS gain did not exceed null p95"}
    return {
        "promote": True,
        "macro_raw_pds_delta": macro_pds,
        "macro_directional_precision_delta": macro_precision,
        "permutation_p95": null_p95,
    }


def assert_authorized(protocol: Mapping[str, Any], key: str) -> None:
    if protocol["current_authorization"].get(key) is not True:
        raise ContractError(f"{key}=false in the frozen protocol")


def verify_path_record(record: Mapping[str, Any], label: str) -> Path:
    path_value = record.get("path")
    expected_hash = record.get("sha256")
    if not isinstance(path_value, str) or not isinstance(expected_hash, str):
        raise ContractError(f"{label} requires path and sha256")
    path = Path(path_value).expanduser().resolve()
    if not path.is_file():
        raise ContractError(f"{label} is missing: {path}")
    if sha256_file(path) != expected_hash:
        raise ContractError(f"{label} hash drifted")
    return path


def load_incumbent(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as archive:
        required = {
            "targets",
            "genes",
            "effects",
            "incumbent_fallback",
            "source_carrier_count",
            "expected_source_count",
        }
        if set(archive.files) != required:
            raise ContractError(
                f"{path}: incumbent NPZ keys must be {sorted(required)}"
            )
        output = {name: archive[name] for name in required}
    targets = output["targets"].astype(str)
    genes = output["genes"].astype(str)
    effects = np.asarray(output["effects"])
    fallback = np.asarray(output["incumbent_fallback"], dtype=bool)
    carriers = np.asarray(output["source_carrier_count"], dtype=np.int16)
    expected = np.asarray(output["expected_source_count"], dtype=np.int16)
    if targets.ndim != 1 or genes.ndim != 1:
        raise ContractError("incumbent axes must be vectors")
    if len(set(targets.tolist())) != len(targets):
        raise ContractError("incumbent targets contain duplicates")
    if len(set(genes.tolist())) != len(genes):
        raise ContractError("incumbent genes contain duplicates")
    if effects.shape != (len(targets), len(genes)):
        raise ContractError("incumbent effects do not align")
    if fallback.shape != targets.shape or carriers.shape != targets.shape:
        raise ContractError("incumbent eligibility metadata does not align")
    if expected.shape not in {(), targets.shape}:
        raise ContractError(
            "expected_source_count must be scalar or one value per target"
        )
    expected_full = (
        np.full(len(targets), int(expected), dtype=np.int16)
        if expected.shape == ()
        else expected
    )
    if (
        np.any(carriers < 0)
        or np.any(carriers > expected_full)
        or np.any(expected_full <= 0)
    ):
        raise ContractError("incumbent carrier counts are invalid")
    return {
        "targets": targets,
        "genes": genes,
        "effects": effects,
        "incumbent_fallback": fallback,
        "source_carrier_count": carriers,
        "expected_source_count": expected_full,
    }


def derive_surface_eligibility(incumbent: Mapping[str, np.ndarray]) -> dict[str, Any]:
    targets = incumbent["targets"].astype(str)
    fallback = np.asarray(incumbent["incumbent_fallback"], dtype=bool)
    carriers = np.asarray(incumbent["source_carrier_count"], dtype=np.int16)
    expected = np.asarray(incumbent["expected_source_count"], dtype=np.int16)
    low = (carriers < expected) & ~fallback
    eligible = fallback | low
    return {
        "incumbent_uncovered_targets": targets[fallback].tolist(),
        "low_confidence_targets": targets[low].tolist(),
        "eligible_targets": targets[eligible].tolist(),
        "incumbent_uncovered_ordered_sha256": sha256_lines(targets[fallback]),
        "low_confidence_ordered_sha256": sha256_lines(targets[low]),
        "eligible_ordered_sha256": sha256_lines(targets[eligible]),
    }


def load_freeze_boundary() -> dict[str, Any]:
    boundary = load_json(SOURCE_CANDIDATE_FREEZE_BOUNDARY_PATH)
    validate_source_candidate_freeze_boundary(boundary)
    return boundary


def approved_role_for_manifest(
    manifest_path: Path,
    boundary: Mapping[str, Any],
) -> tuple[str, Mapping[str, Any]]:
    resolved = require_under_experiment(manifest_path)
    reject_symlink_components(resolved, "role manifest")
    approved = boundary.get("approved_roles", {})
    for role, record in approved.items():
        if not isinstance(record, dict):
            raise ContractError(f"approved role record is invalid: {role}")
        required = {
            "manifest_path",
            "manifest_sha256",
            "surface",
            "role_directory",
            "output_directory",
            "incumbent_endpoint",
            "source_effect_space",
            "allowed_source_ids",
            "forbidden_source_ids",
            "held_surface_source_id",
        }
        if set(record) != required:
            raise ContractError(f"approved role fields drifted: {role}")
        allowed = record["allowed_source_ids"]
        forbidden = record["forbidden_source_ids"]
        held_source = record["held_surface_source_id"]
        if (
            not isinstance(allowed, list)
            or not allowed
            or len(allowed) != len(set(allowed))
            or not isinstance(forbidden, list)
            or len(forbidden) != len(set(forbidden))
            or held_source not in forbidden
            or set(allowed) & set(forbidden)
        ):
            raise ContractError(f"approved role source policy is invalid: {role}")
        expected_path = resolve_repository_path(
            str(record.get("manifest_path", "")), "role manifest"
        )
        if expected_path == resolved:
            if sha256_file(resolved) != record.get("manifest_sha256"):
                raise ContractError(f"approved role manifest hash drifted: {role}")
            return str(role), record
    raise ContractError("role manifest is not hash-pinned in the freeze boundary")


def validate_role_manifest(
    value: Mapping[str, Any],
    *,
    role: str,
    approved: Mapping[str, Any],
) -> Path:
    expected = {
        "schema_version": 1,
        "experiment": "114-depmap-coessentiality",
        "stage": "source_only_role_inputs",
        "role": role,
        "surface": approved.get("surface"),
        "contains_held_truth": False,
        "source_only": True,
        "sealed_before_held_truth": True,
    }
    for key, expected_value in expected.items():
        if value.get(key) != expected_value:
            raise ContractError(f"role manifest drifted: {key}")
    directory_value = approved.get("role_directory")
    if not isinstance(directory_value, str):
        raise ContractError("approved role lacks role_directory")
    role_directory = resolve_repository_path(directory_value, "role directory")
    if role_directory != (HERE / "work" / "roles" / role).resolve():
        raise ContractError("approved role directory is not canonical")
    if value.get("role_directory") != directory_value:
        raise ContractError("role manifest directory drifted")
    expected_role_source = {
        "path": ROLE_SOURCE_MANIFEST_PATH.relative_to(ROOT).as_posix(),
        "sha256": EXPECTED_ROLE_SOURCE_MANIFEST_SHA256,
    }
    if value.get("role_source_manifest") != expected_role_source:
        raise ContractError("role manifest source-plan binding drifted")
    verify_fixed_record(expected_role_source, f"{role} role source manifest")
    for name in ("incumbent", "source_bank"):
        group = value.get(name)
        if not isinstance(group, dict) or set(group) != {"artifact", "receipt"}:
            raise ContractError(f"role manifest lacks strict {name} records")
        verify_fixed_record(
            group["artifact"],
            f"{role} {name} artifact",
            required_under=role_directory,
        )
        verify_fixed_record(
            group["receipt"],
            f"{role} {name} receipt",
            required_under=role_directory,
        )
    return role_directory


def validate_role_receipt(
    receipt: Mapping[str, Any],
    *,
    stage: str,
    role: str,
    surface: str,
    artifact_path: Path,
    artifact_sha256: str,
) -> None:
    expected = {
        "schema_version": 1,
        "experiment": "114-depmap-coessentiality",
        "stage": stage,
        "role": role,
        "surface": surface,
        "source_only": True,
        "contains_held_truth": False,
        "sealed_before_held_truth": True,
        "artifact_path": artifact_path.relative_to(ROOT.resolve()).as_posix(),
        "artifact_sha256": artifact_sha256,
    }
    for key, value in expected.items():
        if receipt.get(key) != value:
            raise ContractError(f"{role} {stage} receipt drifted: {key}")
    expected_role_source = {
        "path": ROLE_SOURCE_MANIFEST_PATH.relative_to(ROOT).as_posix(),
        "sha256": EXPECTED_ROLE_SOURCE_MANIFEST_SHA256,
    }
    if receipt.get("role_source_manifest") != expected_role_source:
        raise ContractError(f"{role} {stage} role source binding drifted")
    verify_fixed_record(expected_role_source, f"{role} {stage} role source manifest")
    builder = receipt.get("builder")
    if (
        not isinstance(builder, dict)
        or not isinstance(builder.get("path"), str)
        or not isinstance(builder.get("sha256"), str)
    ):
        raise ContractError(f"{role} {stage} receipt lacks a pinned builder")
    verify_fixed_record(builder, f"{role} {stage} builder")
    source_inputs = receipt.get("source_inputs")
    if not isinstance(source_inputs, list) or not source_inputs:
        raise ContractError(f"{role} {stage} receipt lacks source input provenance")
    for index, record in enumerate(source_inputs):
        if (
            not isinstance(record, dict)
            or not isinstance(record.get("path"), str)
            or not isinstance(record.get("sha256"), str)
            or record.get("source_only") is not True
            or record.get("contains_held_truth") is not False
        ):
            raise ContractError(
                f"{role} {stage} source input {index} lacks strict provenance"
            )
        verify_fixed_record(record, f"{role} {stage} source input {index}")


def validate_array_hashes(
    receipt: Mapping[str, Any],
    arrays: Mapping[str, np.ndarray],
    names: Sequence[str],
    label: str,
) -> None:
    hashes = receipt.get("array_sha256")
    if not isinstance(hashes, dict) or set(hashes) != set(names):
        raise ContractError(f"{label} receipt array hash set drifted")
    for name in names:
        if hashes[name] != sha256_array(np.asarray(arrays[name])):
            raise ContractError(f"{label} array hash drifted: {name}")


def validate_exact_eligibility(
    incumbent: Mapping[str, np.ndarray],
    protocol: Mapping[str, Any],
    source_manifest: Mapping[str, Any],
) -> dict[str, Any]:
    targets = incumbent["targets"].astype(str)
    eligibility = protocol["eligibility"]
    fallback_targets = eligibility["incumbent_uncovered_targets"]
    low_targets = eligibility["low_confidence_targets"]
    eligible_targets = eligibility["eligible_targets"]
    fallback = np.isin(targets, np.asarray(fallback_targets))
    low = np.isin(targets, np.asarray(low_targets))
    expected = np.asarray(incumbent["expected_source_count"], dtype=np.int16)
    carriers = np.asarray(incumbent["source_carrier_count"], dtype=np.int16)
    donor_paths = {
        "k562_gwps": ROOT / "data/donors/k562_gwps/perts.txt",
        "orion_hct116": ROOT / "data/donors/orion_hct116/perts.txt",
        "orion_hek293t": ROOT / "data/donors/orion_hek293t/perts.txt",
    }
    required_donors = tuple(protocol["incumbent_contract"]["real_donors"])
    if tuple(donor_paths) != required_donors:
        raise ContractError("production donor order drifted")
    donor_sets: dict[str, set[str]] = {}
    for name, path in donor_paths.items():
        require_local_input_hash(source_manifest, path)
        donor_sets[name] = set(read_lines(path))
    perturbmatch_path = ROOT / "data/external/perturbmatch/perts.txt"
    require_local_input_hash(source_manifest, perturbmatch_path)
    perturbmatch_targets = set(read_lines(perturbmatch_path))
    expected_carriers = np.asarray(
        [
            sum(target in donor_sets[name] for name in required_donors)
            for target in targets
        ],
        dtype=np.int16,
    )
    expected_fallback = np.asarray(
        [
            target not in donor_sets["k562_gwps"] and target not in perturbmatch_targets
            for target in targets
        ],
        dtype=bool,
    )
    if not np.array_equal(fallback, expected_fallback):
        raise ContractError("frozen fallback list disagrees with source target lists")
    if not np.array_equal(
        np.asarray(incumbent["incumbent_fallback"], dtype=bool), fallback
    ):
        raise ContractError("incumbent fallback flags differ from frozen eligibility")
    if not np.all(expected == 3):
        raise ContractError("incumbent expected_source_count must be exactly three")
    if not np.array_equal(carriers, expected_carriers):
        raise ContractError("incumbent carrier counts differ from source target lists")
    if not np.array_equal((carriers < expected) & ~fallback, low):
        raise ContractError("incumbent carrier counts differ from frozen eligibility")
    observed = derive_surface_eligibility(incumbent)
    exact = {
        "incumbent_uncovered_targets": fallback_targets,
        "low_confidence_targets": low_targets,
        "eligible_targets": eligible_targets,
        "incumbent_uncovered_ordered_sha256": eligibility[
            "incumbent_uncovered_ordered_sha256"
        ],
        "low_confidence_ordered_sha256": eligibility["low_confidence_ordered_sha256"],
        "eligible_ordered_sha256": eligibility["eligible_ordered_sha256"],
    }
    if observed != exact:
        raise ContractError("incumbent eligibility differs from exact preregistration")
    return observed


def validate_input_manifest(
    value: Mapping[str, Any],
    *,
    role: str,
    approved: Mapping[str, Any],
) -> None:
    if value.get("schema_version") != 1:
        raise ContractError("unsupported input-manifest schema")
    if value.get("experiment") != "114-depmap-coessentiality":
        raise ContractError("input manifest belongs to another experiment")
    if value.get("contains_held_truth") is not False:
        raise ContractError(
            "freeze input manifest must declare contains_held_truth=false"
        )
    if value.get("source_only") is not True:
        raise ContractError("freeze input manifest must declare source_only=true")
    if value.get("role") != role or value.get("surface") != approved.get("surface"):
        raise ContractError("freeze input manifest role or surface drifted")
    if value.get("sealed_before_held_truth") is not True:
        raise ContractError("freeze input manifest was not sealed before held truth")
    if value.get("stage") != "source_only_role_inputs":
        raise ContractError("freeze input manifest stage drifted")
    if set(value) != {
        "schema_version",
        "experiment",
        "stage",
        "role",
        "surface",
        "role_directory",
        "contains_held_truth",
        "source_only",
        "sealed_before_held_truth",
        "role_source_manifest",
        "incumbent",
        "source_bank",
    }:
        raise ContractError("freeze input manifest has missing or unreviewed fields")


def validate_source_id_policy(
    source_bank: SourceBank,
    source_receipt: Mapping[str, Any],
    approved: Mapping[str, Any],
) -> None:
    allowed_source_ids = tuple(map(str, approved["allowed_source_ids"]))
    forbidden_source_ids = set(map(str, approved["forbidden_source_ids"]))
    held_surface_source_id = str(approved["held_surface_source_id"])
    observed_source_ids = set(source_bank.source_ids.tolist())
    if held_surface_source_id not in forbidden_source_ids:
        raise ContractError("held-surface source is not explicitly forbidden")
    if observed_source_ids != set(allowed_source_ids):
        raise ContractError("source-bank source IDs differ from the role allowlist")
    if observed_source_ids & forbidden_source_ids:
        raise ContractError("source bank contains a forbidden held-surface source")
    if source_receipt.get("source_ids") != list(allowed_source_ids):
        raise ContractError("source-bank receipt source IDs drifted")
    if source_receipt.get("forbidden_source_ids") != approved["forbidden_source_ids"]:
        raise ContractError("source-bank receipt forbidden-source policy drifted")
    if source_receipt.get("held_surface_source_id") != held_surface_source_id:
        raise ContractError("source-bank receipt held-surface exclusion drifted")


def load_validated_role_artifacts(
    manifest: Mapping[str, Any],
    *,
    role: str,
    approved: Mapping[str, Any],
    role_directory: Path,
    canonical_targets: np.ndarray,
    canonical_genes: np.ndarray,
    protocol: Mapping[str, Any],
    source_manifest: Mapping[str, Any],
) -> tuple[dict[str, np.ndarray], SourceBank, dict[str, Any]]:
    incumbent_record = manifest["incumbent"]
    source_record = manifest["source_bank"]
    incumbent_path = verify_fixed_record(
        incumbent_record["artifact"],
        f"{role} incumbent artifact",
        required_under=role_directory,
    )
    incumbent_receipt_path = verify_fixed_record(
        incumbent_record["receipt"],
        f"{role} incumbent receipt",
        required_under=role_directory,
    )
    source_path = verify_fixed_record(
        source_record["artifact"],
        f"{role} source-bank artifact",
        required_under=role_directory,
    )
    source_receipt_path = verify_fixed_record(
        source_record["receipt"],
        f"{role} source-bank receipt",
        required_under=role_directory,
    )
    incumbent_receipt = load_json(incumbent_receipt_path)
    source_receipt = load_json(source_receipt_path)
    surface = str(approved["surface"])
    validate_role_receipt(
        incumbent_receipt,
        stage="role_incumbent_source_only",
        role=role,
        surface=surface,
        artifact_path=incumbent_path,
        artifact_sha256=incumbent_record["artifact"]["sha256"],
    )
    validate_role_receipt(
        source_receipt,
        stage="role_source_bank_source_only",
        role=role,
        surface=surface,
        artifact_path=source_path,
        artifact_sha256=source_record["artifact"]["sha256"],
    )
    if incumbent_receipt.get("endpoint") != approved["incumbent_endpoint"]:
        raise ContractError("incumbent endpoint differs from approved role")
    if source_receipt.get("effect_space") != approved["source_effect_space"]:
        raise ContractError("source-bank effect space differs from approved role")

    incumbent = load_incumbent(incumbent_path)
    source_bank = load_source_bank(source_path)
    incumbent_array_names = (
        "targets",
        "genes",
        "effects",
        "incumbent_fallback",
        "source_carrier_count",
        "expected_source_count",
    )
    source_array_names = ("targets", "genes", "effects", "covered", "source_ids")
    validate_array_hashes(
        incumbent_receipt, incumbent, incumbent_array_names, f"{role} incumbent"
    )
    validate_array_hashes(
        source_receipt,
        {
            "targets": source_bank.targets,
            "genes": source_bank.genes,
            "effects": source_bank.effects,
            "covered": source_bank.covered,
            "source_ids": source_bank.source_ids,
        },
        source_array_names,
        f"{role} source bank",
    )

    if not np.array_equal(incumbent["targets"], canonical_targets):
        raise ContractError(
            "incumbent target axis is not the canonical 300-target panel"
        )
    if not np.array_equal(incumbent["genes"], canonical_genes):
        raise ContractError("incumbent gene axis is not canonical")
    if not np.array_equal(source_bank.genes, canonical_genes):
        raise ContractError("source-bank gene axis is not canonical")
    if not np.isfinite(np.asarray(incumbent["effects"], dtype=np.float64)).all():
        raise ContractError("incumbent effects contain nonfinite values")

    panel_positions = np.isin(canonical_genes, canonical_targets)
    if int(panel_positions.sum()) != len(canonical_targets):
        raise ContractError("canonical gene axis does not contain every panel target")
    if np.any(incumbent["effects"][:, panel_positions] != 0):
        raise ContractError("incumbent retains VCC target response coordinates")
    if np.any(source_bank.effects[:, panel_positions] != 0):
        raise ContractError("source bank retains VCC target response coordinates")
    if np.any(source_bank.covered[:, panel_positions]):
        raise ContractError("source bank marks VCC target response coordinates covered")

    eligibility = validate_exact_eligibility(incumbent, protocol, source_manifest)
    validate_source_id_policy(source_bank, source_receipt, approved)
    return (
        incumbent,
        source_bank,
        {
            "eligibility": eligibility,
            "incumbent_receipt_sha256": sha256_file(incumbent_receipt_path),
            "source_bank_receipt_sha256": sha256_file(source_receipt_path),
        },
    )


def freeze_from_manifest(
    input_manifest_path: Path,
    output_dir: Path,
    *,
    protocol: Mapping[str, Any],
) -> dict[str, Any]:
    assert_authorized(protocol, "source_candidate_freeze")
    boundary = load_freeze_boundary()
    if boundary.get("source_candidate_freeze_authorized") is not True:
        raise ContractError("source_candidate_freeze=false in the reviewed boundary")
    role, approved = approved_role_for_manifest(input_manifest_path, boundary)
    input_manifest = load_json(input_manifest_path)
    validate_input_manifest(input_manifest, role=role, approved=approved)
    role_directory = validate_role_manifest(
        input_manifest,
        role=role,
        approved=approved,
    )
    output_dir = require_under_experiment(output_dir)
    expected_output = resolve_repository_path(
        str(approved["output_directory"]), "freeze output"
    )
    if output_dir != expected_output:
        raise ContractError("freeze output path differs from approved role")
    if output_dir.exists():
        raise ContractError(f"freeze output already exists: {output_dir}")

    source_manifest = load_json(SOURCE_MANIFEST_PATH)
    validate_source_manifest(source_manifest)
    canonical_targets, canonical_genes = validate_canonical_axes(
        boundary, source_manifest
    )
    depmap = validate_depmap_feature_artifact(boundary, protocol, source_manifest)
    string = validate_string_feature_artifact(boundary, source_manifest)
    incumbent, source_bank, role_audit = load_validated_role_artifacts(
        input_manifest,
        role=role,
        approved=approved,
        role_directory=role_directory,
        canonical_targets=canonical_targets,
        canonical_genes=canonical_genes,
        protocol=protocol,
        source_manifest=source_manifest,
    )
    consensus = build_consensus(source_bank)
    if not np.array_equal(incumbent["genes"], consensus.genes):
        raise ContractError("incumbent/source axes differ")
    eligibility = role_audit["eligibility"]
    eligible = set(eligibility["eligible_targets"])
    fallback = eligibility["incumbent_uncovered_targets"]
    low = eligibility["low_confidence_targets"]
    panel_targets = set(canonical_targets.tolist())
    blended = blend_features(string, depmap, 1.0)
    feature_by_arm = {
        "string_neighbor": string,
        "depmap_neighbor": depmap,
        PRIMARY_ARM: blended,
    }
    effects: dict[str, np.ndarray] = {
        "identity": incumbent["effects"].copy(),
        "current_fallback": incumbent["effects"].copy(),
    }
    confidence: dict[str, np.ndarray] = {
        "identity": np.abs(np.asarray(incumbent["effects"], dtype=np.float64)),
        "current_fallback": np.abs(np.asarray(incumbent["effects"], dtype=np.float64)),
    }
    arm_audit: dict[str, Any] = {}
    for arm, features in feature_by_arm.items():
        effects[arm], confidence[arm], arm_audit[arm] = build_arm(
            incumbent_targets=incumbent["targets"],
            incumbent_genes=incumbent["genes"],
            incumbent_effects=incumbent["effects"],
            consensus=consensus,
            features=features,
            eligible_targets=eligible,
            panel_targets=panel_targets,
        )
    if not np.array_equal(effects["identity"], effects["current_fallback"]):
        raise AssertionError("identity/current_fallback differ")

    eligible_rows = np.flatnonzero(
        np.isin(incumbent["targets"], np.asarray(eligibility["eligible_targets"]))
    )
    null_rows = np.empty(
        (NULL_DRAWS, len(eligible_rows), len(incumbent["genes"])),
        dtype=incumbent["effects"].dtype,
    )
    null_confidence = np.empty_like(null_rows, dtype=np.float64)
    query_feature_targets = np.empty(
        (NULL_DRAWS, len(eligible_rows)),
        dtype=f"U{max(map(len, eligibility['eligible_targets']))}",
    )
    null_mappings: list[dict[str, str]] = []
    for draw in range(NULL_DRAWS):
        mapping = permuted_query_labels(fallback, low, draw)
        null_effect, null_conf, _ = build_arm(
            incumbent_targets=incumbent["targets"],
            incumbent_genes=incumbent["genes"],
            incumbent_effects=incumbent["effects"],
            consensus=consensus,
            features=blended,
            eligible_targets=eligible,
            panel_targets=panel_targets,
            query_feature_target=mapping,
        )
        null_rows[draw] = null_effect[eligible_rows]
        null_confidence[draw] = null_conf[eligible_rows]
        query_feature_targets[draw] = [
            mapping[target] for target in incumbent["targets"][eligible_rows]
        ]
        null_mappings.append(mapping)

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent)
    )
    try:
        arms_path = staging / "arms.npz"
        nulls_path = staging / "permutation_nulls.npz"
        atomic_npz(
            arms_path,
            targets=incumbent["targets"],
            genes=incumbent["genes"],
            **{f"effect__{name}": effects[name] for name in ARMS},
            **{f"confidence__{name}": confidence[name] for name in ARMS},
        )
        atomic_npz(
            nulls_path,
            eligible_rows=eligible_rows,
            effects=null_rows,
            confidence=null_confidence,
            query_feature_targets=query_feature_targets,
        )
        frozen_manifest = {
            "schema_version": 1,
            "experiment": "114-depmap-coessentiality",
            "stage": "source_only_freeze",
            "contains_held_truth": False,
            "role": role,
            "surface": input_manifest["surface"],
            "protocol_sha256": sha256_file(PROTOCOL_PATH),
            "source_manifest_sha256": sha256_file(SOURCE_MANIFEST_PATH),
            "source_schema_amendment_sha256": sha256_file(SOURCE_SCHEMA_AMENDMENT_PATH),
            "freeze_boundary_sha256": sha256_file(
                SOURCE_CANDIDATE_FREEZE_BOUNDARY_PATH
            ),
            "implementation_sha256": sha256_file(Path(__file__)),
            "input_manifest_sha256": sha256_file(input_manifest_path),
            "panel_targets_ordered_sha256": sha256_lines(canonical_targets),
            "gene_axis_ordered_sha256": sha256_lines(canonical_genes),
            "inputs": {
                "role_manifest": sha256_file(input_manifest_path),
                "role_source_manifest": EXPECTED_ROLE_SOURCE_MANIFEST_SHA256,
                "incumbent_artifact": input_manifest["incumbent"]["artifact"]["sha256"],
                "incumbent_receipt": role_audit["incumbent_receipt_sha256"],
                "source_bank_artifact": input_manifest["source_bank"]["artifact"][
                    "sha256"
                ],
                "source_bank_receipt": role_audit["source_bank_receipt_sha256"],
                "depmap_feature_artifact": boundary["features"]["depmap"][
                    "artifact_sha256"
                ],
                "depmap_feature_receipt": boundary["features"]["depmap"][
                    "receipt_sha256"
                ],
                "string_feature_artifact": boundary["features"]["string"][
                    "artifact_sha256"
                ],
            },
            "eligibility": eligibility,
            "arms": {
                name: {
                    "effect_sha256": sha256_array(effects[name]),
                    "confidence_sha256": sha256_array(confidence[name]),
                }
                for name in ARMS
            },
            "arm_audit": arm_audit,
            "identity_current_fallback_bit_identical": bool(
                np.array_equal(effects["identity"], effects["current_fallback"])
            ),
            "permutation_draws": NULL_DRAWS,
            "permutation_mapping_sha256": sha256_bytes(
                json.dumps(
                    null_mappings, sort_keys=True, separators=(",", ":")
                ).encode()
            ),
            "permutation_query_feature_targets_sha256": sha256_array(
                query_feature_targets
            ),
            "permutation_effects_sha256": sha256_array(null_rows),
            "permutation_confidence_sha256": sha256_array(null_confidence),
            "outputs": {
                "arms.npz": sha256_file(arms_path),
                "permutation_nulls.npz": sha256_file(nulls_path),
            },
            "sealed_before_held_truth": True,
        }
        atomic_json(staging / "manifest.json", frozen_manifest)
        os.replace(staging, output_dir)
        return frozen_manifest
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def load_frozen_directory(path: Path) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    directory = require_under_experiment(path)
    manifest_path = directory / "manifest.json"
    arms_path = directory / "arms.npz"
    nulls_path = directory / "permutation_nulls.npz"
    manifest = load_json(manifest_path)
    protocol = load_json(PROTOCOL_PATH)
    source_manifest = load_json(SOURCE_MANIFEST_PATH)
    boundary = load_freeze_boundary()
    role = manifest.get("role")
    approved = boundary.get("approved_roles", {}).get(role)
    if not isinstance(role, str) or not isinstance(approved, dict):
        raise ContractError("frozen role is not approved by the freeze boundary")
    role_manifest_path = resolve_repository_path(
        str(approved.get("manifest_path", "")), "role manifest"
    )
    if sha256_file(role_manifest_path) != approved.get("manifest_sha256"):
        raise ContractError("approved role manifest changed")
    role_manifest = load_json(role_manifest_path)
    validate_input_manifest(role_manifest, role=role, approved=approved)
    role_directory = validate_role_manifest(
        role_manifest,
        role=role,
        approved=approved,
    )
    canonical_targets, canonical_genes = validate_canonical_axes(
        boundary, source_manifest
    )
    validate_depmap_feature_artifact(boundary, protocol, source_manifest)
    validate_string_feature_artifact(boundary, source_manifest)
    _, _, role_audit = load_validated_role_artifacts(
        role_manifest,
        role=role,
        approved=approved,
        role_directory=role_directory,
        canonical_targets=canonical_targets,
        canonical_genes=canonical_genes,
        protocol=protocol,
        source_manifest=source_manifest,
    )
    expected = {
        "experiment": "114-depmap-coessentiality",
        "stage": "source_only_freeze",
        "contains_held_truth": False,
        "sealed_before_held_truth": True,
        "role": role,
        "surface": approved["surface"],
        "protocol_sha256": sha256_file(PROTOCOL_PATH),
        "source_manifest_sha256": sha256_file(SOURCE_MANIFEST_PATH),
        "source_schema_amendment_sha256": sha256_file(SOURCE_SCHEMA_AMENDMENT_PATH),
        "freeze_boundary_sha256": sha256_file(SOURCE_CANDIDATE_FREEZE_BOUNDARY_PATH),
        "input_manifest_sha256": sha256_file(role_manifest_path),
        "implementation_sha256": sha256_file(Path(__file__)),
        "permutation_draws": NULL_DRAWS,
    }
    for key, value in expected.items():
        if manifest.get(key) != value:
            raise ContractError(f"frozen manifest {key} drifted")
    for name, artifact in (
        ("arms.npz", arms_path),
        ("permutation_nulls.npz", nulls_path),
    ):
        if not artifact.is_file() or sha256_file(artifact) != manifest["outputs"].get(
            name
        ):
            raise ContractError(f"frozen artifact changed: {name}")
    with np.load(arms_path, allow_pickle=False) as archive:
        required = {
            "targets",
            "genes",
            *(f"effect__{name}" for name in ARMS),
            *(f"confidence__{name}" for name in ARMS),
        }
        if set(archive.files) != required:
            raise ContractError("frozen arm keys drifted")
        arrays = {name: archive[name] for name in required}
    with np.load(nulls_path, allow_pickle=False) as archive:
        if set(archive.files) != {
            "eligible_rows",
            "effects",
            "confidence",
            "query_feature_targets",
        }:
            raise ContractError("frozen null keys drifted")
        arrays.update(
            {
                "null_eligible_rows": archive["eligible_rows"],
                "null_effects": archive["effects"],
                "null_confidence": archive["confidence"],
                "null_query_feature_targets": archive["query_feature_targets"],
            }
        )
    targets = arrays["targets"].astype(str)
    genes = arrays["genes"].astype(str)
    if not np.array_equal(targets, canonical_targets):
        raise ContractError("frozen target axis is not canonical")
    if not np.array_equal(genes, canonical_genes):
        raise ContractError("frozen gene axis is not canonical")
    if manifest.get("eligibility") != role_audit["eligibility"]:
        raise ContractError("frozen eligibility differs from the approved role")
    shape = (len(targets), len(genes))
    for name in ARMS:
        if arrays[f"effect__{name}"].shape != shape:
            raise ContractError(f"frozen effect shape drifted for {name}")
        if arrays[f"confidence__{name}"].shape != shape:
            raise ContractError(f"frozen confidence shape drifted for {name}")
        if (
            sha256_array(arrays[f"effect__{name}"])
            != manifest["arms"][name]["effect_sha256"]
        ):
            raise ContractError(f"frozen effect hash drifted for {name}")
        if (
            sha256_array(arrays[f"confidence__{name}"])
            != manifest["arms"][name]["confidence_sha256"]
        ):
            raise ContractError(f"frozen confidence hash drifted for {name}")
        expected_confidence = np.abs(
            np.asarray(arrays[f"effect__{name}"], dtype=np.float64)
        )
        if not np.array_equal(arrays[f"confidence__{name}"], expected_confidence):
            raise ContractError(f"frozen confidence is not abs(effect) for {name}")
    if not np.array_equal(
        arrays["effect__identity"], arrays["effect__current_fallback"]
    ):
        raise ContractError("frozen identity and current_fallback differ")
    eligible_rows = np.asarray(arrays["null_eligible_rows"], dtype=np.int64)
    expected_eligible_rows = np.flatnonzero(
        np.isin(targets, np.asarray(protocol["eligibility"]["eligible_targets"]))
    )
    if not np.array_equal(eligible_rows, expected_eligible_rows):
        raise ContractError("frozen null eligible rows differ from preregistration")
    if arrays["null_effects"].shape != (
        NULL_DRAWS,
        len(eligible_rows),
        len(genes),
    ):
        raise ContractError("frozen null effect shape drifted")
    if arrays["null_confidence"].shape != arrays["null_effects"].shape:
        raise ContractError("frozen null confidence shape drifted")
    if arrays["null_query_feature_targets"].shape != (
        NULL_DRAWS,
        len(eligible_rows),
    ):
        raise ContractError("frozen null query-label shape drifted")
    if sha256_array(arrays["null_effects"]) != manifest["permutation_effects_sha256"]:
        raise ContractError("frozen null effect hash drifted")
    if (
        sha256_array(arrays["null_confidence"])
        != manifest["permutation_confidence_sha256"]
    ):
        raise ContractError("frozen null confidence hash drifted")
    if (
        sha256_array(arrays["null_query_feature_targets"])
        != manifest["permutation_query_feature_targets_sha256"]
    ):
        raise ContractError("frozen null query-label hash drifted")
    if not np.array_equal(
        arrays["null_confidence"],
        np.abs(np.asarray(arrays["null_effects"], dtype=np.float64)),
    ):
        raise ContractError("frozen null confidence is not abs(effect)")
    fallback = protocol["eligibility"]["incumbent_uncovered_targets"]
    low = protocol["eligibility"]["low_confidence_targets"]
    regenerated_mappings = [
        permuted_query_labels(fallback, low, draw) for draw in range(NULL_DRAWS)
    ]
    if manifest.get("permutation_mapping_sha256") != sha256_bytes(
        json.dumps(regenerated_mappings, sort_keys=True, separators=(",", ":")).encode()
    ):
        raise ContractError("frozen null mapping hash does not regenerate exactly")
    expected_query_targets = np.asarray(
        [
            [mapping[target] for target in targets[eligible_rows]]
            for mapping in regenerated_mappings
        ]
    )
    if not np.array_equal(
        arrays["null_query_feature_targets"].astype(str), expected_query_targets
    ):
        raise ContractError("frozen null query labels do not regenerate exactly")
    return manifest, arrays


def validate_truth_manifest(
    value: Mapping[str, Any], frozen: Mapping[str, Any]
) -> Path:
    if value.get("schema_version") != 1:
        raise ContractError("unsupported truth-manifest schema")
    if value.get("experiment") != "114-depmap-coessentiality":
        raise ContractError("truth manifest belongs to another experiment")
    if value.get("contains_held_truth") is not True:
        raise ContractError("truth manifest must declare contains_held_truth=true")
    if value.get("surface") != frozen.get("surface"):
        raise ContractError("truth and frozen surfaces differ")
    if value.get("clean") is not True:
        raise ContractError("only a clean held surface may be scored")
    if value.get("direction_valid") is not True:
        raise ContractError("directional precision requires direction_valid=true")
    truth_record = value.get("truth_npz")
    if not isinstance(truth_record, dict):
        raise ContractError("truth manifest lacks truth_npz")
    return verify_path_record(truth_record, "held truth")


def score_from_manifest(
    frozen_dir: Path,
    truth_manifest_path: Path,
    output_path: Path,
    *,
    protocol: Mapping[str, Any],
) -> dict[str, Any]:
    assert_authorized(protocol, "held_truth_access")
    assert_authorized(protocol, "evaluation")

    # Authenticate the complete source-only freeze before accepting or opening
    # any held-truth path.
    frozen, arrays = load_frozen_directory(frozen_dir)
    truth_manifest = load_json(truth_manifest_path)
    truth_path = validate_truth_manifest(truth_manifest, frozen)
    with np.load(truth_path, allow_pickle=False) as archive:
        required = {"targets", "genes", "effects", "direction_sign"}
        if set(archive.files) != required:
            raise ContractError(f"truth NPZ keys must be {sorted(required)}")
        truth_targets = archive["targets"].astype(str)
        truth_genes = archive["genes"].astype(str)
        truth_effect = np.asarray(archive["effects"], dtype=np.float64)
        truth_sign = np.asarray(archive["direction_sign"], dtype=np.int8)

    targets = arrays["targets"].astype(str)
    genes = arrays["genes"].astype(str)
    if not np.array_equal(truth_targets, targets):
        raise ContractError("held truth target order differs from frozen predictions")
    if not np.array_equal(truth_genes, genes):
        raise ContractError("held truth gene order differs from frozen predictions")
    if truth_effect.shape != (len(targets), len(genes)):
        raise ContractError("held truth effect shape is invalid")
    if truth_sign.shape != truth_effect.shape:
        raise ContractError("held truth direction-sign shape is invalid")
    if not np.isfinite(truth_effect).all():
        raise ContractError("held truth effects contain nonfinite values")
    if not np.isin(truth_sign, [-1, 0, 1]).all():
        raise ContractError("held truth direction signs must be -1, 0, or 1")

    source_manifest = load_json(SOURCE_MANIFEST_PATH)
    panel_path = ROOT / "data/challenge/pert_counts.csv"
    expected_panel_hash = source_manifest["local_inputs"][
        "data/challenge/pert_counts.csv"
    ]
    if sha256_file(panel_path) != expected_panel_hash:
        raise ContractError("challenge panel file drifted")
    panel = read_csv_column(panel_path, "target_gene")
    if sha256_lines(panel) != frozen["panel_targets_ordered_sha256"]:
        raise ContractError("frozen response-coordinate panel drifted")
    panel_set = set(panel)
    allowed = np.asarray([gene not in panel_set for gene in genes], dtype=bool)

    arms = {name: arrays[f"effect__{name}"] for name in ARMS}
    confidence = {name: arrays[f"confidence__{name}"] for name in ARMS}
    real = evaluate_surface(arms, confidence, truth_effect, truth_sign, allowed)
    fallback_pds = float(real["current_fallback"]["raw_pds"])
    eligible_rows = np.asarray(arrays["null_eligible_rows"], dtype=np.int64)
    null_deltas = np.empty(NULL_DRAWS, dtype=np.float64)
    for draw in range(NULL_DRAWS):
        null = np.asarray(arms["current_fallback"]).copy()
        null[eligible_rows] = arrays["null_effects"][draw]
        null_pds, _ = raw_uncentered_pds(null, truth_effect, allowed)
        null_deltas[draw] = null_pds - fallback_pds

    result = {
        "schema_version": 1,
        "experiment": "114-depmap-coessentiality",
        "stage": "held_surface_score",
        "surface": frozen["surface"],
        "clean": True,
        "frozen_manifest_sha256": sha256_file(
            require_under_experiment(frozen_dir) / "manifest.json"
        ),
        "truth_manifest_sha256": sha256_file(truth_manifest_path),
        "truth_npz_sha256": sha256_file(truth_path),
        "response_coordinates": {
            "allowed": int(allowed.sum()),
            "removed_vcc_targets": int((~allowed).sum()),
            "allowed_mask_sha256": sha256_array(allowed.astype(np.uint8)),
        },
        "arms": real,
        "primary_delta": {
            "raw_pds": float(real[PRIMARY_ARM]["raw_pds"]) - fallback_pds,
            "directional_precision": (
                float(real[PRIMARY_ARM]["directional_precision"])
                - float(real["current_fallback"]["directional_precision"])
            ),
        },
        "permutation_macro_delta_not_applicable_single_surface": True,
        "permutation_raw_pds_deltas": null_deltas.tolist(),
        "permutation_p95": float(np.quantile(null_deltas, 0.95, method="higher")),
        "no_post_truth_fit": True,
    }
    atomic_json(output_path, result)
    return result


def run_manifest() -> dict[str, Any]:
    protocol = load_json(PROTOCOL_PATH)
    source = load_json(SOURCE_MANIFEST_PATH)
    validate_protocol(protocol)
    validate_source_manifest(source)
    boundary = load_freeze_boundary()
    depmap = validate_depmap_feature_artifact(boundary, protocol, source)
    string = validate_string_feature_artifact(boundary, source)
    audit = derive_production_eligibility(
        protocol,
        source,
        panel_path=ROOT / "data/challenge/pert_counts.csv",
        incumbent_meta_path=ROOT / "submissions/a15-cd4f010.vcc.meta.json",
        donor_pert_paths={
            "k562_gwps": ROOT / "data/donors/k562_gwps/perts.txt",
            "orion_hct116": ROOT / "data/donors/orion_hct116/perts.txt",
            "orion_hek293t": ROOT / "data/donors/orion_hek293t/perts.txt",
        },
        perturbmatch_pert_path=ROOT / "data/external/perturbmatch/perts.txt",
    )
    authorization = protocol["current_authorization"]
    audit.update(
        {
            "protocol_sha256": sha256_file(PROTOCOL_PATH),
            "source_manifest_sha256": sha256_file(SOURCE_MANIFEST_PATH),
            "source_schema_amendment_sha256": sha256_file(SOURCE_SCHEMA_AMENDMENT_PATH),
            "source_candidate_freeze_boundary_sha256": sha256_file(
                SOURCE_CANDIDATE_FREEZE_BOUNDARY_PATH
            ),
            "role_source_manifest_sha256": sha256_file(ROLE_SOURCE_MANIFEST_PATH),
            "depmap_expected_shape": list(EXPECTED_DEPMAP_SHAPE),
            "depmap_feature_shape": list(depmap.matrix.shape),
            "depmap_feature_artifact_sha256": sha256_file(DEPMAP_FEATURES_PATH),
            "depmap_feature_receipt_sha256": sha256_file(DEPMAP_FEATURES_RECEIPT_PATH),
            "string_feature_shape": list(string.matrix.shape),
            "string_feature_artifact_sha256": sha256_file(STRING_FEATURES_PATH),
            "approved_role_count": len(boundary["approved_roles"]),
            "freeze_boundary_authorized": boundary[
                "source_candidate_freeze_authorized"
            ],
            "depmap_large_download_authorized": authorization["depmap_large_download"],
            "source_feature_build_authorized": authorization["source_feature_build"],
            "source_candidate_freeze_authorized": authorization[
                "source_candidate_freeze"
            ],
            "held_truth_access_authorized": authorization["held_truth_access"],
            "evaluation_authorized": authorization["evaluation"],
        }
    )
    return audit


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("manifest", help="read-only frozen-contract audit")
    download = subparsers.add_parser(
        "download-depmap",
        help="download the one allowlisted DepMap source into the frozen work path",
    )
    download.add_argument("--confirm-large-download", action="store_true")
    build = subparsers.add_parser(
        "build-depmap",
        help="build the exact, provenance-receipted 24Q4 source features",
    )
    build.add_argument("--csv", type=Path, default=DEPMAP_CSV_PATH)
    build.add_argument("--out", type=Path, default=DEPMAP_FEATURES_PATH)
    freeze = subparsers.add_parser(
        "freeze",
        help="reserved source-only candidate freeze; currently locked",
    )
    freeze.add_argument("--input-manifest", type=Path, required=True)
    freeze.add_argument("--out", type=Path, required=True)
    score = subparsers.add_parser(
        "score",
        help="reserved held-truth evaluation; currently locked",
    )
    score.add_argument("--frozen", type=Path, required=True)
    score.add_argument("--truth-manifest", type=Path, required=True)
    score.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    protocol = load_json(PROTOCOL_PATH)
    source = load_json(SOURCE_MANIFEST_PATH)
    validate_protocol(protocol)
    validate_source_manifest(source)
    if args.command == "manifest":
        print(json.dumps(run_manifest(), indent=2, sort_keys=True))
        return
    if args.command == "download-depmap":
        receipt = download_depmap(
            protocol=protocol,
            source_manifest=source,
            confirmed=args.confirm_large_download,
        )
        print(json.dumps(receipt, indent=2, sort_keys=True))
        return
    if args.command == "build-depmap":
        receipt = build_depmap_features(
            args.csv,
            args.out,
            protocol=protocol,
            source_manifest=source,
        )
        print(json.dumps(receipt, indent=2, sort_keys=True))
        return
    if args.command == "freeze":
        result = freeze_from_manifest(
            args.input_manifest,
            args.out,
            protocol=protocol,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
        return
    if args.command == "score":
        result = score_from_manifest(
            args.frozen,
            args.truth_manifest,
            args.out,
            protocol=protocol,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
        return
    raise AssertionError(args.command)


if __name__ == "__main__":
    main()
