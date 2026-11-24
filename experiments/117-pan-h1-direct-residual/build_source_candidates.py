#!/usr/bin/env python3
"""Build source-only Pan2025 direct-residual candidates for frozen a14 roles.

The reviewed Experiment 115 compact donor and this implementation are pinned.
The build remains source-only: held truth, evaluation, counts, and submission
stay prohibited. Array-level functions remain usable without any outcome.
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
import zipfile
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
SOURCE_MANIFEST_PATH = HERE / "SOURCE_INPUT_MANIFEST.json"
FREEZE_BOUNDARY_PATH = HERE / "SOURCE_CANDIDATE_FREEZE_BOUNDARY.json"
DEFAULT_OUTPUT_ROOT = HERE / "work" / "source-candidates"

EXPECTED_PROTOCOL_SHA256 = (
    "2adf7f9e25cd98653a523292cd57538dac7daa441c68b12e7458983fa836326b"
)
EXPECTED_SOURCE_MANIFEST_SHA256 = (
    "d237df1b48d1448d1c2a3b361e15ed5a6937ca7a9fd326fa4e3560db4d6de61d"
)

EXPERIMENT = "117-pan-h1-direct-residual"
ROLE_STATUS = "frozen_before_experiment_107_nadig_outcome_access"
FRACTIONS = (0.05, 0.10, 0.20)
PRIMARY_FRACTION = 0.10
PERMUTATION_REPLICATES = 99
PERMUTATION_NAMESPACE = "experiment-117-pan-direct-target-null-v1"
BOOTSTRAP_REPLICATES = 20_000
BOOTSTRAP_NAMESPACE = "experiment-117-pan-direct-cluster-bootstrap-v1"
NADIG_ROLES = tuple(f"nadig.gallery_{index}" for index in range(1, 7))
CROSS_CONTEXT_ROLES = ("orion_hek293t", *NADIG_ROLES)
FIXED_ZIP_TIMESTAMP = (1980, 1, 1, 0, 0, 0)
GAMMA = 0.65
LAMBDA = 0.0
OFFICIAL_GENE_COUNT = 18_533
PANEL_TARGET_COUNT = 300
PAN_TARGET_COUNT = 2_354
PAN_PANEL_COUNT = 64
PAN_PANEL_NEW_COUNT = 57
COMPACT_ARTIFACT_FILES = {
    "b_ctrl.npy",
    "delta.npy",
    "genes.txt",
    "measured_mask.npy",
    "meta.json",
    "perts.txt",
    "rows.tsv",
    "source_covered_gene_mask.npy",
    "source_gene_map.tsv",
    "vcc_target_mask.npy",
    "vcc_targets.txt",
}
REQUIRED_COMPACT_FILES = COMPACT_ARTIFACT_FILES | {
    "artifact-manifest.json",
    "compact-donor-receipt.json",
    "compact-donor-receipt.sha256",
}


class ContractError(RuntimeError):
    """A source-only provenance, geometry, or freeze invariant failed."""


@dataclass(frozen=True)
class CompactDonor:
    directory: Path
    genes: tuple[str, ...]
    targets: tuple[str, ...]
    delta: np.ndarray
    measured: np.ndarray
    receipt_sha256: str
    artifact_manifest_sha256: str


@dataclass(frozen=True)
class RoleBundle:
    name: str
    targets: np.ndarray
    genes: np.ndarray
    score_indices: np.ndarray
    incumbent: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    operator_columns: np.ndarray
    input_bindings: tuple[dict[str, Any], ...]


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256(
        "".join(f"{value}\n" for value in values).encode("utf-8")
    ).hexdigest()


def sha256_array(values: np.ndarray) -> str:
    """Match the logical-array binding used by Experiment 107."""
    array = np.ascontiguousarray(np.asarray(values))
    if array.dtype.hasobject:
        raise ContractError("object arrays cannot be bound")
    digest = hashlib.sha256()
    digest.update(b"exp107-array-v1\0")
    digest.update(array.dtype.str.encode("ascii"))
    digest.update(b"\0")
    digest.update(json.dumps(list(array.shape), separators=(",", ":")).encode("ascii"))
    digest.update(b"\0")
    if array.nbytes:
        digest.update(memoryview(array).cast("B"))
    return digest.hexdigest()


def array_binding(values: np.ndarray) -> dict[str, Any]:
    array = np.ascontiguousarray(np.asarray(values))
    return {
        "dtype": array.dtype.str,
        "shape": list(array.shape),
        "c_order": True,
        "sha256_array": sha256_array(array),
    }


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ContractError(f"cannot read valid JSON: {path}") from error
    if not isinstance(value, dict):
        raise ContractError(f"JSON root must be an object: {path}")
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def _is_sha256(value: object) -> bool:
    text = str(value)
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def deterministic_seed(namespace: str, *parts: object) -> int:
    """Derive one stable PCG64 seed from an explicit SHA-256 namespace."""
    payload = "\0".join((namespace, *(str(part) for part in parts))).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def _safe_repo_file(relative: str, label: str) -> Path:
    path = (ROOT / relative).resolve()
    _require(path.is_relative_to(ROOT.resolve()), f"{label} escapes repository")
    _require(path.is_file() and not path.is_symlink(), f"{label} missing or unsafe")
    return path


def _safe_output(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    _require(
        resolved != HERE and resolved.is_relative_to(HERE),
        f"output must be a child of {HERE}",
    )
    return resolved


def _read_lines(
    path: Path,
    *,
    expected: int | None = None,
    require_final_newline: bool = False,
) -> tuple[str, ...]:
    if require_final_newline:
        _require(path.read_bytes().endswith(b"\n"), f"{path.name} lacks final newline")
    values = tuple(path.read_text(encoding="utf-8").splitlines())
    _require(
        all(value and value == value.strip() for value in values),
        f"{path.name} has invalid labels",
    )
    _require(len(values) == len(set(values)), f"{path.name} has duplicate labels")
    if expected is not None:
        _require(len(values) == expected, f"{path.name} count changed")
    return values


def _read_csv_axis(
    path: Path, column: str, *, expected: int | None = None
) -> tuple[str, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        _require(reader.fieldnames == [column], f"{path.name} header changed")
        values = tuple(str(row[column]) for row in reader)
    _require(len(values) == len(set(values)), f"{path.name} has duplicate labels")
    if expected is not None:
        _require(len(values) == expected, f"{path.name} count changed")
    return values


def validate_scaffold() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Validate immutable declarations without opening any role or outcome array."""
    protocol = load_json(PROTOCOL_PATH)
    source = load_json(SOURCE_MANIFEST_PATH)
    boundary = load_json(FREEZE_BOUNDARY_PATH)
    _require(protocol.get("experiment") == EXPERIMENT, "protocol experiment changed")
    _require(
        source.get("experiment") == EXPERIMENT, "source manifest experiment changed"
    )
    _require(
        boundary.get("experiment") == EXPERIMENT, "freeze boundary experiment changed"
    )
    if EXPECTED_PROTOCOL_SHA256 != "TO_BE_PINNED":
        _require(
            sha256_file(PROTOCOL_PATH) == EXPECTED_PROTOCOL_SHA256,
            "preregistered protocol changed",
        )
    if EXPECTED_SOURCE_MANIFEST_SHA256 != "TO_BE_PINNED":
        _require(
            sha256_file(SOURCE_MANIFEST_PATH) == EXPECTED_SOURCE_MANIFEST_SHA256,
            "source input manifest changed",
        )
    _require(
        source["protocol"]["sha256"] == sha256_file(PROTOCOL_PATH),
        "source manifest does not bind the protocol",
    )
    _require(
        boundary["protocol"]["sha256"] == sha256_file(PROTOCOL_PATH),
        "freeze boundary does not bind the protocol",
    )
    _require(
        boundary["source_manifest"]["sha256"] == sha256_file(SOURCE_MANIFEST_PATH),
        "freeze boundary does not bind the source manifest",
    )
    for label, record in source["canonical_axes"].items():
        path = _safe_repo_file(record["path"], label)
        _require(sha256_file(path) == record["sha256"], f"{label} hash drifted")
        if label == "genes":
            values = _read_csv_axis(path, "gene_name", expected=record["count"])
        elif label == "panel":
            values = _read_csv_axis(path, "target_gene", expected=record["count"])
        else:
            values = _read_lines(path, expected=record["count"])
        if "ordered_sha256" in record:
            _require(
                sha256_lines(values) == record["ordered_sha256"],
                f"{label} ordered axis changed",
            )
    for label in (
        "source_manifest",
        "acquisition_amendment",
        "download_receipt",
        "extractor",
    ):
        record = source["pan_source"][label]
        path = _safe_repo_file(record["path"], f"Pan {label}")
        _require(
            sha256_file(path) == record["sha256"],
            f"Pan {label} hash drifted",
        )
    for label in ("strict_contract", "validator"):
        record = source["upstream_geometry"][label]
        path = _safe_repo_file(record["path"], f"a14 {label}")
        _require(
            sha256_file(path) == record["sha256"],
            f"a14 {label} hash drifted",
        )
    _require(
        source["upstream_geometry"]["a14"]
        == {
            "alpha": 0.982,
            "donor_weights": "equal",
            "gene_norm": False,
            "empirical_bayes": False,
            "perturbmatch_fraction": 0.5,
            "cd4_fraction": 0.05,
            "whiten_gamma": 0.65,
            "whiten_lambda": 0.0,
            "coordinate_exclusion": "VCC_targets_union_role_held_targets",
        },
        "exact a14 geometry declaration changed",
    )
    for role, record in source["roles"].items():
        expected_pan = record.get("expected_pan_overlap")
        expected_novel = record.get("expected_novel_overlap")
        _require(
            isinstance(expected_pan, int)
            and isinstance(expected_novel, int)
            and 0 <= expected_novel <= expected_pan <= record["target_count"],
            f"{role} expected overlap counts are invalid",
        )
        for kind in ("incumbent", "operator"):
            path = _safe_repo_file(record[f"{kind}_manifest"], f"{role} {kind}")
            _require(
                sha256_file(path) == record[f"{kind}_manifest_sha256"],
                f"{role} {kind} manifest hash drifted",
            )
            manifest = load_json(path)
            _require(
                manifest.get("status") == ROLE_STATUS, f"{role} {kind} status changed"
            )
            _require(
                manifest.get("npz", {}).get("sha256") == record[f"{kind}_npz_sha256"],
                f"{role} {kind} NPZ declaration changed",
            )
    return protocol, source, boundary


def assert_real_build_authorized(
    source: Mapping[str, Any], boundary: Mapping[str, Any]
) -> None:
    pan = source["pan_source"]
    missing = [
        key
        for key in (
            "receipt_sha256",
            "artifact_manifest_sha256",
            "delta_sha256",
            "measured_mask_sha256",
        )
        if not _is_sha256(pan.get(key))
    ]
    _require(not missing, f"compact donor hashes are not pinned: {missing}")
    _require(pan.get("real_source_authorized") is True, "real Pan source is locked")
    _require(
        boundary.get("source_candidate_build_authorized") is True,
        "source candidate build is locked",
    )
    _require(
        boundary.get("held_truth_access_authorized") is False,
        "source build boundary must not authorize held truth",
    )
    for key in (
        "evaluation_authorized",
        "cell_eval2_authorized",
        "count_generation_authorized",
        "challenge_build_authorized",
        "submission_authorized",
    ):
        _require(
            boundary.get(key) is False,
            f"source build boundary must keep {key}=false",
        )
    authorization = source.get("authorization", {})
    for key in (
        "held_truth_access",
        "evaluation",
        "cell_eval2",
        "count_generation",
        "submission",
    ):
        _require(
            authorization.get(key) is False,
            f"source manifest must keep {key}=false",
        )
    builder_hash = boundary.get("builder", {}).get("sha256")
    _require(_is_sha256(builder_hash), "builder hash is not pinned")
    _require(
        builder_hash == sha256_file(Path(__file__).resolve()),
        "builder hash differs from freeze boundary",
    )
    roles = source["roles"]
    approved = boundary.get("approved_roles")
    _require(
        isinstance(approved, dict) and set(approved) == set(roles),
        "every role must be independently approved",
    )
    _require(all(value is True for value in approved.values()), "a role is unapproved")


def validate_compact_hash_pins(
    pan: Mapping[str, Any],
    *,
    receipt_sha256: str,
    artifact_manifest_sha256: str,
    artifact_hashes: Mapping[str, str],
) -> None:
    """Fail closed if any independently reviewed compact hash drifts."""
    _require(
        receipt_sha256 == pan.get("receipt_sha256"),
        "compact donor receipt hash drifted",
    )
    _require(
        artifact_manifest_sha256 == pan.get("artifact_manifest_sha256"),
        "compact donor artifact manifest hash drifted",
    )
    _require(
        artifact_hashes.get("delta.npy") == pan.get("delta_sha256"),
        "delta.npy pinned hash drifted",
    )
    _require(
        artifact_hashes.get("measured_mask.npy")
        == pan.get("measured_mask_sha256"),
        "measured_mask.npy pinned hash drifted",
    )


def authenticate_compact_donor(source: Mapping[str, Any]) -> CompactDonor:
    """Authenticate Experiment 115's complete receipt chain before array load."""
    pan = source["pan_source"]
    directory = (ROOT / pan["directory"]).resolve()
    _require(directory.is_dir() and not directory.is_symlink(), "compact donor missing")
    observed = {path.name for path in directory.iterdir() if path.is_file()}
    _require(
        REQUIRED_COMPACT_FILES.issubset(observed),
        "compact donor file set is incomplete",
    )
    receipt_path = directory / "compact-donor-receipt.json"
    checksum_path = directory / "compact-donor-receipt.sha256"
    manifest_path = directory / "artifact-manifest.json"
    receipt_sha256 = sha256_file(receipt_path)
    _require(receipt_sha256 == pan["receipt_sha256"], "compact donor receipt hash drifted")
    checksum = checksum_path.read_text(encoding="utf-8").rstrip("\n").split("  ")
    _require(
        checksum == [pan["receipt_sha256"], "compact-donor-receipt.json"],
        "compact donor receipt checksum file changed",
    )
    receipt = load_json(receipt_path)
    _require(
        receipt.get("status") == pan["required_receipt_status"],
        "compact donor receipt status changed",
    )
    _require(
        receipt.get("axes")
        == {
            "rows": pan["expected_rows"],
            "canonical_genes": pan["expected_genes"],
            "vcc_targets": pan["expected_vcc_targets"],
        },
        "compact donor receipt axes changed",
    )
    _require(
        receipt.get("implementation")
        == {
            "path": pan["extractor"]["path"],
            "sha256": pan["extractor"]["sha256"],
        },
        "compact donor extractor binding changed",
    )
    _require(
        receipt.get("source")
        == {
            "source_manifest_sha256": pan["source_manifest"]["sha256"],
            "acquisition_amendment_sha256": pan["acquisition_amendment"]["sha256"],
            "download_receipt_sha256": pan["download_receipt"]["sha256"],
            "h5ad_sha256": pan["source_asset"]["sha256"],
            "h5ad_size_bytes": pan["source_asset"]["size_bytes"],
        },
        "compact donor upstream source binding changed",
    )
    download_receipt = load_json(ROOT / pan["download_receipt"]["path"])
    _require(
        download_receipt.get("status") == "verified_public_pseudobulk_download",
        "Pan download receipt status changed",
    )
    _require(
        download_receipt.get("source_manifest_sha256")
        == pan["source_manifest"]["sha256"],
        "Pan download receipt source binding changed",
    )
    _require(
        download_receipt.get("acquisition_amendment_sha256")
        == pan["acquisition_amendment"]["sha256"],
        "Pan download receipt acquisition binding changed",
    )
    _require(
        download_receipt.get("asset", {}).get("sha256")
        == pan["source_asset"]["sha256"]
        and download_receipt.get("asset", {}).get("size_bytes")
        == pan["source_asset"]["size_bytes"],
        "Pan download receipt asset binding changed",
    )
    artifact_manifest_sha256 = sha256_file(manifest_path)
    _require(
        artifact_manifest_sha256 == pan["artifact_manifest_sha256"],
        "compact donor artifact manifest hash drifted",
    )
    _require(
        receipt.get("artifact_manifest")
        == {
            "path": "artifact-manifest.json",
            "size_bytes": manifest_path.stat().st_size,
            "sha256": pan["artifact_manifest_sha256"],
        },
        "compact donor receipt-to-manifest binding changed",
    )
    artifact_manifest = load_json(manifest_path)
    files = artifact_manifest.get("files")
    _require(
        isinstance(files, dict) and set(files) == COMPACT_ARTIFACT_FILES,
        "compact artifact inventory changed",
    )
    observed_hashes: dict[str, str] = {}
    for name in sorted(COMPACT_ARTIFACT_FILES):
        path = directory / name
        _require(path.is_file() and not path.is_symlink(), f"{name} missing or unsafe")
        observed_hashes[name] = sha256_file(path)
        expected = files[name]
        _require(
            expected.get("sha256") == observed_hashes[name]
            and expected.get("size_bytes") == path.stat().st_size,
            f"{name} artifact binding changed",
        )
    validate_compact_hash_pins(
        pan,
        receipt_sha256=receipt_sha256,
        artifact_manifest_sha256=artifact_manifest_sha256,
        artifact_hashes=observed_hashes,
    )
    genes = _read_lines(
        directory / "genes.txt",
        expected=OFFICIAL_GENE_COUNT,
        require_final_newline=True,
    )
    targets = _read_lines(
        directory / "perts.txt",
        expected=PAN_TARGET_COUNT,
        require_final_newline=True,
    )
    _require(
        genes
        == _read_csv_axis(
            ROOT / source["canonical_axes"]["genes"]["path"],
            "gene_name",
            expected=OFFICIAL_GENE_COUNT,
        ),
        "compact donor gene axis is not canonical",
    )
    panel_targets = _read_csv_axis(
        ROOT / source["canonical_axes"]["panel"]["path"],
        "target_gene",
        expected=PANEL_TARGET_COUNT,
    )
    expected_vcc = tuple(target for target in targets if target in set(panel_targets))
    observed_vcc = _read_lines(
        directory / "vcc_targets.txt",
        expected=PAN_PANEL_COUNT,
        require_final_newline=True,
    )
    _require(observed_vcc == expected_vcc, "compact donor VCC target axis changed")
    _require(
        sha256_lines(observed_vcc)
        == source["coverage"]["vcc_panel_supported_ordered_sha256"],
        "compact donor VCC target-axis hash changed",
    )
    delta = np.load(directory / "delta.npy", mmap_mode="r", allow_pickle=False)
    measured = np.load(
        directory / "measured_mask.npy", mmap_mode="r", allow_pickle=False
    )
    _require(
        delta.shape == (PAN_TARGET_COUNT, OFFICIAL_GENE_COUNT),
        "Pan delta shape changed",
    )
    _require(delta.dtype == np.dtype("<f4"), "Pan delta dtype changed")
    _require(measured.shape == delta.shape, "Pan measured mask shape changed")
    _require(measured.dtype == np.dtype("|b1"), "Pan measured mask dtype changed")
    vcc_mask = np.load(
        directory / "vcc_target_mask.npy", mmap_mode="r", allow_pickle=False
    )
    _require(vcc_mask.shape == (PAN_TARGET_COUNT,), "Pan VCC mask shape changed")
    _require(vcc_mask.dtype == np.dtype("|b1"), "Pan VCC mask dtype changed")
    _require(
        np.array_equal(
            np.flatnonzero(vcc_mask),
            np.asarray([targets.index(target) for target in observed_vcc]),
        ),
        "Pan VCC mask does not match its target axis",
    )
    meta = load_json(directory / "meta.json")
    _require(
        meta.get("safety")
        == {
            "held_truth_opened": False,
            "model_fit": False,
            "evaluation_run": False,
            "counts_generated": False,
            "submission_built": False,
        },
        "compact donor safety declaration changed",
    )
    return CompactDonor(
        directory=directory,
        genes=genes,
        targets=targets,
        delta=delta,
        measured=measured,
        receipt_sha256=pan["receipt_sha256"],
        artifact_manifest_sha256=pan["artifact_manifest_sha256"],
    )


def _validate_bound_array(
    archive: Mapping[str, np.ndarray],
    bindings: Mapping[str, Any],
    name: str,
) -> np.ndarray:
    _require(name in archive and name in bindings, f"role array missing: {name}")
    value = np.asarray(archive[name])
    _require(array_binding(value) == bindings[name], f"role array drifted: {name}")
    return value


def load_role_bundle(name: str, record: Mapping[str, Any]) -> RoleBundle:
    """Authenticate and load only frozen prediction/operator arrays, never truth."""
    incumbent_manifest_path = _safe_repo_file(
        record["incumbent_manifest"], f"{name} incumbent manifest"
    )
    operator_manifest_path = _safe_repo_file(
        record["operator_manifest"], f"{name} operator manifest"
    )
    _require(
        sha256_file(incumbent_manifest_path) == record["incumbent_manifest_sha256"],
        f"{name} incumbent manifest hash drifted",
    )
    _require(
        sha256_file(operator_manifest_path) == record["operator_manifest_sha256"],
        f"{name} operator manifest hash drifted",
    )
    incumbent_manifest = load_json(incumbent_manifest_path)
    operator_manifest = load_json(operator_manifest_path)
    incumbent_path = _safe_repo_file(
        incumbent_manifest["npz"]["path"], f"{name} incumbent NPZ"
    )
    operator_path = _safe_repo_file(
        operator_manifest["npz"]["path"], f"{name} operator NPZ"
    )
    _require(
        sha256_file(incumbent_path) == record["incumbent_npz_sha256"],
        f"{name} incumbent NPZ hash drifted",
    )
    _require(
        sha256_file(operator_path) == record["operator_npz_sha256"],
        f"{name} operator NPZ hash drifted",
    )
    with np.load(incumbent_path, allow_pickle=False) as archive:
        bindings = incumbent_manifest["arrays"]
        targets = _validate_bound_array(archive, bindings, "targets")
        genes = _validate_bound_array(archive, bindings, "genes")
        score_indices = _validate_bound_array(archive, bindings, "score_indices")
        incumbent = _validate_bound_array(archive, bindings, "effects")
    with np.load(operator_path, allow_pickle=False) as archive:
        bindings = operator_manifest["arrays"]
        eigenvalues = _validate_bound_array(archive, bindings, "evals")
        eigenvectors = _validate_bound_array(archive, bindings, "evecs")
        columns = _validate_bound_array(archive, bindings, "columns")
    _require(len(targets) == record["target_count"], f"{name} target count changed")
    _require(
        sha256_array(targets) == record["target_array_sha256"],
        f"{name} target axis changed",
    )
    _require(
        len(score_indices) == record["score_gene_count"],
        f"{name} score-axis count changed",
    )
    _require(
        sha256_array(score_indices) == record["score_index_array_sha256"],
        f"{name} score indices changed",
    )
    _require(
        len(columns) == record["operator_column_count"],
        f"{name} operator-column count changed",
    )
    _require(
        sha256_array(columns) == record["operator_column_array_sha256"],
        f"{name} operator columns changed",
    )
    _require(
        incumbent.shape == (len(targets), len(genes)), f"{name} incumbent shape changed"
    )
    _require(
        eigenvectors.shape == (len(columns), len(columns)),
        f"{name} eigensystem shape changed",
    )
    _require(eigenvalues.shape == (len(columns),), f"{name} eigenvalue shape changed")
    return RoleBundle(
        name=name,
        targets=targets.astype(str),
        genes=genes.astype(str),
        score_indices=np.asarray(score_indices, dtype=np.int64),
        incumbent=np.asarray(incumbent, dtype=np.float64),
        eigenvalues=np.asarray(eigenvalues, dtype=np.float64),
        eigenvectors=np.asarray(eigenvectors, dtype=np.float64),
        operator_columns=np.asarray(columns, dtype=np.int64),
        input_bindings=(
            {
                "path": record["incumbent_manifest"],
                "sha256": record["incumbent_manifest_sha256"],
            },
            {
                "path": incumbent_manifest["npz"]["path"],
                "sha256": record["incumbent_npz_sha256"],
            },
            {
                "path": record["operator_manifest"],
                "sha256": record["operator_manifest_sha256"],
            },
            {
                "path": operator_manifest["npz"]["path"],
                "sha256": record["operator_npz_sha256"],
            },
        ),
    )


def apply_row_supported_whitener(
    raw_operator: np.ndarray,
    coordinate_supported: np.ndarray,
    eigenvalues: np.ndarray,
    eigenvectors: np.ndarray,
    *,
    gamma: float = GAMMA,
    lam: float = LAMBDA,
) -> np.ndarray:
    """Apply each row's supported principal block without fitting a new operator."""
    raw = np.asarray(raw_operator, dtype=np.float64)
    supported = np.asarray(coordinate_supported, dtype=bool)
    values = np.asarray(eigenvalues, dtype=np.float64)
    vectors = np.asarray(eigenvectors, dtype=np.float64)
    _require(
        raw.ndim == 2 and supported.shape == raw.shape, "source support does not align"
    )
    _require(
        vectors.shape == (len(values), len(values)),
        "operator eigensystem does not align",
    )
    _require(raw.shape[1] == len(values), "source and operator dimensions differ")
    _require(
        np.all(np.isfinite(values)) and np.all(values > 0),
        "operator eigenvalues are invalid",
    )
    _require(np.all(np.isfinite(raw)), "source effects contain nonfinite values")
    masked = np.where(supported, raw, 0.0)
    weights = (values + float(lam)) ** (-float(gamma) / 2.0)
    transformed = ((masked @ vectors) * weights) @ vectors.T
    transformed[~supported] = 0.0
    return transformed


def blend_row_supported(
    incumbent: np.ndarray,
    source: np.ndarray,
    coordinate_supported: np.ndarray,
    target_supported: np.ndarray,
    fraction: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Norm-match and blend with exact row/coordinate fallback."""
    if float(fraction) == 0.0:
        rows = np.asarray(target_supported, dtype=bool)
        return np.asarray(incumbent).copy(), np.zeros_like(rows)
    residual, effective, effective_support = normalized_residual_from_assignment(
        incumbent,
        source,
        coordinate_supported,
        target_supported,
    )
    return candidate_from_normalized_residual(
        incumbent,
        residual,
        effective_support,
        effective,
        fraction,
    ), effective


def normalized_residual_from_assignment(
    incumbent: np.ndarray,
    source: np.ndarray,
    coordinate_supported: np.ndarray,
    eligible_rows: np.ndarray,
    source_row_for_destination: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build scaled-source minus destination-incumbent residuals.

    A row is effective only when it has a valid assigned source row, at least
    one supported coordinate, finite values, positive source/incumbent norms,
    and a finite nonzero realized residual. Invalid rows are exact fallback.
    """
    base = np.asarray(incumbent)
    added = np.asarray(source, dtype=np.float64)
    support = np.asarray(coordinate_supported, dtype=bool)
    rows = np.asarray(eligible_rows, dtype=bool)
    _require(base.ndim == 2 and added.shape == base.shape, "blend arrays do not align")
    _require(support.shape == base.shape, "coordinate support does not align")
    _require(rows.shape == (base.shape[0],), "target support does not align")
    if source_row_for_destination is None:
        assignment = np.arange(len(rows), dtype=np.int64)
    else:
        assignment = np.asarray(source_row_for_destination, dtype=np.int64)
        _require(assignment.shape == rows.shape, "source-row assignment does not align")
        _require(
            np.all((assignment >= -1) & (assignment < len(rows))),
            "source-row assignment is out of bounds",
        )
    residual = np.zeros_like(added)
    effective = np.zeros_like(rows)
    effective_support = np.zeros_like(support)
    for destination in np.flatnonzero(rows):
        source_index = int(assignment[destination])
        if source_index < 0:
            continue
        mask = support[source_index]
        if not np.any(mask):
            continue
        source_row = added[source_index, mask]
        incumbent_row = np.asarray(base[destination, mask], dtype=np.float64)
        if not np.all(np.isfinite(source_row)) or not np.all(
            np.isfinite(incumbent_row)
        ):
            continue
        source_norm = float(np.linalg.norm(source_row))
        incumbent_norm = float(np.linalg.norm(incumbent_row))
        if (
            not math.isfinite(source_norm)
            or source_norm <= 0.0
            or not math.isfinite(incumbent_norm)
            or incumbent_norm <= 0.0
        ):
            continue
        scaled = source_row * (incumbent_norm / source_norm)
        row_residual = scaled - incumbent_row
        residual_norm = float(np.linalg.norm(row_residual))
        if (
            not np.all(np.isfinite(row_residual))
            or not math.isfinite(residual_norm)
            or residual_norm <= 0.0
        ):
            continue
        residual[destination, mask] = row_residual
        effective[destination] = True
        effective_support[destination, mask] = True
    return residual, effective, effective_support


def candidate_from_normalized_residual(
    incumbent: np.ndarray,
    normalized_residual: np.ndarray,
    effective_coordinate_supported: np.ndarray,
    effective_rows: np.ndarray,
    fraction: float,
) -> np.ndarray:
    """Apply a frozen normalized residual while preserving fallback bytes."""
    base = np.asarray(incumbent)
    residual = np.asarray(normalized_residual, dtype=np.float64)
    support = np.asarray(effective_coordinate_supported, dtype=bool)
    rows = np.asarray(effective_rows, dtype=bool)
    _require(base.ndim == 2 and residual.shape == base.shape, "blend arrays do not align")
    _require(support.shape == base.shape, "coordinate support does not align")
    _require(rows.shape == (base.shape[0],), "effective rows do not align")
    _require(
        not np.any(support[~rows]),
        "an ineffective row has effective supported coordinates",
    )
    _require(float(fraction) in (0.0, *FRACTIONS), "fraction is outside frozen grid")
    candidate = base.copy()
    if fraction == 0.0:
        return candidate
    candidate[support] = (
        np.asarray(base[support], dtype=np.float64)
        + float(fraction) * residual[support]
    ).astype(base.dtype, copy=False)
    fallback = (~rows)[:, None] | ~support
    _require(
        np.array_equal(
            np.ascontiguousarray(candidate[fallback]).view(np.uint8),
            np.ascontiguousarray(base[fallback]).view(np.uint8),
        ),
        "unsupported target or coordinate changed",
    )
    return candidate


def reconstruct_target_label_null(
    incumbent: np.ndarray,
    source_effects: np.ndarray,
    coordinate_supported: np.ndarray,
    target_novel: np.ndarray,
    destination_to_source: np.ndarray,
    *,
    population: str = "novel",
    fraction: float = PRIMARY_FRACTION,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Rebuild one label-null arm from the assigned source rows.

    The assigned source direction and its support move together. Norm matching
    is then recomputed against the destination incumbent; precomputed observed
    residuals are never permuted.
    """
    mapping = np.asarray(destination_to_source, dtype=np.int64)
    novel = np.asarray(target_novel, dtype=bool)
    _require(mapping.shape == novel.shape, "null assignment does not align")
    _require(population in {"novel", "all"}, "unknown null population")
    assigned = mapping >= 0
    if np.any(assigned):
        destinations = np.flatnonzero(assigned)
        _require(
            np.array_equal(novel[destinations], novel[mapping[destinations]]),
            "null assignment crosses novelty strata",
        )
    eligible = assigned & (novel if population == "novel" else True)
    residual, effective, effective_support = normalized_residual_from_assignment(
        incumbent,
        source_effects,
        coordinate_supported,
        eligible,
        mapping,
    )
    candidate = candidate_from_normalized_residual(
        incumbent,
        residual,
        effective_support,
        effective,
        fraction,
    )
    return candidate, residual, effective, effective_support


def deterministic_target_label_permutations(
    target_supported: np.ndarray,
    target_novel: np.ndarray,
    *,
    role: str,
    replicates: int = PERMUTATION_REPLICATES,
) -> np.ndarray:
    """Freeze destination-to-source row maps, stratified by novelty."""
    supported = np.asarray(target_supported, dtype=bool)
    novel = np.asarray(target_novel, dtype=bool)
    _require(supported.shape == novel.shape, "permutation strata do not align")
    result = np.full((replicates, len(supported)), -1, dtype=np.int32)
    for draw in range(replicates):
        seed = deterministic_seed(PERMUTATION_NAMESPACE, role, draw)
        rng = np.random.Generator(np.random.PCG64(seed))
        for novelty in (False, True):
            rows = np.flatnonzero(supported & (novel == novelty))
            if len(rows) == 0:
                continue
            if len(rows) == 1:
                permuted = rows
            else:
                offset = int(rng.integers(1, len(rows)))
                permuted = np.roll(rng.permutation(rows), offset)
                if np.any(permuted == rows):
                    permuted = np.roll(rows, 1 + draw % (len(rows) - 1))
            result[draw, rows] = permuted.astype(np.int32)
    return result


def reciprocal_knn_clusters(
    incumbent: np.ndarray,
    score_indices: np.ndarray,
    *,
    k: int = 5,
    minimum_cosine: float = 0.80,
) -> np.ndarray:
    values = np.asarray(incumbent[:, score_indices], dtype=np.float64)
    _require(values.ndim == 2 and len(values) >= 2, "cluster input is invalid")
    norms = np.linalg.norm(values, axis=1)
    unit = np.divide(
        values, norms[:, None], out=np.zeros_like(values), where=norms[:, None] > 0
    )
    similarity = unit @ unit.T
    np.fill_diagonal(similarity, -np.inf)
    count = min(max(1, int(k)), len(values) - 1)
    neighbors = np.argsort(-similarity, axis=1, kind="stable")[:, :count]
    adjacency = [set() for _ in range(len(values))]
    for left in range(len(values)):
        for right in neighbors[left]:
            right = int(right)
            if left in neighbors[right] and similarity[left, right] >= minimum_cosine:
                adjacency[left].add(right)
                adjacency[right].add(left)
    labels = np.full(len(values), -1, dtype=np.int32)
    next_label = 0
    for start in range(len(values)):
        if labels[start] >= 0:
            continue
        stack = [start]
        labels[start] = next_label
        while stack:
            node = stack.pop()
            for neighbor in sorted(adjacency[node]):
                if labels[neighbor] < 0:
                    labels[neighbor] = next_label
                    stack.append(neighbor)
        next_label += 1
    return labels


def directional_rank_indices(
    normalized_residual: np.ndarray,
    coordinate_supported: np.ndarray,
    effective_rows: np.ndarray,
) -> np.ndarray:
    effects = np.asarray(normalized_residual, dtype=np.float64)
    support = np.asarray(coordinate_supported, dtype=bool)
    effective = np.asarray(effective_rows, dtype=bool)
    _require(effects.shape == support.shape, "directional ranking arrays do not align")
    _require(effective.shape == (len(effects),), "directional rows do not align")
    active_support = support & effective[:, None]
    _require(
        np.all(np.isfinite(effects[active_support])),
        "directional residual contains nonfinite values",
    )
    maximum = int(np.max(np.sum(active_support, axis=1), initial=0))
    rankings = np.full((len(effects), maximum), -1, dtype=np.int32)
    for row in np.flatnonzero(effective):
        indices = np.flatnonzero(active_support[row])
        order = np.lexsort((indices, -np.abs(effects[row, indices])))
        rankings[row, : len(indices)] = indices[order].astype(np.int32)
    return rankings


def deterministic_cross_context_aggregate(
    surface_values: Mapping[str, float],
) -> dict[str, float]:
    """Aggregate six Jurkat galleries, then HEK293T/Jurkat with equal context weight."""
    _require(
        set(surface_values) == set(CROSS_CONTEXT_ROLES),
        "cross-context surfaces are incomplete or contain extras",
    )
    values = {name: float(surface_values[name]) for name in CROSS_CONTEXT_ROLES}
    _require(
        all(math.isfinite(value) for value in values.values()),
        "cross-context surface value is nonfinite",
    )
    jurkat = math.fsum(values[name] for name in NADIG_ROLES) / len(NADIG_ROLES)
    cross_context = math.fsum((values["orion_hek293t"], jurkat)) / 2.0
    return {
        "orion_hek293t": values["orion_hek293t"],
        "nadig_jurkat": jurkat,
        "cross_context": cross_context,
    }


def target_cluster_bootstrap_weights(
    cluster_labels: np.ndarray,
    *,
    role: str,
    replicate_zero_based: int,
) -> np.ndarray:
    """Return deterministic whole-cluster multiplicities for one role draw."""
    labels = np.asarray(cluster_labels, dtype=np.int64)
    _require(labels.ndim == 1 and labels.size > 0, "bootstrap labels are invalid")
    _require(np.all(labels >= 0), "bootstrap labels must be nonnegative")
    _require(
        0 <= replicate_zero_based < BOOTSTRAP_REPLICATES,
        "bootstrap replicate is outside the frozen range",
    )
    unique = np.unique(labels)
    seed = deterministic_seed(
        BOOTSTRAP_NAMESPACE,
        role,
        replicate_zero_based,
    )
    rng = np.random.Generator(np.random.PCG64(seed))
    sampled = rng.choice(unique, size=len(unique), replace=True)
    counts = np.bincount(sampled, minlength=int(unique[-1]) + 1)
    return counts[labels].astype(np.int64, copy=False)


def type7_quantile(values: np.ndarray, probability: float) -> float:
    """Freeze NumPy's linear quantile, equivalent to Hyndman-Fan type 7."""
    array = np.asarray(values, dtype=np.float64)
    _require(
        array.ndim == 1 and array.size > 0 and np.all(np.isfinite(array)),
        "quantile input is invalid",
    )
    _require(0.0 <= probability <= 1.0, "quantile probability is invalid")
    return float(np.quantile(array, probability, method="linear"))


def build_candidate_arrays(
    role: RoleBundle,
    compact: CompactDonor,
    panel_targets: Sequence[str],
    old_h1_targets: Sequence[str],
    *,
    enforce_production_coverage: bool = True,
    expected_pan_overlap: int | None = None,
    expected_novel_overlap: int | None = None,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Construct every frozen arm using source and incumbent predictions only."""
    _require(
        tuple(role.genes.astype(str)) == compact.genes, "role and Pan gene axes differ"
    )
    source_index = {target: index for index, target in enumerate(compact.targets)}
    panel_set = set(map(str, panel_targets))
    old_h1_set = set(map(str, old_h1_targets))
    held_target_set = set(role.targets.astype(str))
    excluded = np.asarray(
        [gene in panel_set or gene in held_target_set for gene in compact.genes],
        dtype=bool,
    )
    operator_mask = np.zeros(len(compact.genes), dtype=bool)
    operator_mask[role.operator_columns] = True
    _require(
        not np.any(operator_mask & excluded),
        f"{role.name}: operator includes excluded panel or held-target coordinates",
    )
    _require(
        not np.any(role.incumbent[:, excluded]),
        f"{role.name}: incumbent excluded coordinates are nonzero",
    )
    raw_operator = np.zeros(
        (len(role.targets), len(role.operator_columns)), dtype=np.float64
    )
    support_operator = np.zeros_like(raw_operator, dtype=bool)
    target_supported = np.zeros(len(role.targets), dtype=bool)
    target_novel = np.zeros(len(role.targets), dtype=bool)
    for row, target in enumerate(role.targets.astype(str)):
        source_row = source_index.get(target)
        if source_row is None:
            continue
        target_supported[row] = True
        target_novel[row] = target not in old_h1_set
        measured = np.asarray(compact.measured[source_row], dtype=bool)
        allowed = measured[role.operator_columns]
        if not np.any(allowed):
            continue
        support_operator[row] = allowed
        values = np.asarray(
            compact.delta[source_row, role.operator_columns], dtype=np.float64
        )
        if not np.all(np.isfinite(values[allowed])):
            continue
        raw_operator[row, allowed] = values[allowed]
    transformed_operator = apply_row_supported_whitener(
        raw_operator,
        support_operator,
        role.eigenvalues,
        role.eigenvectors,
    )
    source_effects = np.zeros_like(role.incumbent)
    coordinate_supported = np.zeros_like(role.incumbent, dtype=bool)
    source_effects[:, role.operator_columns] = transformed_operator
    coordinate_supported[:, role.operator_columns] = support_operator
    _require(
        not np.any(coordinate_supported[:, excluded]),
        "excluded coordinate became supported",
    )
    _require(
        not np.any(source_effects[:, excluded]), "source changed an excluded coordinate"
    )

    normalized_residual, effective_all, effective_support = (
        normalized_residual_from_assignment(
            role.incumbent,
            source_effects,
            coordinate_supported,
            target_supported,
        )
    )
    effective_novel = effective_all & target_novel
    novel_residual = np.where(effective_novel[:, None], normalized_residual, 0.0)
    novel_support = effective_support & effective_novel[:, None]

    arrays: dict[str, np.ndarray] = {
        "targets": role.targets.astype("U"),
        "genes": role.genes.astype("U"),
        "score_indices": role.score_indices.astype(np.int64),
        "incumbent_effects": role.incumbent.astype(np.float64),
        "source_effects": source_effects,
        "coordinate_supported": coordinate_supported,
        "target_supported": target_supported,
        "target_novel_vs_vcc2025_h1": target_novel,
        "normalized_residual": normalized_residual,
        "target_effective_all": effective_all,
        "target_effective_novel": effective_novel,
        "identity_f000": role.incumbent.copy(),
    }
    effective_counts: dict[str, int] = {}
    for population, residual, effective, support in (
        ("novel", novel_residual, effective_novel, novel_support),
        ("all", normalized_residual, effective_all, effective_support),
    ):
        for fraction in FRACTIONS:
            suffix = f"f{int(round(fraction * 1000)):03d}"
            name = f"pan_{population}_{suffix}"
            candidate = candidate_from_normalized_residual(
                role.incumbent,
                residual,
                support,
                effective,
                fraction,
            )
            arrays[name] = candidate
            effective_counts[name] = int(effective.sum())
    arrays["directional_rank_indices"] = directional_rank_indices(
        novel_residual,
        novel_support,
        effective_novel,
    )
    arrays["directional_rank_indices_all"] = directional_rank_indices(
        normalized_residual,
        effective_support,
        effective_all,
    )
    arrays["target_cluster_labels"] = reciprocal_knn_clusters(
        role.incumbent, role.score_indices
    )
    arrays["target_label_null_permutations"] = deterministic_target_label_permutations(
        effective_all,
        target_novel,
        role=role.name,
    )
    panel_supported = tuple(
        target for target in panel_targets if target in source_index
    )
    panel_new = tuple(target for target in panel_supported if target not in old_h1_set)
    if enforce_production_coverage:
        _require(
            isinstance(expected_pan_overlap, int)
            and isinstance(expected_novel_overlap, int),
            f"{role.name}: expected role overlap counts were not supplied",
        )
        _require(
            int(target_supported.sum()) == expected_pan_overlap,
            f"{role.name}: Pan target overlap changed",
        )
        _require(
            int((target_supported & target_novel).sum()) == expected_novel_overlap,
            f"{role.name}: novel Pan target overlap changed",
        )
        _require(len(panel_supported) == PAN_PANEL_COUNT, "Pan panel coverage changed")
        _require(
            len(panel_new) == PAN_PANEL_NEW_COUNT,
            "Pan novel panel coverage changed",
        )
    audit = {
        "role": role.name,
        "targets": len(role.targets),
        "genes": len(role.genes),
        "operator_columns": len(role.operator_columns),
        "pan_supported_targets": int(target_supported.sum()),
        "pan_novel_targets": int((target_supported & target_novel).sum()),
        "expected_pan_overlap": expected_pan_overlap,
        "expected_novel_overlap": expected_novel_overlap,
        "effective_all_targets": int(effective_all.sum()),
        "effective_novel_targets": int(effective_novel.sum()),
        "excluded_coordinates": int(excluded.sum()),
        "supported_coordinate_entries": int(coordinate_supported.sum()),
        "effective_counts": effective_counts,
        "panel_supported_count": len(panel_supported),
        "panel_supported_ordered_sha256": sha256_lines(panel_supported),
        "panel_new_count": len(panel_new),
        "panel_new_ordered_sha256": sha256_lines(panel_new),
        "held_truth_opened": False,
    }
    return arrays, audit


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def write_deterministic_npz(path: Path, arrays: Mapping[str, np.ndarray]) -> None:
    """Write a byte-deterministic compressed NPZ without giant byte copies."""
    with zipfile.ZipFile(
        path,
        mode="w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=6,
        allowZip64=True,
    ) as archive:
        for name in sorted(arrays):
            value = np.ascontiguousarray(np.asarray(arrays[name]))
            _require(not value.dtype.hasobject, f"object array forbidden: {name}")
            member = zipfile.ZipInfo(
                filename=f"{name}.npy",
                date_time=FIXED_ZIP_TIMESTAMP,
            )
            member.compress_type = zipfile.ZIP_DEFLATED
            member.create_system = 3
            member.external_attr = 0o600 << 16
            with archive.open(member, mode="w", force_zip64=True) as handle:
                np.lib.format.write_array(handle, value, allow_pickle=False)


def _write_role_output(
    output_root: Path,
    role: RoleBundle,
    arrays: Mapping[str, np.ndarray],
    audit: Mapping[str, Any],
    compact: CompactDonor,
) -> dict[str, Any]:
    directory = output_root / role.name.replace(".", "/")
    _require(not directory.exists(), f"refusing to overwrite role output: {directory}")
    directory.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(dir=directory.parent, prefix=f".{directory.name}.")
    )
    try:
        npz_path = temporary / "source-candidates.npz"
        write_deterministic_npz(npz_path, arrays)
        manifest = {
            "schema_version": 1,
            "experiment": EXPERIMENT,
            "status": "frozen_before_experiment_117_held_outcome_access",
            "role": role.name,
            "npz": {
                "path": "source-candidates.npz",
                "bytes": npz_path.stat().st_size,
                "sha256": sha256_file(npz_path),
            },
            "arrays": {name: array_binding(value) for name, value in arrays.items()},
            "audit": dict(audit),
            "inputs": [
                *role.input_bindings,
                {
                    "path": str(
                        (compact.directory / "compact-donor-receipt.json").relative_to(
                            ROOT
                        )
                    ),
                    "sha256": compact.receipt_sha256,
                },
                {
                    "path": str(
                        (compact.directory / "artifact-manifest.json").relative_to(ROOT)
                    ),
                    "sha256": compact.artifact_manifest_sha256,
                },
            ],
            "safety": {
                "held_truth_opened": False,
                "cell_eval2_run": False,
                "counts_generated": False,
                "submission_built": False,
            },
        }
        _atomic_json(temporary / "source-candidates.json", manifest)
        os.replace(temporary, directory)
        return manifest
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def build(output_root: Path = DEFAULT_OUTPUT_ROOT) -> dict[str, Any]:
    protocol, source, boundary = validate_scaffold()
    del protocol
    assert_real_build_authorized(source, boundary)
    output_root = _safe_output(output_root)
    _require(not output_root.exists(), f"refusing to overwrite {output_root}")
    compact = authenticate_compact_donor(source)
    panel_targets = _read_csv_axis(
        ROOT / source["canonical_axes"]["panel"]["path"],
        "target_gene",
        expected=PANEL_TARGET_COUNT,
    )
    old_h1_targets = _read_lines(
        ROOT / source["canonical_axes"]["vcc2025_h1_targets"]["path"],
        expected=source["canonical_axes"]["vcc2025_h1_targets"]["count"],
    )
    output_root.mkdir(parents=True)
    role_manifests: dict[str, Any] = {}
    try:
        for name, record in source["roles"].items():
            role = load_role_bundle(name, record)
            arrays, audit = build_candidate_arrays(
                role,
                compact,
                panel_targets,
                old_h1_targets,
                expected_pan_overlap=record["expected_pan_overlap"],
                expected_novel_overlap=record["expected_novel_overlap"],
            )
            role_manifests[name] = _write_role_output(
                output_root, role, arrays, audit, compact
            )
        top = {
            "schema_version": 1,
            "experiment": EXPERIMENT,
            "status": "source_candidates_frozen_before_held_outcome_access",
            "protocol_sha256": sha256_file(PROTOCOL_PATH),
            "source_manifest_sha256": sha256_file(SOURCE_MANIFEST_PATH),
            "freeze_boundary_sha256": sha256_file(FREEZE_BOUNDARY_PATH),
            "builder_sha256": sha256_file(Path(__file__).resolve()),
            "compact_receipt_sha256": compact.receipt_sha256,
            "roles": {
                role: {
                    "manifest": str(
                        (
                            output_root
                            / role.replace(".", "/")
                            / "source-candidates.json"
                        ).relative_to(ROOT)
                    ),
                    "sha256": sha256_file(
                        output_root / role.replace(".", "/") / "source-candidates.json"
                    ),
                }
                for role in role_manifests
            },
            "held_truth_opened": False,
        }
        _atomic_json(output_root / "source-candidate-inputs.json", top)
        return top
    except BaseException:
        shutil.rmtree(output_root, ignore_errors=True)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("verify-scaffold", "build"))
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "verify-scaffold":
        protocol, source, boundary = validate_scaffold()
        print(
            json.dumps(
                {
                    "experiment": protocol["experiment"],
                    "status": boundary["status"],
                    "roles": list(source["roles"]),
                    "real_build_authorized": boundary[
                        "source_candidate_build_authorized"
                    ],
                    "held_truth_opened": False,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return
    print(json.dumps(build(args.output_root), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
