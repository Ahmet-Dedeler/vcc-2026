#!/usr/bin/env python3
"""Fail-closed provenance audit and selective ZIP-range acquisition for exp 111."""

from __future__ import annotations

import argparse
import binascii
import csv
import fcntl
import hashlib
import io
import itertools
import json
import math
import os
import re
import shutil
import stat
import struct
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import zlib
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path, PurePosixPath
from typing import Any, BinaryIO, Iterator


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_MANIFEST = HERE / "SOURCE_MANIFEST.json"
DEFAULT_PROTOCOL = HERE / "preregistered_protocol.json"
DEFAULT_AUTHORIZATION = HERE / "EXECUTION_AUTHORIZATION.v1.json"
DEFAULT_SCHEMA_RECOVERY = HERE / "ACQUISITION_SCHEMA_RECOVERY.v1.json"
DEFAULT_EXTRACTED = REPO_ROOT / "data/external/depmap26q1/extracted"
DEFAULT_AUDIT = HERE / "audit"
RECEIPT_NAME = "acquisition_receipt.json"
RECEIPT_SCHEMA_VERSION = 1

ALLOWED_MEMBERS = (
    "metadata/cell_line_metadata.csv",
    "processed/pseudobulk_sum_by_ko.csv",
    "processed/lfc_matrix.csv",
    "processed/fdr_matrix.csv",
    "processed/crispr_table.csv",
    "processed/cell_quality_covariates.csv",
)
PROCESSED_MEMBERS = ALLOWED_MEMBERS[1:]
FORBIDDEN_ARCHIVES = (
    "data.zip",
    "downloaded.zip",
    "h2h_data_formatted.zip",
    "single_cell_data.zip",
)
AUTHORIZED_ACTIONS = {
    "selective_member_acquisition": True,
    "metadata_extraction": True,
    "source_schema_validation": True,
    "source_only_model_fitting": True,
    "complete_selection_nulls": True,
}
FORBIDDEN_ACTIONS = {
    "full_processed_zip_download": False,
    "single_cell_archive_download": False,
    "unlisted_member_access": False,
    "held_truth_access": False,
    "cell_eval2": False,
    "count_generation": False,
    "production_change": False,
    "challenge_build": False,
    "submission": False,
    "commit": False,
    "push": False,
}
EXPECTED_DATA_FILES = {
    67643355: ("data.zip", False),
    67643532: ("h2h_data_formatted.zip", False),
    67643358: ("metadata.zip", True),
    67643565: ("processed.zip", "range_only_selected_members"),
    67643568: ("single_cell_data.zip", False),
    67643802: ("downloaded.zip", False),
}
METADATA_COLUMNS = ("cell_line", "arxspan_id", "batch", "OncotreeLineage")
CRISPR_COLUMNS = (
    "cell_line",
    "guide",
    "gene_effect",
    "dependency_strength_class",
    "is_dependent",
)
QUALITY_COLUMNS = (
    "cell_line",
    "avg_umis_per_ko",
    "avg_genes_per_ko",
    "n_self_downregulated",
    "frac_self_downregulated",
    "total_umis_per_ko",
    "n_cells_per_ko",
)
HASH_RE = re.compile(r"^[0-9a-f]+$")
CONTENT_RANGE_RE = re.compile(r"^bytes ([0-9]+)-([0-9]+)/([0-9]+)$")
EOCD_SIGNATURE = b"PK\x05\x06"
CENTRAL_SIGNATURE = b"PK\x01\x02"
LOCAL_SIGNATURE = b"PK\x03\x04"
EOCD_STRUCT = struct.Struct("<4s4H2LH")
CENTRAL_STRUCT = struct.Struct("<4s6H3L5H2L")
LOCAL_STRUCT = struct.Struct("<4s5H3L2H")


class AcquisitionError(RuntimeError):
    """A source, archive, extraction, receipt, or schema contract failed."""


@dataclass(frozen=True)
class ZipEntry:
    """The central-directory fields needed for safe selective extraction."""

    name: str
    version_made_by: int
    version_needed: int
    flag_bits: int
    compression_method: int
    crc32: int
    compressed_size: int
    uncompressed_size: int
    disk_start: int
    internal_attributes: int
    external_attributes: int
    local_header_offset: int


@dataclass(frozen=True)
class LocalMember:
    """A validated local header plus its exact compressed byte interval."""

    entry: ZipEntry
    data_offset: int
    data_end: int


@dataclass
class RemoteAudit:
    """Small remote artifacts and verified archive metadata used by acquisition."""

    report: dict[str, Any]
    metadata_member: bytes
    processed_source: HttpsRangeSource
    processed_members: dict[str, LocalMember]


def _is_lower_hash(value: Any, length: int) -> bool:
    return (
        isinstance(value, str)
        and len(value) == length
        and HASH_RE.fullmatch(value) is not None
    )


def _require_hash(value: Any, length: int, label: str) -> str:
    if not _is_lower_hash(value, length):
        raise AcquisitionError(f"{label} must be a lowercase {length}-character hash")
    return value


def _require_https(url: Any, label: str) -> str:
    if not isinstance(url, str) or urllib.parse.urlsplit(url).scheme != "https":
        raise AcquisitionError(f"{label} must use HTTPS")
    return url


def _require_exact_keys(
    value: Any,
    expected: set[str],
    label: str,
) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise AcquisitionError(f"{label} must be an object")
    actual = set(value)
    if actual != expected:
        raise AcquisitionError(
            f"{label} keys drifted: missing={sorted(expected - actual)}, "
            f"extra={sorted(actual - expected)}"
        )
    return value


def _require_positive_int(value: Any, label: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise AcquisitionError(f"{label} must be a positive integer")
    return value


def _require_nonnegative_int(value: Any, label: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise AcquisitionError(f"{label} must be a nonnegative integer")
    return value


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def crc32_file(path: Path, block_bytes: int = 8 << 20) -> str:
    checksum = 0
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            checksum = binascii.crc32(block, checksum)
    return f"{checksum & 0xFFFFFFFF:08x}"


def canonical_json_sha256(payload: Any) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return sha256_bytes(encoded)


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def require_under(path: Path, root: Path, label: str) -> Path:
    resolved_root = root.resolve()
    resolved = path.resolve()
    if not resolved.is_relative_to(resolved_root):
        raise AcquisitionError(f"{label} escapes {resolved_root}: {resolved}")
    return resolved


def require_under_experiment(path: Path) -> Path:
    return require_under(path, HERE, "path")


def require_exact_source_root(path: Path) -> Path:
    expected = (REPO_ROOT / "data/external/depmap26q1/extracted").resolve()
    resolved = path.resolve()
    if resolved != expected:
        raise AcquisitionError(
            f"source root must remain exactly {expected}, received {resolved}"
        )
    return resolved


def require_safe_extraction_root(path: Path) -> Path:
    resolved = path.resolve()
    authorized = (REPO_ROOT / "data/external/depmap26q1/extracted").resolve()
    if resolved == authorized or resolved.is_relative_to(HERE.resolve()):
        return resolved
    raise AcquisitionError(
        f"extraction root must remain under {HERE.resolve()} or equal {authorized}: "
        f"{resolved}"
    )


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AcquisitionError(f"cannot read {label} {path}: {error}") from error
    if not isinstance(value, dict):
        raise AcquisitionError(f"{label} must contain a JSON object")
    return value


def _validate_figshare_file(record: Any, label: str) -> dict[str, Any]:
    if not isinstance(record, dict):
        raise AcquisitionError(f"{label} must be an object")
    _require_positive_int(record.get("id"), f"{label}.id")
    name = record.get("name")
    if not isinstance(name, str) or Path(name).name != name:
        raise AcquisitionError(f"{label}.name must be a plain basename")
    _require_positive_int(record.get("size_bytes"), f"{label}.size_bytes")
    _require_hash(record.get("md5"), 32, f"{label}.md5")
    _require_https(record.get("download_url"), f"{label}.download_url")
    sha256 = record.get("sha256")
    if sha256 is not None:
        _require_hash(sha256, 64, f"{label}.sha256")
    return record


def _validate_central_spec(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise AcquisitionError(f"{label} must be an object")
    _require_positive_int(value.get("entry_count"), f"{label}.entry_count")
    _require_nonnegative_int(value.get("offset"), f"{label}.offset")
    _require_positive_int(value.get("size_bytes"), f"{label}.size_bytes")
    _require_hash(value.get("sha256"), 64, f"{label}.sha256")
    return value


def _validate_member_spec(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise AcquisitionError(f"{label} must be an object")
    safe_zip_member(str(value.get("name", "")))
    _require_nonnegative_int(
        value.get("local_header_offset"),
        f"{label}.local_header_offset",
    )
    _require_positive_int(value.get("compressed_size"), f"{label}.compressed_size")
    _require_positive_int(
        value.get("uncompressed_size"),
        f"{label}.uncompressed_size",
    )
    _require_nonnegative_int(value.get("flag_bits"), f"{label}.flag_bits")
    method = value.get("compression_method")
    if method not in {0, 8}:
        raise AcquisitionError(f"{label}.compression_method must be stored or deflate")
    _require_hash(value.get("crc32"), 8, f"{label}.crc32")
    for hash_name, length in (("content_md5", 32), ("content_sha256", 64)):
        if hash_name in value:
            _require_hash(value.get(hash_name), length, f"{label}.{hash_name}")
    return value


def validate_manifest(manifest: dict[str, Any]) -> dict[str, Any]:
    if manifest.get("schema_version") != 1:
        raise AcquisitionError("unsupported source-manifest schema")
    if manifest.get("dataset") != "broad-depmap-perturb-seq-2026":
        raise AcquisitionError("unexpected source-manifest dataset")
    if manifest.get("preprint", {}).get("doi") != "10.64898/2026.08.24.746802":
        raise AcquisitionError("preprint DOI drifted")

    figshare = manifest.get("figshare")
    if not isinstance(figshare, dict):
        raise AcquisitionError("manifest has no Figshare source block")
    data_article = figshare.get("data_article")
    code_article = figshare.get("code_article")
    if not isinstance(data_article, dict) or data_article.get("id") != 33273600:
        raise AcquisitionError("Figshare data article must be 33273600")
    if not isinstance(code_article, dict) or code_article.get("id") != 33282171:
        raise AcquisitionError("Figshare code article must be 33282171")
    _require_https(data_article.get("api_url"), "data article API URL")
    _require_https(code_article.get("api_url"), "code article API URL")

    data_files = data_article.get("files")
    if not isinstance(data_files, list) or len(data_files) != len(EXPECTED_DATA_FILES):
        raise AcquisitionError("Figshare data-file inventory drifted")
    seen_ids: set[int] = set()
    for index, raw_file in enumerate(data_files):
        file_record = _validate_figshare_file(raw_file, f"data file {index}")
        file_id = file_record["id"]
        if file_id in seen_ids:
            raise AcquisitionError(f"duplicate Figshare data-file id {file_id}")
        seen_ids.add(file_id)
        if file_id not in EXPECTED_DATA_FILES:
            raise AcquisitionError(f"unexpected Figshare data-file id {file_id}")
        expected_name, expected_permission = EXPECTED_DATA_FILES[file_id]
        if (
            file_record["name"] != expected_name
            or file_record.get("permitted") != expected_permission
        ):
            raise AcquisitionError(f"permission/name drift for {expected_name}")
    if seen_ids != set(EXPECTED_DATA_FILES):
        raise AcquisitionError("Figshare data-file ids are incomplete")
    metadata_asset = next(item for item in data_files if item["name"] == "metadata.zip")
    processed_asset = next(
        item for item in data_files if item["name"] == "processed.zip"
    )
    if metadata_asset.get("sha256") is None:
        raise AcquisitionError("metadata.zip must have a pinned SHA-256")
    if (
        processed_asset.get("sha256") is not None
        or processed_asset.get("sha256_status")
        != "not_computable_without_forbidden_full_download"
    ):
        raise AcquisitionError("processed.zip full-archive SHA-256 status drifted")

    code_files = code_article.get("files")
    if not isinstance(code_files, list) or len(code_files) != 1:
        raise AcquisitionError("Figshare code article must contain exactly one file")
    code_file = _validate_figshare_file(code_files[0], "Figshare code archive")
    if code_file["id"] != 67644690 or code_file["name"] != "perturb-seq-depmap.zip":
        raise AcquisitionError("Figshare code archive identity drifted")
    _require_hash(
        code_file.get("normalized_tree_sha256"),
        64,
        "Figshare normalized tree SHA-256",
    )

    github = manifest.get("github")
    if not isinstance(github, dict):
        raise AcquisitionError("manifest has no GitHub source block")
    _require_hash(github.get("commit"), 40, "GitHub commit")
    github_archive = github.get("commit_archive")
    if not isinstance(github_archive, dict):
        raise AcquisitionError("manifest has no GitHub commit archive")
    _require_https(github_archive.get("url"), "GitHub archive URL")
    _require_positive_int(github_archive.get("size_bytes"), "GitHub archive size")
    _require_hash(github_archive.get("md5"), 32, "GitHub archive MD5")
    _require_hash(github_archive.get("sha256"), 64, "GitHub archive SHA-256")
    _require_hash(
        github_archive.get("normalized_tree_sha256"),
        64,
        "GitHub normalized tree SHA-256",
    )
    if github_archive["normalized_tree_sha256"] != code_file["normalized_tree_sha256"]:
        raise AcquisitionError("pinned GitHub/Figshare normalized trees disagree")
    if github.get("normalized_tree_algorithm") != (
        "sha256 over sorted safe relative files as path_utf8 NUL bytes NUL; "
        "strip one root; ignore __MACOSX, .DS_Store, and AppleDouble"
    ):
        raise AcquisitionError("normalized-tree algorithm is missing or drifted")

    archives = manifest.get("archives")
    if not isinstance(archives, dict) or set(archives) != {
        "metadata_zip",
        "processed_zip",
    }:
        raise AcquisitionError("archive contract is incomplete or unexpected")
    metadata_zip = archives["metadata_zip"]
    processed_zip = archives["processed_zip"]
    if metadata_zip.get("asset_file_id") != 67643358:
        raise AcquisitionError("metadata archive file id drifted")
    if processed_zip.get("asset_file_id") != 67643565:
        raise AcquisitionError("processed archive file id drifted")
    _validate_central_spec(metadata_zip.get("central_directory"), "metadata central")
    processed_central = _validate_central_spec(
        processed_zip.get("central_directory"),
        "processed central",
    )
    if processed_central.get("tail_probe_bytes") != 131072:
        raise AcquisitionError("processed ZIP tail probe must remain 131072 bytes")

    metadata_member = _validate_member_spec(
        metadata_zip.get("selected_member"),
        "metadata selected member",
    )
    if metadata_member["name"] != ALLOWED_MEMBERS[0]:
        raise AcquisitionError("metadata selected member drifted")
    selected = processed_zip.get("selected_members")
    if not isinstance(selected, list):
        raise AcquisitionError("processed selected_members must be a list")
    selected_by_name: dict[str, dict[str, Any]] = {}
    for index, raw_member in enumerate(selected):
        member = _validate_member_spec(raw_member, f"processed member {index}")
        name = member["name"]
        if name in selected_by_name:
            raise AcquisitionError(f"duplicate selected processed member {name}")
        selected_by_name[name] = member
    if tuple(selected_by_name) != PROCESSED_MEMBERS:
        raise AcquisitionError("processed selective-member allowlist/order drifted")

    contract = manifest.get("dataset_contract")
    if not isinstance(contract, dict):
        raise AcquisitionError("manifest has no dataset contract")
    lines = contract.get("cell_lines")
    if not isinstance(lines, list) or len(lines) != 16:
        raise AcquisitionError("dataset contract must pin exactly 16 cell lines")
    if contract.get("expected_cell_line_count") != 16:
        raise AcquisitionError("expected_cell_line_count must be 16")
    line_names: list[str] = []
    arxspan_ids: list[str] = []
    batches: set[str] = set()
    lineages: set[str] = set()
    for index, line in enumerate(lines):
        required = {"cell_line", "arxspan_id", "batch", "lineage"}
        _require_exact_keys(line, required, f"cell line {index}")
        if any(not isinstance(line[key], str) or not line[key] for key in required):
            raise AcquisitionError(f"cell line {index} contains an empty field")
        line_names.append(line["cell_line"])
        arxspan_ids.append(line["arxspan_id"])
        batches.add(line["batch"])
        lineages.add(line["lineage"])
    if len(set(line_names)) != 16 or len(set(arxspan_ids)) != 16:
        raise AcquisitionError("cell-line names and IDs must be unique")
    if len(batches) != 6 or contract.get("expected_batch_count") != 6:
        raise AcquisitionError("dataset contract must pin exactly six batches")
    if len(lineages) != 11 or contract.get("expected_lineage_count") != 11:
        raise AcquisitionError("dataset contract must pin exactly 11 lineages")

    targets = contract.get("targets")
    controls = contract.get("control_targets")
    if (
        not isinstance(targets, list)
        or len(targets) != 100
        or len(set(targets)) != 100
        or not all(isinstance(target, str) and target for target in targets)
    ):
        raise AcquisitionError("dataset contract must pin 100 unique targets")
    if (
        not isinstance(controls, list)
        or len(controls) != 10
        or len(set(controls)) != 10
        or not set(controls).issubset(targets)
    ):
        raise AcquisitionError("dataset contract must pin ten target controls")
    if contract.get("expected_target_count") != 100:
        raise AcquisitionError("expected_target_count must be 100")
    if contract.get("expected_control_target_count") != 10:
        raise AcquisitionError("expected_control_target_count must be 10")
    if contract.get("expected_analysis_target_count") != 90:
        raise AcquisitionError("expected_analysis_target_count must be 90")
    if contract.get("matrix_missing_values") != {
        "encoding": "empty_csv_field",
        "require_lfc_fdr_mask_identity": True,
    }:
        raise AcquisitionError("matrix missing-value contract drifted")

    pseudobulk = contract.get("pseudobulk")
    missing_pairs = pseudobulk.get("missing_context_target_pairs", [])
    normalized_missing = [tuple(pair) for pair in missing_pairs]
    if len(normalized_missing) != 7 or len(set(normalized_missing)) != 7:
        raise AcquisitionError("pseudobulk must pin seven unique missing pairs")
    valid_pairs = set((line, target) for line in line_names for target in targets)
    if not set(normalized_missing).issubset(valid_pairs):
        raise AcquisitionError("pseudobulk missing pair is outside the source grid")
    expected_pseudobulk_columns = 16 * 100 - 7
    if pseudobulk.get("expected_column_count") != expected_pseudobulk_columns:
        raise AcquisitionError("pseudobulk expected column count drifted")
    if contract.get("sceptre_matrices", {}).get("expected_column_count") != 16 * 90:
        raise AcquisitionError("SCEPTRE expected column count drifted")

    challenge = manifest.get("challenge_contract")
    if not isinstance(challenge, dict):
        raise AcquisitionError("manifest has no challenge-panel contract")
    if (
        challenge.get("expected_panel_overlap") != 0
        or challenge.get("expected_panel_target_count") != 300
    ):
        raise AcquisitionError("challenge-panel count/overlap contract drifted")
    panel_path = challenge.get("panel_path")
    if panel_path != "data/challenge/pert_counts.csv":
        raise AcquisitionError("challenge panel path drifted")
    _require_hash(challenge.get("panel_sha256"), 64, "challenge panel SHA-256")

    safety = manifest.get("safety")
    expected_safety = {
        "download_root": "downloads",
        "extraction_root": "extracted",
        "forbidden_file_ids": [67643355, 67643532, 67643568, 67643802],
        "max_full_download_bytes": 10000000,
        "max_range_chunk_bytes": 8388608,
        "max_selective_compressed_bytes": 350000000,
        "never_download_full_processed_zip": True,
        "never_download_single_cell_data": True,
        "only_https": True,
        "selected_processed_member_count": 5,
    }
    if safety != expected_safety:
        raise AcquisitionError("source-manifest safety contract drifted")
    compressed_total = sum(item["compressed_size"] for item in selected)
    if compressed_total > safety["max_selective_compressed_bytes"]:
        raise AcquisitionError("selected compressed members exceed the safety ceiling")
    return manifest


def validate_protocol(protocol: dict[str, Any]) -> dict[str, Any]:
    if protocol.get("schema_version") != 1:
        raise AcquisitionError("unsupported preregistration schema")
    if protocol.get("experiment") != "111-depmap-context-program-gain":
        raise AcquisitionError("unexpected preregistered experiment")
    source = protocol.get("source")
    if not isinstance(source, dict) or source.get("manifest") != "SOURCE_MANIFEST.json":
        raise AcquisitionError("preregistration source-manifest pointer drifted")
    if tuple(source.get("allowed_members", [])) != ALLOWED_MEMBERS:
        raise AcquisitionError("preregistered selective-member allowlist drifted")
    if set(source.get("forbidden_archives", [])) != set(FORBIDDEN_ARCHIVES):
        raise AcquisitionError("preregistered forbidden-archive set drifted")
    if source.get("forbidden_primary_representation") != "zscore_matrix.csv":
        raise AcquisitionError("z-score primary-use prohibition drifted")
    authorization = protocol.get("current_authorization")
    if not isinstance(authorization, dict):
        raise AcquisitionError("preregistration has no authorization block")
    for key in (
        "model_fitting",
        "held_truth_access",
        "count_generation",
        "production_change",
        "challenge_build",
        "submission",
        "commit",
        "push",
    ):
        if authorization.get(key) is not False:
            raise AcquisitionError(f"preregistration must keep {key}=false")
    return protocol


def validate_execution_authorization(
    authorization: dict[str, Any],
    *,
    manifest_path: Path = DEFAULT_MANIFEST,
    protocol_path: Path = DEFAULT_PROTOCOL,
) -> dict[str, Any]:
    expected_keys = {
        "schema_version",
        "experiment",
        "date_authorized",
        "status",
        "authorization_basis",
        "bindings",
        "allowed_members",
        "authorized_actions",
        "forbidden_actions",
        "output_policy",
        "fixed_source_model_contract",
        "recovery_amendments",
    }
    _require_exact_keys(authorization, expected_keys, "execution authorization")
    if authorization.get("schema_version") != 1:
        raise AcquisitionError("unsupported execution-authorization schema")
    if authorization.get("experiment") != "111-depmap-context-program-gain":
        raise AcquisitionError("execution-authorization experiment drifted")
    if authorization.get("status") != "authorized_source_only_execution":
        raise AcquisitionError("execution authorization is not active")
    if not isinstance(authorization.get("authorization_basis"), str):
        raise AcquisitionError("execution authorization basis is missing")
    if tuple(authorization.get("allowed_members", [])) != ALLOWED_MEMBERS:
        raise AcquisitionError("execution-authorization allowlist drifted")
    if authorization.get("authorized_actions") != AUTHORIZED_ACTIONS:
        raise AcquisitionError("execution-authorization allowed actions drifted")
    if authorization.get("forbidden_actions") != FORBIDDEN_ACTIONS:
        raise AcquisitionError("execution-authorization forbidden actions drifted")

    bindings = _require_exact_keys(
        authorization.get("bindings"),
        {"source_manifest", "preregistered_protocol"},
        "execution authorization bindings",
    )
    expected_bindings = {
        "source_manifest": {
            "path": str(DEFAULT_MANIFEST.relative_to(REPO_ROOT)),
            "sha256": sha256_file(manifest_path),
        },
        "preregistered_protocol": {
            "path": str(DEFAULT_PROTOCOL.relative_to(REPO_ROOT)),
            "sha256": sha256_file(protocol_path),
        },
    }
    if bindings != expected_bindings:
        raise AcquisitionError("execution-authorization source bindings drifted")

    recovery_bindings = authorization.get("recovery_amendments")
    expected_recovery_bindings = [
        {
            "path": str(DEFAULT_SCHEMA_RECOVERY.relative_to(REPO_ROOT)),
            "sha256": sha256_file(DEFAULT_SCHEMA_RECOVERY),
        }
    ]
    if recovery_bindings != expected_recovery_bindings:
        raise AcquisitionError("execution-authorization recovery binding drifted")
    validate_schema_recovery(
        load_json(DEFAULT_SCHEMA_RECOVERY, "acquisition schema recovery"),
        manifest_path=manifest_path,
        protocol_path=protocol_path,
    )

    output_policy = _require_exact_keys(
        authorization.get("output_policy"),
        {
            "source_root",
            "source_result_root",
            "refuse_source_overwrite",
            "refuse_result_overwrite",
        },
        "execution authorization output policy",
    )
    if output_policy != {
        "source_root": "data/external/depmap26q1/extracted",
        "source_result_root": (
            "experiments/111-depmap-context-program-gain/work/source-only-v1"
        ),
        "refuse_source_overwrite": True,
        "refuse_result_overwrite": True,
    }:
        raise AcquisitionError("execution-authorization output policy drifted")

    fixed_contract = _require_exact_keys(
        authorization.get("fixed_source_model_contract"),
        {
            "source_model_sha256",
            "source_validation_sha256",
            "source_promotion_gates_sha256",
            "held_truth_stops_before_access",
        },
        "fixed source-model contract",
    )
    protocol = load_json(protocol_path, "preregistration")
    expected_contract = {
        "source_model_sha256": canonical_json_sha256(protocol["source_model"]),
        "source_validation_sha256": canonical_json_sha256(
            protocol["source_validation"]
        ),
        "source_promotion_gates_sha256": canonical_json_sha256(
            protocol["source_promotion_gates"]
        ),
        "held_truth_stops_before_access": True,
    }
    if fixed_contract != expected_contract:
        raise AcquisitionError("fixed source-model contract drifted")
    return authorization


def validate_schema_recovery(
    recovery: dict[str, Any],
    *,
    manifest_path: Path = DEFAULT_MANIFEST,
    protocol_path: Path = DEFAULT_PROTOCOL,
) -> dict[str, Any]:
    expected = {
        "schema_version": 1,
        "experiment": "111-depmap-context-program-gain",
        "date_recorded": "2026-08-30",
        "status": "authorized_exact_source_schema_recovery",
        "bindings": {
            "source_manifest_sha256": sha256_file(manifest_path),
            "preregistered_protocol_sha256": sha256_file(protocol_path),
            "failed_execution_authorization_sha256": (
                "33fa3e55f1df461fd9e1ae9ffb42e27b75ded36141cc15b67897cfb7e5e456b4"
            ),
        },
        "crispr_table_exception": {
            "member_size_bytes": 76367,
            "member_crc32": "ff7d6f0c",
            "member_sha256": (
                "cf10cfedf2036478298a0b595a4cdd82dc6682d07aa42caa553618f59b26ea41"
            ),
            "row_count": 1600,
            "missing_row_count": 70,
            "sorted_cell_line_guide_pairs_sha256": (
                "cf7d62365d6898ad27856a6773b3ca2c2533703b677a6ee75cc0c9012feb480a"
            ),
            "missing_fields_exactly": [
                "gene_effect",
                "dependency_strength_class",
            ],
            "missing_rows_must_be_control_targets": True,
            "is_dependent_must_equal": "False",
            "all_other_fields_nonempty": True,
        },
        "scope": {
            "only_crispr_table_control_metadata_missingness": True,
            "matrix_empty_field_contract_unchanged": True,
            "analysis_target_rows_must_remain_complete": True,
            "source_model_grid_and_gates_unchanged": True,
            "held_truth_access": False,
            "cell_eval2": False,
            "counts": False,
            "production_change": False,
            "challenge_build": False,
            "submission": False,
            "commit": False,
            "push": False,
        },
    }
    if recovery != expected:
        raise AcquisitionError("acquisition schema-recovery contract drifted")
    return recovery


def load_contracts(
    manifest_path: Path = DEFAULT_MANIFEST,
    protocol_path: Path = DEFAULT_PROTOCOL,
) -> tuple[dict[str, Any], dict[str, Any]]:
    manifest_path = require_under_experiment(manifest_path)
    protocol_path = require_under_experiment(protocol_path)
    manifest = validate_manifest(load_json(manifest_path, "source manifest"))
    protocol = validate_protocol(load_json(protocol_path, "preregistration"))
    validate_challenge_panel(manifest)
    return manifest, protocol


def load_execution_contracts(
    manifest_path: Path = DEFAULT_MANIFEST,
    protocol_path: Path = DEFAULT_PROTOCOL,
    authorization_path: Path = DEFAULT_AUTHORIZATION,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    manifest, protocol = load_contracts(manifest_path, protocol_path)
    authorization_path = require_under_experiment(authorization_path)
    authorization = validate_execution_authorization(
        load_json(authorization_path, "execution authorization"),
        manifest_path=manifest_path,
        protocol_path=protocol_path,
    )
    return manifest, protocol, authorization


def validate_challenge_panel(manifest: dict[str, Any]) -> dict[str, Any]:
    challenge = manifest["challenge_contract"]
    panel_path = require_under(
        REPO_ROOT / challenge["panel_path"],
        REPO_ROOT,
        "challenge panel path",
    )
    if panel_path.is_symlink() or not panel_path.is_file():
        raise AcquisitionError(f"challenge panel is missing or symlinked: {panel_path}")
    actual_hash = sha256_file(panel_path)
    if actual_hash != challenge["panel_sha256"]:
        raise AcquisitionError(
            f"challenge panel SHA-256 {actual_hash} != pinned "
            f"{challenge['panel_sha256']}"
        )
    with panel_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["target_gene"]:
            raise AcquisitionError("challenge panel header must be target_gene")
        panel = [row["target_gene"] for row in reader]
    if len(panel) != 300 or len(set(panel)) != 300 or any(not item for item in panel):
        raise AcquisitionError("challenge panel must contain 300 unique targets")
    overlap = sorted(set(panel) & set(manifest["dataset_contract"]["targets"]))
    if len(overlap) != challenge["expected_panel_overlap"]:
        raise AcquisitionError(f"source/VCC target overlap drifted: {overlap}")
    return {
        "path": str(panel_path.relative_to(REPO_ROOT)),
        "sha256": actual_hash,
        "target_count": len(panel),
        "source_overlap": overlap,
    }


def safe_zip_member(name: str) -> PurePosixPath:
    if not name or "\x00" in name or "\\" in name or name.startswith("/"):
        raise AcquisitionError(f"unsafe ZIP member path: {name!r}")
    trimmed = name[:-1] if name.endswith("/") else name
    if not trimmed or "//" in trimmed:
        raise AcquisitionError(f"unsafe ZIP member path: {name!r}")
    path = PurePosixPath(trimmed)
    if any(part in {"", ".", ".."} for part in path.parts):
        raise AcquisitionError(f"unsafe ZIP member path: {name!r}")
    if path.parts and re.match(r"^[A-Za-z]:$", path.parts[0]):
        raise AcquisitionError(f"unsafe ZIP member drive path: {name!r}")
    return path


def _decode_zip_name(raw_name: bytes, flag_bits: int) -> str:
    encoding = "utf-8" if flag_bits & 0x800 else "cp437"
    try:
        return raw_name.decode(encoding)
    except UnicodeDecodeError as error:
        raise AcquisitionError("ZIP member name is not decodable") from error


def parse_eocd(tail: bytes, tail_start: int, archive_size: int) -> dict[str, int]:
    position = tail.rfind(EOCD_SIGNATURE)
    if position < 0 or position + EOCD_STRUCT.size > len(tail):
        raise AcquisitionError("ZIP EOCD was not found in the pinned tail probe")
    (
        signature,
        disk_number,
        central_disk,
        entries_on_disk,
        entry_count,
        central_size,
        central_offset,
        comment_length,
    ) = EOCD_STRUCT.unpack_from(tail, position)
    if signature != EOCD_SIGNATURE:
        raise AcquisitionError("invalid ZIP EOCD signature")
    if position + EOCD_STRUCT.size + comment_length != len(tail):
        raise AcquisitionError("ZIP EOCD/comment does not end at archive EOF")
    if disk_number != 0 or central_disk != 0 or entries_on_disk != entry_count:
        raise AcquisitionError("multi-disk ZIP archives are forbidden")
    if (
        entry_count == 0xFFFF
        or central_size == 0xFFFFFFFF
        or central_offset == 0xFFFFFFFF
    ):
        raise AcquisitionError("ZIP64 central directories are not authorized")
    absolute_eocd = tail_start + position
    if central_offset + central_size != absolute_eocd:
        raise AcquisitionError("ZIP central directory is not contiguous with EOCD")
    if absolute_eocd + EOCD_STRUCT.size + comment_length != archive_size:
        raise AcquisitionError("ZIP EOCD absolute position does not match archive size")
    return {
        "entry_count": entry_count,
        "offset": central_offset,
        "size_bytes": central_size,
        "comment_length": comment_length,
    }


def _entry_is_symlink(entry: ZipEntry) -> bool:
    creator_system = entry.version_made_by >> 8
    unix_mode = (entry.external_attributes >> 16) & 0xFFFF
    return creator_system == 3 and stat.S_ISLNK(unix_mode)


def parse_central_directory(
    payload: bytes,
    expected_entry_count: int,
) -> dict[str, ZipEntry]:
    entries: dict[str, ZipEntry] = {}
    cursor = 0
    for index in range(expected_entry_count):
        if cursor + CENTRAL_STRUCT.size > len(payload):
            raise AcquisitionError(f"central directory truncated at entry {index}")
        fields = CENTRAL_STRUCT.unpack_from(payload, cursor)
        (
            signature,
            version_made_by,
            version_needed,
            flag_bits,
            compression_method,
            _mod_time,
            _mod_date,
            crc32,
            compressed_size,
            uncompressed_size,
            name_length,
            extra_length,
            comment_length,
            disk_start,
            internal_attributes,
            external_attributes,
            local_header_offset,
        ) = fields
        if signature != CENTRAL_SIGNATURE:
            raise AcquisitionError(
                f"invalid central-directory signature at entry {index}"
            )
        end = cursor + CENTRAL_STRUCT.size + name_length + extra_length + comment_length
        if end > len(payload):
            raise AcquisitionError(f"central-directory entry {index} is truncated")
        raw_name = payload[
            cursor + CENTRAL_STRUCT.size : cursor + CENTRAL_STRUCT.size + name_length
        ]
        name = _decode_zip_name(raw_name, flag_bits)
        safe_zip_member(name)
        if name in entries:
            raise AcquisitionError(f"duplicate ZIP member path: {name}")
        if flag_bits & 0x1:
            raise AcquisitionError(f"encrypted ZIP member is forbidden: {name}")
        if flag_bits & 0x40:
            raise AcquisitionError(f"strong-encryption ZIP flag is forbidden: {name}")
        if disk_start != 0:
            raise AcquisitionError(f"multi-disk ZIP member is forbidden: {name}")
        if (
            compressed_size == 0xFFFFFFFF
            or uncompressed_size == 0xFFFFFFFF
            or local_header_offset == 0xFFFFFFFF
        ):
            raise AcquisitionError(f"ZIP64 member is not authorized: {name}")
        entry = ZipEntry(
            name=name,
            version_made_by=version_made_by,
            version_needed=version_needed,
            flag_bits=flag_bits,
            compression_method=compression_method,
            crc32=crc32,
            compressed_size=compressed_size,
            uncompressed_size=uncompressed_size,
            disk_start=disk_start,
            internal_attributes=internal_attributes,
            external_attributes=external_attributes,
            local_header_offset=local_header_offset,
        )
        if _entry_is_symlink(entry):
            raise AcquisitionError(f"symlink ZIP member is forbidden: {name}")
        entries[name] = entry
        cursor = end
    if cursor != len(payload):
        raise AcquisitionError(
            f"central-directory size drift: parsed {cursor}, received {len(payload)}"
        )
    return entries


class MemoryRangeSource:
    """An exact-range source used for verified small archives and tests."""

    def __init__(self, payload: bytes) -> None:
        self.payload = payload
        self.size_bytes = len(payload)
        self.requests: list[tuple[int, int]] = []

    def read_range(self, start: int, end: int) -> bytes:
        if start < 0 or end < start or end >= self.size_bytes:
            raise AcquisitionError(f"invalid memory range {start}-{end}")
        self.requests.append((start, end))
        return self.payload[start : end + 1]


class HTTPSOnlyRedirectHandler(urllib.request.HTTPRedirectHandler):
    """Reject any redirect whose destination is not HTTPS."""

    def redirect_request(
        self,
        req: urllib.request.Request,
        fp: BinaryIO,
        code: int,
        msg: str,
        headers: Any,
        newurl: str,
    ) -> urllib.request.Request | None:
        _require_https(newurl, "redirect URL")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


class HttpsClient:
    """Small-file and exact-range HTTPS client with strict response contracts."""

    def __init__(self) -> None:
        self._opener = urllib.request.build_opener(HTTPSOnlyRedirectHandler())

    def _open(self, request: urllib.request.Request) -> Any:
        try:
            return self._opener.open(request, timeout=60)
        except (OSError, urllib.error.URLError) as error:
            raise AcquisitionError(f"HTTPS request failed: {error}") from error

    def fetch_bytes(self, url: str, expected_size: int, maximum_size: int) -> bytes:
        _require_https(url, "source URL")
        if expected_size <= 0 or expected_size > maximum_size:
            raise AcquisitionError("small-file request exceeds its byte ceiling")
        request = urllib.request.Request(
            url,
            headers={
                "Accept-Encoding": "identity",
                "User-Agent": "vcc-exp111-provenance/1",
            },
        )
        with self._open(request) as response:
            if response.status != 200:
                raise AcquisitionError(
                    f"small-file request returned HTTP {response.status}, expected 200"
                )
            _require_https(response.geturl(), "final response URL")
            encoding = response.headers.get("Content-Encoding")
            if encoding not in {None, "", "identity"}:
                raise AcquisitionError(f"unexpected Content-Encoding: {encoding}")
            raw_length = response.headers.get("Content-Length")
            if raw_length is not None:
                try:
                    header_length = int(raw_length)
                except ValueError as error:
                    raise AcquisitionError(
                        "Content-Length is not an integer"
                    ) from error
                if header_length != expected_size:
                    raise AcquisitionError(
                        f"Content-Length {header_length} != expected {expected_size}"
                    )
            payload = response.read(expected_size + 1)
        if len(payload) != expected_size:
            raise AcquisitionError(
                f"small-file body has {len(payload)} bytes, expected {expected_size}"
            )
        return payload

    def fetch_json(self, url: str, maximum_size: int = 2 << 20) -> dict[str, Any]:
        _require_https(url, "JSON API URL")
        request = urllib.request.Request(
            url,
            headers={
                "Accept": "application/json",
                "Accept-Encoding": "identity",
                "User-Agent": "vcc-exp111-provenance/1",
            },
        )
        with self._open(request) as response:
            if response.status != 200:
                raise AcquisitionError(
                    f"JSON API returned HTTP {response.status}, expected 200"
                )
            _require_https(response.geturl(), "final JSON response URL")
            payload = response.read(maximum_size + 1)
        if len(payload) > maximum_size:
            raise AcquisitionError("JSON API response exceeded its byte ceiling")
        try:
            decoded = json.loads(payload)
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise AcquisitionError(f"JSON API response is invalid: {error}") from error
        if not isinstance(decoded, dict):
            raise AcquisitionError("JSON API response must be an object")
        return decoded

    def read_range(self, url: str, start: int, end: int, total_size: int) -> bytes:
        _require_https(url, "range source URL")
        if start < 0 or end < start or end >= total_size:
            raise AcquisitionError(f"invalid HTTPS byte range {start}-{end}")
        request = urllib.request.Request(
            url,
            headers={
                "Accept-Encoding": "identity",
                "Range": f"bytes={start}-{end}",
                "User-Agent": "vcc-exp111-provenance/1",
            },
        )
        response = self._open(request)
        try:
            # Check status before reading. A server that ignores Range must not stream
            # a multi-gigabyte 200 response into this process.
            if response.status != 206:
                raise AcquisitionError(
                    f"range request returned HTTP {response.status}, expected 206"
                )
            _require_https(response.geturl(), "final range response URL")
            encoding = response.headers.get("Content-Encoding")
            if encoding not in {None, "", "identity"}:
                raise AcquisitionError(f"unexpected Content-Encoding: {encoding}")
            expected_length = end - start + 1
            raw_length = response.headers.get("Content-Length")
            if raw_length is None:
                raise AcquisitionError("range response omitted Content-Length")
            try:
                content_length = int(raw_length)
            except ValueError as error:
                raise AcquisitionError(
                    "range Content-Length is not an integer"
                ) from error
            if content_length != expected_length:
                raise AcquisitionError(
                    f"range Content-Length {content_length} != {expected_length}"
                )
            content_range = response.headers.get("Content-Range")
            match = CONTENT_RANGE_RE.fullmatch(content_range or "")
            if match is None:
                raise AcquisitionError("range response has invalid Content-Range")
            observed = tuple(int(value) for value in match.groups())
            if observed != (start, end, total_size):
                raise AcquisitionError(
                    f"Content-Range {observed} != {(start, end, total_size)}"
                )
            payload = response.read(expected_length + 1)
            if len(payload) != expected_length:
                raise AcquisitionError(
                    f"range body has {len(payload)} bytes, expected {expected_length}"
                )
            return payload
        finally:
            response.close()


class HttpsRangeSource:
    """A bound HTTPS archive that can expose only exact byte ranges."""

    def __init__(self, client: HttpsClient, url: str, size_bytes: int) -> None:
        self.client = client
        self.url = _require_https(url, "archive URL")
        self.size_bytes = _require_positive_int(size_bytes, "archive size")

    def read_range(self, start: int, end: int) -> bytes:
        return self.client.read_range(self.url, start, end, self.size_bytes)


def _entry_matches_spec(entry: ZipEntry, spec: dict[str, Any]) -> bool:
    return (
        entry.name == spec["name"]
        and entry.flag_bits == spec["flag_bits"]
        and entry.compression_method == spec["compression_method"]
        and f"{entry.crc32:08x}" == spec["crc32"]
        and entry.compressed_size == spec["compressed_size"]
        and entry.uncompressed_size == spec["uncompressed_size"]
        and entry.local_header_offset == spec["local_header_offset"]
    )


def verify_selected_entries(
    entries: dict[str, ZipEntry],
    member_specs: list[dict[str, Any]],
) -> dict[str, ZipEntry]:
    expected_names = [spec["name"] for spec in member_specs]
    if len(expected_names) != len(set(expected_names)):
        raise AcquisitionError("selected-member specifications contain duplicates")
    selected: dict[str, ZipEntry] = {}
    for spec in member_specs:
        name = spec["name"]
        entry = entries.get(name)
        if entry is None:
            raise AcquisitionError(f"selected ZIP member is missing: {name}")
        if not _entry_matches_spec(entry, spec):
            raise AcquisitionError(f"central-directory metadata drifted for {name}")
        if name.endswith("/"):
            raise AcquisitionError(f"selected member cannot be a directory: {name}")
        selected[name] = entry
    return selected


def inspect_local_header(
    source: Any, entry: ZipEntry, central_offset: int
) -> LocalMember:
    offset = entry.local_header_offset
    fixed = source.read_range(offset, offset + LOCAL_STRUCT.size - 1)
    (
        signature,
        version_needed,
        flag_bits,
        compression_method,
        _mod_time,
        _mod_date,
        local_crc32,
        local_compressed_size,
        local_uncompressed_size,
        name_length,
        extra_length,
    ) = LOCAL_STRUCT.unpack(fixed)
    if signature != LOCAL_SIGNATURE:
        raise AcquisitionError(f"invalid local-header signature for {entry.name}")
    if version_needed != entry.version_needed:
        raise AcquisitionError(f"local version-needed drift for {entry.name}")
    if flag_bits != entry.flag_bits:
        raise AcquisitionError(f"local flag-bit drift for {entry.name}")
    if compression_method != entry.compression_method:
        raise AcquisitionError(f"local compression-method drift for {entry.name}")
    variable_length = name_length + extra_length
    if variable_length <= 0 or variable_length > (1 << 20):
        raise AcquisitionError(f"unsafe local-header variable length for {entry.name}")
    variable = source.read_range(
        offset + LOCAL_STRUCT.size,
        offset + LOCAL_STRUCT.size + variable_length - 1,
    )
    local_name = _decode_zip_name(variable[:name_length], flag_bits)
    if local_name != entry.name:
        raise AcquisitionError(f"local/central member name mismatch for {entry.name}")
    if not flag_bits & 0x8:
        if (
            local_crc32 != entry.crc32
            or local_compressed_size != entry.compressed_size
            or local_uncompressed_size != entry.uncompressed_size
        ):
            raise AcquisitionError(f"local size/CRC drift for {entry.name}")
    else:
        for observed, expected, field in (
            (local_crc32, entry.crc32, "CRC"),
            (local_compressed_size, entry.compressed_size, "compressed size"),
            (local_uncompressed_size, entry.uncompressed_size, "uncompressed size"),
        ):
            if observed not in {0, expected}:
                raise AcquisitionError(f"local {field} drift for {entry.name}")
    data_offset = offset + LOCAL_STRUCT.size + variable_length
    data_end = data_offset + entry.compressed_size - 1
    if data_end >= central_offset:
        raise AcquisitionError(f"member data overlaps central directory: {entry.name}")
    return LocalMember(entry=entry, data_offset=data_offset, data_end=data_end)


def inspect_zip(
    source: Any,
    central_spec: dict[str, Any],
    member_specs: list[dict[str, Any]],
    *,
    tail_probe_bytes: int,
) -> tuple[dict[str, ZipEntry], dict[str, LocalMember]]:
    archive_size = source.size_bytes
    if central_spec["offset"] + central_spec["size_bytes"] >= archive_size:
        raise AcquisitionError("pinned central directory is outside archive bounds")
    probe_size = min(tail_probe_bytes, archive_size)
    tail_start = archive_size - probe_size
    tail = source.read_range(tail_start, archive_size - 1)
    observed_eocd = parse_eocd(tail, tail_start, archive_size)
    for key in ("entry_count", "offset", "size_bytes"):
        if observed_eocd[key] != central_spec[key]:
            raise AcquisitionError(
                f"ZIP EOCD {key} {observed_eocd[key]} != pinned {central_spec[key]}"
            )
    central = source.read_range(
        central_spec["offset"],
        central_spec["offset"] + central_spec["size_bytes"] - 1,
    )
    actual_central_hash = sha256_bytes(central)
    if actual_central_hash != central_spec["sha256"]:
        raise AcquisitionError(
            f"central-directory SHA-256 {actual_central_hash} != pinned "
            f"{central_spec['sha256']}"
        )
    entries = parse_central_directory(central, central_spec["entry_count"])
    selected = verify_selected_entries(entries, member_specs)
    local_members = {
        name: inspect_local_header(source, entry, central_spec["offset"])
        for name, entry in selected.items()
    }
    intervals = sorted(
        (member.data_offset, member.data_end, name)
        for name, member in local_members.items()
    )
    for previous, current in zip(intervals, intervals[1:], strict=False):
        if previous[1] >= current[0]:
            raise AcquisitionError(
                f"selected compressed intervals overlap: {previous[2]} and {current[2]}"
            )
    return entries, local_members


def extract_member(
    source: Any,
    member: LocalMember,
    output: BinaryIO,
    *,
    max_chunk_bytes: int,
) -> dict[str, Any]:
    entry = member.entry
    if entry.compression_method not in {0, 8}:
        raise AcquisitionError(f"unsupported compression method for {entry.name}")
    digest = hashlib.sha256()
    checksum = 0
    output_size = 0
    decompressor = (
        zlib.decompressobj(-zlib.MAX_WBITS) if entry.compression_method == 8 else None
    )

    cursor = member.data_offset
    while cursor <= member.data_end:
        end = min(member.data_end, cursor + max_chunk_bytes - 1)
        compressed = source.read_range(cursor, end)
        cursor = end + 1
        if decompressor is None:
            pending_outputs = [compressed]
        else:
            pending_outputs = []
            pending = compressed
            while pending:
                remaining = entry.uncompressed_size - output_size
                if remaining < 0:
                    raise AcquisitionError(
                        f"uncompressed size overflow for {entry.name}"
                    )
                try:
                    decoded = decompressor.decompress(
                        pending,
                        min(max_chunk_bytes, remaining + 1),
                    )
                except zlib.error as error:
                    raise AcquisitionError(
                        f"invalid deflate stream for {entry.name}: {error}"
                    ) from error
                pending = decompressor.unconsumed_tail
                pending_outputs.append(decoded)
                if not decoded and pending:
                    raise AcquisitionError(f"deflate stream stalled for {entry.name}")
        for decoded in pending_outputs:
            output_size += len(decoded)
            if output_size > entry.uncompressed_size:
                raise AcquisitionError(f"uncompressed size overflow for {entry.name}")
            checksum = binascii.crc32(decoded, checksum)
            digest.update(decoded)
            output.write(decoded)

    if decompressor is not None:
        try:
            final = decompressor.flush()
        except zlib.error as error:
            raise AcquisitionError(
                f"invalid deflate trailer for {entry.name}: {error}"
            ) from error
        output_size += len(final)
        if output_size > entry.uncompressed_size:
            raise AcquisitionError(f"uncompressed size overflow for {entry.name}")
        checksum = binascii.crc32(final, checksum)
        digest.update(final)
        output.write(final)
        if not decompressor.eof:
            raise AcquisitionError(f"truncated deflate stream for {entry.name}")
        if decompressor.unused_data or decompressor.unconsumed_tail:
            raise AcquisitionError(f"trailing deflate data for {entry.name}")

    checksum &= 0xFFFFFFFF
    if output_size != entry.uncompressed_size:
        raise AcquisitionError(
            f"{entry.name}: extracted {output_size} bytes, "
            f"expected {entry.uncompressed_size}"
        )
    if checksum != entry.crc32:
        raise AcquisitionError(
            f"{entry.name}: CRC32 {checksum:08x} != expected {entry.crc32:08x}"
        )
    return {
        "size_bytes": output_size,
        "crc32": f"{checksum:08x}",
        "sha256": digest.hexdigest(),
    }


def extract_member_bytes(
    source: Any,
    member: LocalMember,
    *,
    max_chunk_bytes: int,
) -> tuple[bytes, dict[str, Any]]:
    output = io.BytesIO()
    report = extract_member(
        source,
        member,
        output,
        max_chunk_bytes=max_chunk_bytes,
    )
    return output.getvalue(), report


def _figshare_observed(article: dict[str, Any]) -> dict[str, Any]:
    license_value = article.get("license")
    license_name = (
        license_value.get("name") if isinstance(license_value, dict) else license_value
    )
    files = article.get("files")
    if not isinstance(files, list):
        raise AcquisitionError("Figshare API article has no file list")
    normalized_files = []
    for file_record in files:
        if not isinstance(file_record, dict):
            raise AcquisitionError("Figshare API file record is malformed")
        normalized_files.append(
            {
                "id": file_record.get("id"),
                "name": file_record.get("name"),
                "size_bytes": file_record.get("size"),
                "md5": file_record.get("computed_md5"),
                "download_url": file_record.get("download_url"),
            }
        )
    return {
        "id": article.get("id"),
        "doi": article.get("doi"),
        "title": article.get("title"),
        "published_date": article.get("published_date"),
        "modified_date": article.get("modified_date"),
        "license": license_name,
        "files": sorted(normalized_files, key=lambda item: item["id"]),
    }


def verify_figshare_article(
    observed_article: dict[str, Any],
    pinned_article: dict[str, Any],
) -> dict[str, Any]:
    observed = _figshare_observed(observed_article)
    expected = {
        "id": pinned_article["id"],
        "doi": pinned_article["doi"],
        "title": pinned_article["title"],
        "published_date": pinned_article["published_date"],
        "modified_date": pinned_article["modified_date"],
        "license": pinned_article.get("license"),
        "files": sorted(
            [
                {
                    "id": item["id"],
                    "name": item["name"],
                    "size_bytes": item["size_bytes"],
                    "md5": item["md5"],
                    "download_url": item["download_url"],
                }
                for item in pinned_article["files"]
            ],
            key=lambda item: item["id"],
        ),
    }
    if observed != expected:
        raise AcquisitionError("live Figshare article metadata differs from the pin")
    return {
        "article_id": observed["id"],
        "normalized_metadata_sha256": canonical_json_sha256(observed),
        "file_count": len(observed["files"]),
    }


def audit_zipfile_members(archive: zipfile.ZipFile) -> list[zipfile.ZipInfo]:
    seen: set[str] = set()
    infos = archive.infolist()
    for info in infos:
        safe_zip_member(info.filename)
        if info.filename in seen:
            raise AcquisitionError(f"duplicate ZIP member path: {info.filename}")
        seen.add(info.filename)
        if info.flag_bits & 0x1:
            raise AcquisitionError(
                f"encrypted ZIP member is forbidden: {info.filename}"
            )
        unix_mode = (info.external_attr >> 16) & 0xFFFF
        if info.create_system == 3 and stat.S_ISLNK(unix_mode):
            raise AcquisitionError(f"symlink ZIP member is forbidden: {info.filename}")
    return infos


def normalized_code_tree_sha256(payload: bytes) -> dict[str, Any]:
    try:
        archive = zipfile.ZipFile(io.BytesIO(payload))
    except zipfile.BadZipFile as error:
        raise AcquisitionError(f"code archive is not a valid ZIP: {error}") from error
    with archive:
        infos = audit_zipfile_members(archive)
        normalized: dict[str, zipfile.ZipInfo] = {}
        total_uncompressed = 0
        roots: set[str] = set()
        for info in infos:
            path = safe_zip_member(info.filename)
            if "__MACOSX" in path.parts or path.name == ".DS_Store":
                continue
            if path.name.startswith("._"):
                continue
            if info.is_dir():
                continue
            if len(path.parts) < 2:
                raise AcquisitionError("code archive member lacks one common root")
            roots.add(path.parts[0])
            relative = "/".join(path.parts[1:])
            safe_zip_member(relative)
            if relative in normalized:
                raise AcquisitionError(f"duplicate normalized code path: {relative}")
            if info.file_size > 5 << 20:
                raise AcquisitionError(f"code member exceeds 5 MiB: {relative}")
            total_uncompressed += info.file_size
            if total_uncompressed > 20 << 20:
                raise AcquisitionError("code archive exceeds 20 MiB uncompressed")
            normalized[relative] = info
        if len(roots) != 1 or not normalized:
            raise AcquisitionError("code archive must have one nonempty source root")
        digest = hashlib.sha256()
        for relative in sorted(normalized):
            digest.update(relative.encode("utf-8"))
            digest.update(b"\0")
            try:
                with archive.open(normalized[relative]) as handle:
                    while block := handle.read(1 << 20):
                        digest.update(block)
            except (OSError, RuntimeError, zipfile.BadZipFile) as error:
                raise AcquisitionError(
                    f"cannot read code archive member {relative}: {error}"
                ) from error
            digest.update(b"\0")
    return {
        "sha256": digest.hexdigest(),
        "file_count": len(normalized),
        "uncompressed_size_bytes": total_uncompressed,
    }


def _asset_by_name(manifest: dict[str, Any], name: str) -> dict[str, Any]:
    matches = [
        record
        for record in manifest["figshare"]["data_article"]["files"]
        if record["name"] == name
    ]
    if len(matches) != 1:
        raise AcquisitionError(f"manifest does not contain exactly one {name}")
    return matches[0]


def _verify_small_archive(
    payload: bytes,
    expected_size: int,
    expected_md5: str,
    expected_sha256: str,
    label: str,
) -> dict[str, Any]:
    actual_md5 = hashlib.md5(payload, usedforsecurity=False).hexdigest()
    actual_sha256 = sha256_bytes(payload)
    if len(payload) != expected_size:
        raise AcquisitionError(f"{label} size drifted")
    if actual_md5 != expected_md5:
        raise AcquisitionError(f"{label} MD5 drifted")
    if actual_sha256 != expected_sha256:
        raise AcquisitionError(f"{label} SHA-256 drifted")
    return {
        "size_bytes": len(payload),
        "md5": actual_md5,
        "sha256": actual_sha256,
    }


def _validate_cell_line_metadata_bytes(
    payload: bytes,
    manifest: dict[str, Any],
) -> dict[str, Any]:
    try:
        text = io.StringIO(payload.decode("utf-8"), newline="")
    except UnicodeDecodeError as error:
        raise AcquisitionError("cell-line metadata is not UTF-8") from error
    reader = csv.DictReader(text)
    if tuple(reader.fieldnames or ()) != METADATA_COLUMNS:
        raise AcquisitionError("cell-line metadata columns drifted")
    rows = list(reader)
    if len(rows) != 16:
        raise AcquisitionError("cell-line metadata must contain exactly 16 rows")
    observed = [
        {
            "cell_line": row["cell_line"],
            "arxspan_id": row["arxspan_id"],
            "batch": row["batch"],
            "lineage": row["OncotreeLineage"],
        }
        for row in rows
    ]
    if observed != manifest["dataset_contract"]["cell_lines"]:
        raise AcquisitionError("cell-line metadata records/order drifted")
    if len({row["lineage"] for row in observed}) != 11:
        raise AcquisitionError("cell-line metadata must contain 11 lineages")
    if len({row["batch"] for row in observed}) != 6:
        raise AcquisitionError("cell-line metadata must contain six batches")
    return {
        "row_count": len(rows),
        "lineage_count": len({row["lineage"] for row in observed}),
        "batch_count": len({row["batch"] for row in observed}),
        "sha256": sha256_bytes(payload),
    }


def _metadata_from_archive(
    payload: bytes,
    manifest: dict[str, Any],
) -> tuple[bytes, dict[str, Any]]:
    asset = _asset_by_name(manifest, "metadata.zip")
    archive_report = _verify_small_archive(
        payload,
        asset["size_bytes"],
        asset["md5"],
        asset["sha256"],
        "metadata.zip",
    )
    source = MemoryRangeSource(payload)
    archive_spec = manifest["archives"]["metadata_zip"]
    entries, members = inspect_zip(
        source,
        archive_spec["central_directory"],
        [archive_spec["selected_member"]],
        tail_probe_bytes=min(131072, len(payload)),
    )
    selected = members[ALLOWED_MEMBERS[0]]
    metadata_bytes, extracted_report = extract_member_bytes(
        source,
        selected,
        max_chunk_bytes=manifest["safety"]["max_range_chunk_bytes"],
    )
    member_spec = archive_spec["selected_member"]
    if (
        hashlib.md5(metadata_bytes, usedforsecurity=False).hexdigest()
        != member_spec["content_md5"]
    ):
        raise AcquisitionError("cell-line metadata content MD5 drifted")
    if sha256_bytes(metadata_bytes) != member_spec["content_sha256"]:
        raise AcquisitionError("cell-line metadata content SHA-256 drifted")
    schema_report = _validate_cell_line_metadata_bytes(metadata_bytes, manifest)
    return metadata_bytes, {
        "archive": archive_report,
        "central_directory_sha256": archive_spec["central_directory"]["sha256"],
        "entry_count": len(entries),
        "member": {
            **extracted_report,
            "data_offset": selected.data_offset,
            "schema": schema_report,
        },
    }


def _article_file_binding(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "id": record["id"],
        "name": record["name"],
        "size_bytes": record["size_bytes"],
        "md5": record["md5"],
        "download_url": record["download_url"],
        "sha256": record.get("sha256"),
    }


def source_binding(
    manifest: dict[str, Any],
    protocol: dict[str, Any],
    manifest_sha256: str,
    protocol_sha256: str,
) -> dict[str, Any]:
    return {
        "experiment": protocol["experiment"],
        "manifest_sha256": manifest_sha256,
        "protocol_sha256": protocol_sha256,
        "figshare_data_article": {
            "id": manifest["figshare"]["data_article"]["id"],
            "doi": manifest["figshare"]["data_article"]["doi"],
            "modified_date": manifest["figshare"]["data_article"]["modified_date"],
            "files": [
                _article_file_binding(record)
                for record in manifest["figshare"]["data_article"]["files"]
            ],
        },
        "figshare_code_article": {
            "id": manifest["figshare"]["code_article"]["id"],
            "doi": manifest["figshare"]["code_article"]["doi"],
            "modified_date": manifest["figshare"]["code_article"]["modified_date"],
            "files": [
                _article_file_binding(record)
                for record in manifest["figshare"]["code_article"]["files"]
            ],
        },
        "github": manifest["github"],
        "archives": manifest["archives"],
        "challenge_contract": manifest["challenge_contract"],
        "allowed_members": list(ALLOWED_MEMBERS),
    }


def audit_remote(
    manifest: dict[str, Any],
    protocol: dict[str, Any],
    *,
    client: HttpsClient | None = None,
) -> RemoteAudit:
    client = client or HttpsClient()
    data_pin = manifest["figshare"]["data_article"]
    code_pin = manifest["figshare"]["code_article"]
    data_live = client.fetch_json(data_pin["api_url"])
    code_live = client.fetch_json(code_pin["api_url"])
    data_report = verify_figshare_article(data_live, data_pin)
    code_report = verify_figshare_article(code_live, code_pin)

    code_asset = code_pin["files"][0]
    figshare_code = client.fetch_bytes(
        code_asset["download_url"],
        code_asset["size_bytes"],
        manifest["safety"]["max_full_download_bytes"],
    )
    figshare_code_archive = _verify_small_archive(
        figshare_code,
        code_asset["size_bytes"],
        code_asset["md5"],
        code_asset["sha256"],
        "Figshare code archive",
    )
    figshare_tree = normalized_code_tree_sha256(figshare_code)
    if figshare_tree["sha256"] != code_asset["normalized_tree_sha256"]:
        raise AcquisitionError("Figshare code normalized-tree SHA-256 drifted")

    github_pin = manifest["github"]["commit_archive"]
    github_code = client.fetch_bytes(
        github_pin["url"],
        github_pin["size_bytes"],
        manifest["safety"]["max_full_download_bytes"],
    )
    github_code_archive = _verify_small_archive(
        github_code,
        github_pin["size_bytes"],
        github_pin["md5"],
        github_pin["sha256"],
        "GitHub commit archive",
    )
    github_tree = normalized_code_tree_sha256(github_code)
    if github_tree["sha256"] != github_pin["normalized_tree_sha256"]:
        raise AcquisitionError("GitHub code normalized-tree SHA-256 drifted")
    if github_tree["sha256"] != figshare_tree["sha256"]:
        raise AcquisitionError("live GitHub/Figshare normalized code trees disagree")

    metadata_asset = _asset_by_name(manifest, "metadata.zip")
    metadata_zip = client.fetch_bytes(
        metadata_asset["download_url"],
        metadata_asset["size_bytes"],
        manifest["safety"]["max_full_download_bytes"],
    )
    metadata_member, metadata_report = _metadata_from_archive(metadata_zip, manifest)

    processed_asset = _asset_by_name(manifest, "processed.zip")
    processed_source = HttpsRangeSource(
        client,
        processed_asset["download_url"],
        processed_asset["size_bytes"],
    )
    processed_spec = manifest["archives"]["processed_zip"]
    entries, processed_members = inspect_zip(
        processed_source,
        processed_spec["central_directory"],
        processed_spec["selected_members"],
        tail_probe_bytes=processed_spec["central_directory"]["tail_probe_bytes"],
    )
    report = {
        "status": "pass",
        "figshare": {
            "data_article": data_report,
            "code_article": code_report,
        },
        "code": {
            "figshare_archive": figshare_code_archive,
            "github_archive": github_code_archive,
            "normalized_tree": github_tree,
        },
        "metadata_zip": metadata_report,
        "processed_zip": {
            "archive_size_bytes": processed_source.size_bytes,
            "published_md5": processed_asset["md5"],
            "published_md5_verified": False,
            "central_directory_sha256": processed_spec["central_directory"]["sha256"],
            "entry_count": len(entries),
            "selected_members": {
                name: {
                    **asdict(local.entry),
                    "crc32": f"{local.entry.crc32:08x}",
                    "data_offset": local.data_offset,
                    "data_end": local.data_end,
                }
                for name, local in processed_members.items()
            },
        },
        "authorization": {
            "remote_metadata_audit": protocol["current_authorization"][
                "remote_metadata_audit"
            ],
            "selective_member_acquisition": protocol["current_authorization"][
                "selective_member_acquisition"
            ],
        },
    }
    return RemoteAudit(
        report=report,
        metadata_member=metadata_member,
        processed_source=processed_source,
        processed_members=processed_members,
    )


def _iter_csv_rows(path: Path) -> Iterator[list[str]]:
    try:
        handle = path.open("r", encoding="utf-8", newline="")
    except OSError as error:
        raise AcquisitionError(f"cannot open CSV {path}: {error}") from error
    with handle:
        try:
            yield from csv.reader(handle)
        except (csv.Error, UnicodeDecodeError) as error:
            raise AcquisitionError(f"cannot parse CSV {path}: {error}") from error


def _parse_finite(value: str, label: str) -> float:
    if value == "":
        raise AcquisitionError(f"{label} is empty")
    try:
        parsed = float(value)
    except ValueError as error:
        raise AcquisitionError(f"{label} is not numeric: {value!r}") from error
    if not math.isfinite(parsed):
        raise AcquisitionError(f"{label} is not finite")
    return parsed


def _validate_matrix_csv(
    path: Path,
    expected_pairs: list[tuple[str, str]],
    *,
    value_kind: str,
    minimum_response_rows: int = 10000,
) -> dict[str, Any]:
    rows = _iter_csv_rows(path)
    try:
        outer = next(rows)
        inner = next(rows)
    except StopIteration as error:
        raise AcquisitionError(f"{path.name} is missing its two-row header") from error
    expected_width = len(expected_pairs) + 1
    if len(outer) != expected_width or len(inner) != expected_width:
        raise AcquisitionError(f"{path.name} column count drifted")
    observed_pairs = list(zip(outer[1:], inner[1:], strict=True))
    if observed_pairs != expected_pairs:
        raise AcquisitionError(f"{path.name} context-target axis/order drifted")

    try:
        first_data = next(rows)
    except StopIteration as error:
        raise AcquisitionError(f"{path.name} has no response rows") from error
    index_label = None
    if (
        len(first_data) == expected_width
        and first_data[0]
        and all(value == "" for value in first_data[1:])
    ):
        index_label = first_data[0]
        try:
            first_data = next(rows)
        except StopIteration as error:
            raise AcquisitionError(f"{path.name} has no response rows") from error

    response_digest = hashlib.sha256()
    missing_digest = hashlib.sha256()
    response_ids: list[str] = []
    seen: set[str] = set()
    row_count = 0
    missing_count = 0
    for row in itertools.chain((first_data,), rows):
        row_count += 1
        if len(row) != expected_width:
            raise AcquisitionError(f"{path.name} row {row_count} width drifted")
        response_id = row[0]
        if not response_id or response_id in seen:
            raise AcquisitionError(
                f"{path.name} has an empty/duplicate response id at row {row_count}"
            )
        seen.add(response_id)
        response_ids.append(response_id)
        response_digest.update(response_id.encode("utf-8"))
        response_digest.update(b"\n")
        missing_mask = bytearray()
        nonmissing_count = 0
        for column_index, raw_value in enumerate(row[1:], start=1):
            if raw_value == "":
                missing_mask.append(1)
                missing_count += 1
                continue
            missing_mask.append(0)
            nonmissing_count += 1
            value = _parse_finite(
                raw_value,
                f"{path.name} row {row_count} column {column_index}",
            )
            if value_kind == "count" and (value < 0 or not value.is_integer()):
                raise AcquisitionError(
                    f"{path.name} contains a non-integer/negative count"
                )
            if value_kind == "fdr" and not 0.0 <= value <= 1.0:
                raise AcquisitionError(f"{path.name} contains FDR outside [0,1]")
        if nonmissing_count == 0:
            raise AcquisitionError(f"{path.name} has an entirely missing response row")
        missing_digest.update(missing_mask)
    if row_count < minimum_response_rows:
        raise AcquisitionError(
            f"{path.name} has fewer than {minimum_response_rows:,} response genes"
        )
    return {
        "column_count": len(expected_pairs),
        "row_count": row_count,
        "response_axis_sha256": response_digest.hexdigest(),
        "missing_mask_sha256": missing_digest.hexdigest(),
        "missing_value_count": missing_count,
        "response_ids": response_ids,
        "index_label": index_label,
    }


def _validate_metadata_csv(path: Path, manifest: dict[str, Any]) -> dict[str, Any]:
    payload = path.read_bytes()
    spec = manifest["archives"]["metadata_zip"]["selected_member"]
    actual_md5 = hashlib.md5(payload, usedforsecurity=False).hexdigest()
    actual_sha256 = sha256_bytes(payload)
    if len(payload) != spec["uncompressed_size"]:
        raise AcquisitionError("cell-line metadata extracted size drifted")
    if actual_md5 != spec["content_md5"] or actual_sha256 != spec["content_sha256"]:
        raise AcquisitionError("cell-line metadata extracted hash drifted")
    return _validate_cell_line_metadata_bytes(payload, manifest)


def _validate_crispr_csv(
    path: Path,
    manifest: dict[str, Any],
    recovery: dict[str, Any],
) -> dict[str, Any]:
    rows = _iter_csv_rows(path)
    try:
        header = next(rows)
    except StopIteration as error:
        raise AcquisitionError("crispr_table.csv is empty") from error
    if tuple(header) != CRISPR_COLUMNS:
        raise AcquisitionError("crispr_table.csv columns drifted")
    expected_lines = [
        row["cell_line"] for row in manifest["dataset_contract"]["cell_lines"]
    ]
    expected_targets = manifest["dataset_contract"]["targets"]
    expected_pairs = {
        (line, target) for line in expected_lines for target in expected_targets
    }
    observed_pairs: set[tuple[str, str]] = set()
    missing_pairs: list[tuple[str, str]] = []
    controls = set(manifest["dataset_contract"]["control_targets"])
    classes = {"Strong", "Moderate", "Weak", "Nondependent"}
    for row_number, row in enumerate(rows, start=2):
        if len(row) != len(CRISPR_COLUMNS):
            raise AcquisitionError(f"crispr_table.csv row {row_number} width drifted")
        cell_line, guide, raw_effect, strength, raw_dependent = row
        pair = (cell_line, guide)
        if pair not in expected_pairs or pair in observed_pairs:
            raise AcquisitionError(f"crispr_table.csv pair drift at row {row_number}")
        observed_pairs.add(pair)
        if raw_effect == "" or strength == "":
            if (
                raw_effect != ""
                or strength != ""
                or guide not in controls
                or raw_dependent != "False"
            ):
                raise AcquisitionError(
                    "crispr_table.csv has unauthorized metadata missingness "
                    f"at row {row_number}"
                )
            missing_pairs.append(pair)
            continue
        _parse_finite(raw_effect, f"crispr_table.csv row {row_number} gene_effect")
        if strength not in classes:
            raise AcquisitionError(
                f"crispr_table.csv row {row_number} dependency class drifted"
            )
        if raw_dependent not in {"True", "False"}:
            raise AcquisitionError(
                f"crispr_table.csv row {row_number} boolean encoding drifted"
            )
        if (strength == "Strong") != (raw_dependent == "True"):
            raise AcquisitionError(
                f"crispr_table.csv row {row_number} dependency fields disagree"
            )
    if observed_pairs != expected_pairs:
        raise AcquisitionError("crispr_table.csv is not the exact 16x100 grid")
    exception = recovery["crispr_table_exception"]
    missing_digest = hashlib.sha256()
    for cell_line, guide in sorted(missing_pairs):
        missing_digest.update(cell_line.encode("utf-8"))
        missing_digest.update(b"\0")
        missing_digest.update(guide.encode("utf-8"))
        missing_digest.update(b"\n")
    if (
        len(missing_pairs) != exception["missing_row_count"]
        or missing_digest.hexdigest()
        != exception["sorted_cell_line_guide_pairs_sha256"]
    ):
        raise AcquisitionError("crispr_table.csv control missingness pattern drifted")
    if (
        path.stat().st_size != exception["member_size_bytes"]
        or crc32_file(path) != exception["member_crc32"]
        or sha256_file(path) != exception["member_sha256"]
    ):
        raise AcquisitionError("crispr_table.csv recovery member binding drifted")
    return {
        "row_count": len(observed_pairs),
        "cell_line_count": len({pair[0] for pair in observed_pairs}),
        "target_count": len({pair[1] for pair in observed_pairs}),
        "control_metadata_missing_row_count": len(missing_pairs),
        "control_metadata_missing_pairs_sha256": missing_digest.hexdigest(),
    }


def _validate_quality_csv(path: Path, manifest: dict[str, Any]) -> dict[str, Any]:
    rows = _iter_csv_rows(path)
    try:
        header = next(rows)
    except StopIteration as error:
        raise AcquisitionError("cell_quality_covariates.csv is empty") from error
    if tuple(header) != QUALITY_COLUMNS:
        raise AcquisitionError("cell_quality_covariates.csv columns drifted")
    expected_lines = sorted(
        record["cell_line"] for record in manifest["dataset_contract"]["cell_lines"]
    )
    observed_lines: list[str] = []
    for row_number, row in enumerate(rows, start=2):
        if len(row) != len(QUALITY_COLUMNS):
            raise AcquisitionError(
                f"cell_quality_covariates.csv row {row_number} width drifted"
            )
        line = row[0]
        if line in observed_lines:
            raise AcquisitionError(f"duplicate cell-quality row for {line}")
        observed_lines.append(line)
        values = [
            _parse_finite(
                raw_value,
                f"cell_quality_covariates.csv row {row_number} column {index}",
            )
            for index, raw_value in enumerate(row[1:], start=1)
        ]
        if any(value < 0 for value in values):
            raise AcquisitionError("cell-quality covariates must be nonnegative")
        if not 0.0 <= values[3] <= 1.0:
            raise AcquisitionError("frac_self_downregulated must be in [0,1]")
    if observed_lines != expected_lines:
        raise AcquisitionError("cell-quality cell-line axis/order drifted")
    return {"row_count": len(observed_lines)}


def _expected_matrix_pairs(
    manifest: dict[str, Any],
) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    contract = manifest["dataset_contract"]
    lines = [record["cell_line"] for record in contract["cell_lines"]]
    targets = contract["targets"]
    controls = set(contract["control_targets"])
    missing = {
        tuple(pair) for pair in contract["pseudobulk"]["missing_context_target_pairs"]
    }
    pseudobulk = [
        (line, target)
        for line in lines
        for target in targets
        if (line, target) not in missing
    ]
    sceptre = [
        (line, target) for line in lines for target in targets if target not in controls
    ]
    return pseudobulk, sceptre


def _audit_extracted_tree(root: Path) -> dict[str, Path]:
    root = require_safe_extraction_root(root)
    if root.is_symlink() or not root.is_dir():
        raise AcquisitionError(f"extraction root is missing or symlinked: {root}")
    found: dict[str, Path] = {}
    for current, directories, files in os.walk(root, followlinks=False):
        current_path = Path(current)
        for directory in directories:
            child = current_path / directory
            if child.is_symlink():
                raise AcquisitionError(f"symlinked extraction directory: {child}")
        for filename in files:
            child = current_path / filename
            if child.is_symlink() or not child.is_file():
                raise AcquisitionError(f"non-regular extracted path: {child}")
            relative = child.relative_to(root).as_posix()
            safe_zip_member(relative)
            if relative in found:
                raise AcquisitionError(f"duplicate extracted path: {relative}")
            found[relative] = child
    if set(found) != set(ALLOWED_MEMBERS):
        raise AcquisitionError(
            "extracted file allowlist drifted: "
            f"missing={sorted(set(ALLOWED_MEMBERS) - set(found))}, "
            f"extra={sorted(set(found) - set(ALLOWED_MEMBERS))}"
        )
    return found


def validate_csv_schemas(
    root: Path,
    manifest: dict[str, Any],
    recovery: dict[str, Any],
) -> dict[str, Any]:
    files = _audit_extracted_tree(root)
    metadata = _validate_metadata_csv(files[ALLOWED_MEMBERS[0]], manifest)
    crispr = _validate_crispr_csv(
        files["processed/crispr_table.csv"],
        manifest,
        recovery,
    )
    quality = _validate_quality_csv(
        files["processed/cell_quality_covariates.csv"],
        manifest,
    )
    pseudobulk_pairs, sceptre_pairs = _expected_matrix_pairs(manifest)
    pseudobulk = _validate_matrix_csv(
        files["processed/pseudobulk_sum_by_ko.csv"],
        pseudobulk_pairs,
        value_kind="count",
    )
    lfc = _validate_matrix_csv(
        files["processed/lfc_matrix.csv"],
        sceptre_pairs,
        value_kind="lfc",
    )
    fdr = _validate_matrix_csv(
        files["processed/fdr_matrix.csv"],
        sceptre_pairs,
        value_kind="fdr",
    )
    if lfc["response_ids"] != fdr["response_ids"]:
        raise AcquisitionError("LFC/FDR response-gene axes are not bit-identical")
    if lfc["missing_mask_sha256"] != fdr["missing_mask_sha256"]:
        raise AcquisitionError("LFC/FDR missing-value masks are not bit-identical")
    del pseudobulk["response_ids"]
    del lfc["response_ids"]
    del fdr["response_ids"]
    return {
        "metadata": metadata,
        "crispr_table": crispr,
        "cell_quality_covariates": quality,
        "pseudobulk_sum_by_ko": pseudobulk,
        "lfc_matrix": lfc,
        "fdr_matrix": fdr,
    }


def _member_receipt_binding(
    name: str,
    manifest: dict[str, Any],
) -> dict[str, Any]:
    if name == ALLOWED_MEMBERS[0]:
        asset = _asset_by_name(manifest, "metadata.zip")
        member = manifest["archives"]["metadata_zip"]["selected_member"]
        central = manifest["archives"]["metadata_zip"]["central_directory"]
    else:
        asset = _asset_by_name(manifest, "processed.zip")
        members = {
            item["name"]: item
            for item in manifest["archives"]["processed_zip"]["selected_members"]
        }
        member = members[name]
        central = manifest["archives"]["processed_zip"]["central_directory"]
    return {
        "member_name": name,
        "archive_file_id": asset["id"],
        "archive_name": asset["name"],
        "archive_size_bytes": asset["size_bytes"],
        "archive_published_md5": asset["md5"],
        "archive_sha256": asset.get("sha256"),
        "central_directory_sha256": central["sha256"],
        "pinned_member": member,
    }


def build_receipt(
    root: Path,
    manifest: dict[str, Any],
    protocol: dict[str, Any],
    remote_report: dict[str, Any],
    local_members: dict[str, LocalMember],
    manifest_sha256: str,
    protocol_sha256: str,
    authorization_sha256: str,
    authorization: dict[str, Any],
) -> dict[str, Any]:
    files = _audit_extracted_tree(root)
    members: dict[str, Any] = {}
    for name in ALLOWED_MEMBERS:
        path = files[name]
        binding = _member_receipt_binding(name, manifest)
        if name == ALLOWED_MEMBERS[0]:
            data_offset = remote_report["metadata_zip"]["member"]["data_offset"]
        else:
            data_offset = local_members[name].data_offset
        members[name] = {
            "binding": binding,
            "relative_path": name,
            "data_offset": data_offset,
            "size_bytes": path.stat().st_size,
            "crc32": crc32_file(path),
            "sha256": sha256_file(path),
        }
    return {
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "status": "complete",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "source_binding": source_binding(
            manifest,
            protocol,
            manifest_sha256,
            protocol_sha256,
        ),
        "execution_authorization": {
            "path": str(DEFAULT_AUTHORIZATION.relative_to(REPO_ROOT)),
            "sha256": authorization_sha256,
            "status": authorization["status"],
            "authorized_actions": authorization["authorized_actions"],
            "forbidden_actions": authorization["forbidden_actions"],
            "output_policy": authorization["output_policy"],
        },
        "remote_audit": remote_report,
        "members": members,
    }


def verify_receipt(
    receipt_path: Path,
    extracted_root: Path,
    manifest: dict[str, Any],
    protocol: dict[str, Any],
    manifest_sha256: str,
    protocol_sha256: str,
    authorization_sha256: str,
    authorization: dict[str, Any],
) -> dict[str, Any]:
    receipt_path = require_under_experiment(receipt_path)
    if receipt_path.is_symlink() or not receipt_path.is_file():
        raise AcquisitionError(f"receipt is missing or symlinked: {receipt_path}")
    receipt = load_json(receipt_path, "acquisition receipt")
    if receipt.get("schema_version") != RECEIPT_SCHEMA_VERSION:
        raise AcquisitionError("unsupported acquisition-receipt schema")
    if receipt.get("status") != "complete":
        raise AcquisitionError("acquisition receipt is not complete")
    expected_binding = source_binding(
        manifest,
        protocol,
        manifest_sha256,
        protocol_sha256,
    )
    if receipt.get("source_binding") != expected_binding:
        raise AcquisitionError("receipt source binding is stale or tampered")
    expected_authorization = {
        "path": str(DEFAULT_AUTHORIZATION.relative_to(REPO_ROOT)),
        "sha256": authorization_sha256,
        "status": authorization["status"],
        "authorized_actions": authorization["authorized_actions"],
        "forbidden_actions": authorization["forbidden_actions"],
        "output_policy": authorization["output_policy"],
    }
    if receipt.get("execution_authorization") != expected_authorization:
        raise AcquisitionError(
            "receipt execution-authorization binding is stale or tampered"
        )
    members = receipt.get("members")
    if not isinstance(members, dict) or set(members) != set(ALLOWED_MEMBERS):
        raise AcquisitionError("receipt selective-member allowlist drifted")
    files = _audit_extracted_tree(extracted_root)
    for name in ALLOWED_MEMBERS:
        record = members[name]
        if not isinstance(record, dict):
            raise AcquisitionError(f"receipt member record is malformed: {name}")
        if record.get("binding") != _member_receipt_binding(name, manifest):
            raise AcquisitionError(f"receipt member binding drifted: {name}")
        if record.get("relative_path") != name:
            raise AcquisitionError(f"receipt member path drifted: {name}")
        path = files[name]
        actual_size = path.stat().st_size
        actual_crc = crc32_file(path)
        actual_sha256 = sha256_file(path)
        if (
            record.get("size_bytes") != actual_size
            or record.get("crc32") != actual_crc
            or record.get("sha256") != actual_sha256
        ):
            raise AcquisitionError(f"receipt/file hash or size mismatch: {name}")
        pinned = record["binding"]["pinned_member"]
        if actual_size != pinned["uncompressed_size"]:
            raise AcquisitionError(f"receipt/member uncompressed size drifted: {name}")
        if actual_crc != pinned["crc32"]:
            raise AcquisitionError(f"receipt/member CRC32 drifted: {name}")
        if not isinstance(record.get("data_offset"), int) or record["data_offset"] <= 0:
            raise AcquisitionError(f"receipt member data offset is invalid: {name}")
    return receipt


@contextmanager
def acquisition_lock() -> Iterator[None]:
    downloads = HERE / "downloads"
    downloads.mkdir(parents=True, exist_ok=True)
    if downloads.is_symlink():
        raise AcquisitionError("downloads directory cannot be a symlink")
    lock_path = downloads / ".acquisition.lock"
    with lock_path.open("a+", encoding="utf-8") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise AcquisitionError(
                "another Experiment 111 acquisition is active"
            ) from error
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            lock_path.unlink(missing_ok=True)


def _publish_extraction(
    remote: RemoteAudit,
    manifest: dict[str, Any],
    protocol: dict[str, Any],
    *,
    extracted_root: Path,
    audit_root: Path,
    manifest_sha256: str,
    protocol_sha256: str,
    authorization_sha256: str,
    authorization: dict[str, Any],
    recovery: dict[str, Any],
) -> dict[str, Any]:
    extracted_root = require_exact_source_root(extracted_root)
    audit_root = require_under_experiment(audit_root)
    if extracted_root.exists() or extracted_root.is_symlink():
        raise AcquisitionError(
            f"refusing to overwrite existing extraction root: {extracted_root}"
        )
    if audit_root.is_symlink():
        raise AcquisitionError("audit root cannot be a symlink")
    stage = Path(
        tempfile.mkdtemp(
            dir=HERE,
            prefix=".extracted.stage-",
        )
    )
    try:
        metadata_path = stage / ALLOWED_MEMBERS[0]
        metadata_path.parent.mkdir(parents=True, exist_ok=True)
        metadata_path.write_bytes(remote.metadata_member)

        member_reports: dict[str, Any] = {}
        chunk_bytes = manifest["safety"]["max_range_chunk_bytes"]
        for name in PROCESSED_MEMBERS:
            destination = stage / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            with destination.open("xb") as output:
                member_reports[name] = extract_member(
                    remote.processed_source,
                    remote.processed_members[name],
                    output,
                    max_chunk_bytes=chunk_bytes,
                )
                output.flush()
                os.fsync(output.fileno())
        schema_report = validate_csv_schemas(stage, manifest, recovery)
        os.replace(stage, extracted_root)
        receipt = build_receipt(
            extracted_root,
            manifest,
            protocol,
            remote.report,
            remote.processed_members,
            manifest_sha256,
            protocol_sha256,
            authorization_sha256,
            authorization,
        )
        atomic_write_json(audit_root / RECEIPT_NAME, receipt)
        return {
            "status": "acquired",
            "extracted_root": str(extracted_root),
            "receipt": str(audit_root / RECEIPT_NAME),
            "members": member_reports,
            "schema": schema_report,
        }
    finally:
        if stage.exists():
            shutil.rmtree(stage)


def run_manifest(
    manifest_path: Path = DEFAULT_MANIFEST,
    protocol_path: Path = DEFAULT_PROTOCOL,
) -> dict[str, Any]:
    manifest, protocol = load_contracts(manifest_path, protocol_path)
    authorization_status = False
    authorization_sha256 = None
    if DEFAULT_AUTHORIZATION.exists():
        authorization = validate_execution_authorization(
            load_json(DEFAULT_AUTHORIZATION, "execution authorization"),
            manifest_path=manifest_path,
            protocol_path=protocol_path,
        )
        authorization_status = authorization["authorized_actions"][
            "selective_member_acquisition"
        ]
        authorization_sha256 = sha256_file(DEFAULT_AUTHORIZATION)
    processed = manifest["archives"]["processed_zip"]["selected_members"]
    return {
        "status": "pass",
        "experiment": protocol["experiment"],
        "manifest_sha256": sha256_file(manifest_path),
        "protocol_sha256": sha256_file(protocol_path),
        "cell_lines": manifest["dataset_contract"]["expected_cell_line_count"],
        "lineages": manifest["dataset_contract"]["expected_lineage_count"],
        "targets": manifest["dataset_contract"]["expected_target_count"],
        "vcc_panel_overlap": manifest["challenge_contract"]["expected_panel_overlap"],
        "selected_members": list(ALLOWED_MEMBERS),
        "selected_compressed_bytes": sum(
            member["compressed_size"] for member in processed
        ),
        "selective_member_acquisition_authorized": protocol["current_authorization"][
            "selective_member_acquisition"
        ],
        "execution_amendment_selective_acquisition_authorized": authorization_status,
        "execution_authorization_sha256": authorization_sha256,
    }


def run_remote(
    manifest_path: Path = DEFAULT_MANIFEST,
    protocol_path: Path = DEFAULT_PROTOCOL,
) -> dict[str, Any]:
    manifest, protocol = load_contracts(manifest_path, protocol_path)
    if not protocol["current_authorization"]["remote_metadata_audit"]:
        raise AcquisitionError(
            "remote metadata audit is not preregistered as authorized"
        )
    return audit_remote(manifest, protocol).report


def run_acquire(
    *,
    confirm_selective_download: bool,
    dry_run: bool,
    manifest_path: Path = DEFAULT_MANIFEST,
    protocol_path: Path = DEFAULT_PROTOCOL,
    authorization_path: Path = DEFAULT_AUTHORIZATION,
    extracted_root: Path = DEFAULT_EXTRACTED,
    audit_root: Path = DEFAULT_AUDIT,
) -> dict[str, Any]:
    manifest, protocol, authorization = load_execution_contracts(
        manifest_path,
        protocol_path,
        authorization_path,
    )
    processed = manifest["archives"]["processed_zip"]["selected_members"]
    plan = {
        "status": "dry_run",
        "network_requests_per_processed_member": "one local-header pair plus "
        "ceil(compressed_size / 8388608) payload ranges",
        "selected_members": [
            {
                "name": member["name"],
                "local_header_offset": member["local_header_offset"],
                "compressed_size": member["compressed_size"],
                "uncompressed_size": member["uncompressed_size"],
            }
            for member in processed
        ],
        "selected_compressed_bytes": sum(
            member["compressed_size"] for member in processed
        ),
        "full_processed_archive_downloaded": False,
        "single_cell_data_downloaded": False,
        "original_preregistration_authorized": protocol["current_authorization"][
            "selective_member_acquisition"
        ],
        "execution_amendment_authorized": authorization["authorized_actions"][
            "selective_member_acquisition"
        ],
    }
    if dry_run:
        return plan
    if not confirm_selective_download:
        raise AcquisitionError(
            "selective acquisition requires --confirm-selective-download"
        )
    if not authorization["authorized_actions"]["selective_member_acquisition"]:
        raise AcquisitionError(
            "execution amendment sets selective_member_acquisition=false"
        )
    manifest_sha256 = sha256_file(manifest_path)
    protocol_sha256 = sha256_file(protocol_path)
    authorization_sha256 = sha256_file(authorization_path)
    recovery = validate_schema_recovery(
        load_json(DEFAULT_SCHEMA_RECOVERY, "acquisition schema recovery"),
        manifest_path=manifest_path,
        protocol_path=protocol_path,
    )
    with acquisition_lock():
        remote = audit_remote(manifest, protocol)
        return _publish_extraction(
            remote,
            manifest,
            protocol,
            extracted_root=extracted_root,
            audit_root=audit_root,
            manifest_sha256=manifest_sha256,
            protocol_sha256=protocol_sha256,
            authorization_sha256=authorization_sha256,
            authorization=authorization,
            recovery=recovery,
        )


def run_local(
    manifest_path: Path = DEFAULT_MANIFEST,
    protocol_path: Path = DEFAULT_PROTOCOL,
    authorization_path: Path = DEFAULT_AUTHORIZATION,
    extracted_root: Path = DEFAULT_EXTRACTED,
    audit_root: Path = DEFAULT_AUDIT,
) -> dict[str, Any]:
    manifest, protocol, authorization = load_execution_contracts(
        manifest_path,
        protocol_path,
        authorization_path,
    )
    manifest_sha256 = sha256_file(manifest_path)
    protocol_sha256 = sha256_file(protocol_path)
    authorization_sha256 = sha256_file(authorization_path)
    receipt = verify_receipt(
        audit_root / RECEIPT_NAME,
        extracted_root,
        manifest,
        protocol,
        manifest_sha256,
        protocol_sha256,
        authorization_sha256,
        authorization,
    )
    recovery = validate_schema_recovery(
        load_json(DEFAULT_SCHEMA_RECOVERY, "acquisition schema recovery"),
        manifest_path=manifest_path,
        protocol_path=protocol_path,
    )
    schema = validate_csv_schemas(extracted_root, manifest, recovery)
    return {
        "status": "pass",
        "receipt_created_at_utc": receipt.get("created_at_utc"),
        "schema": schema,
    }


def run_recover_receipt(
    *,
    confirm_receipt_recovery: bool,
    manifest_path: Path = DEFAULT_MANIFEST,
    protocol_path: Path = DEFAULT_PROTOCOL,
    authorization_path: Path = DEFAULT_AUTHORIZATION,
    extracted_root: Path = DEFAULT_EXTRACTED,
    audit_root: Path = DEFAULT_AUDIT,
) -> dict[str, Any]:
    if not confirm_receipt_recovery:
        raise AcquisitionError("receipt recovery requires --confirm-receipt-recovery")
    manifest, protocol, authorization = load_execution_contracts(
        manifest_path,
        protocol_path,
        authorization_path,
    )
    extracted_root = require_exact_source_root(extracted_root)
    audit_root = require_under_experiment(audit_root)
    receipt_path = audit_root / RECEIPT_NAME
    if receipt_path.exists() or receipt_path.is_symlink():
        raise AcquisitionError(
            f"refusing to overwrite existing receipt: {receipt_path}"
        )
    recovery = validate_schema_recovery(
        load_json(DEFAULT_SCHEMA_RECOVERY, "acquisition schema recovery"),
        manifest_path=manifest_path,
        protocol_path=protocol_path,
    )
    with acquisition_lock():
        schema = validate_csv_schemas(extracted_root, manifest, recovery)
        remote = audit_remote(manifest, protocol)
        receipt = build_receipt(
            extracted_root,
            manifest,
            protocol,
            remote.report,
            remote.processed_members,
            sha256_file(manifest_path),
            sha256_file(protocol_path),
            sha256_file(authorization_path),
            authorization,
        )
        atomic_write_json(receipt_path, receipt)
        verify_receipt(
            receipt_path,
            extracted_root,
            manifest,
            protocol,
            sha256_file(manifest_path),
            sha256_file(protocol_path),
            sha256_file(authorization_path),
            authorization,
        )
    return {
        "status": "recovered",
        "receipt": str(receipt_path),
        "selected_payload_redownloaded": False,
        "schema": schema,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate Experiment 111 provenance or selectively range-extract its "
            "six authorized files."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("manifest", help="validate only local frozen contracts")
    subparsers.add_parser(
        "remote",
        help="audit live metadata, small archives, and processed ZIP structure",
    )
    acquire_parser = subparsers.add_parser(
        "acquire",
        help="selectively range-extract the five processed members plus metadata",
    )
    acquire_parser.add_argument("--dry-run", action="store_true")
    acquire_parser.add_argument("--confirm-selective-download", action="store_true")
    subparsers.add_parser(
        "local",
        help="verify an existing receipt, extracted files, and CSV schemas",
    )
    recover_parser = subparsers.add_parser(
        "recover-receipt",
        help="validate already-published exact files and create a missing receipt",
    )
    recover_parser.add_argument("--confirm-receipt-recovery", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        if arguments.command == "manifest":
            result = run_manifest()
        elif arguments.command == "remote":
            result = run_remote()
        elif arguments.command == "acquire":
            result = run_acquire(
                confirm_selective_download=arguments.confirm_selective_download,
                dry_run=arguments.dry_run,
            )
        elif arguments.command == "local":
            result = run_local()
        elif arguments.command == "recover-receipt":
            result = run_recover_receipt(
                confirm_receipt_recovery=arguments.confirm_receipt_recovery,
            )
        else:
            raise AcquisitionError(f"unsupported command: {arguments.command}")
    except AcquisitionError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
