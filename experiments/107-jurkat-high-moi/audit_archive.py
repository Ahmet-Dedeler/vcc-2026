#!/usr/bin/env python3
"""Provenance, inventory, and sparse schema audit for GSE249595."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import math
import os
import re
import shutil
import sys
import tarfile
import tempfile
from array import array
from collections import Counter
from pathlib import Path, PurePosixPath
from typing import Any, BinaryIO

from acquire import (
    DEFAULT_DOWNLOADS,
    DEFAULT_MANIFEST,
    HERE,
    AcquisitionError,
    asset_binding,
    load_manifest,
    load_receipt,
    receipt_entry_matches_binding,
    require_under_experiment,
    sha256_file,
    verify_local_asset,
)


DEFAULT_PROVENANCE = HERE / "PROVENANCE.json"
DEFAULT_REPORT_ROOT = HERE / "audit"
MEMBER_PATTERN = re.compile(
    r"^(?P<accession>GSM\d+)_channel(?P<channel>\d+)_"
    r"(?P<modality>transcriptome|guides|labels)_"
    r"(?P<role>barcode|features|matrix)"
    r"(?P<suffix>\.tsv\.gz|\.mtx\.gz)$"
)


class AuditError(RuntimeError):
    """The archive or provenance violates the frozen acquisition contract."""


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


def nested_value(payload: dict[str, Any], dotted_key: str) -> Any:
    value: Any = payload
    for component in dotted_key.split("."):
        if not isinstance(value, dict) or component not in value:
            raise AuditError(f"missing provenance field: {dotted_key}")
        value = value[component]
    return value


def verify_provenance(
    manifest_path: Path = DEFAULT_MANIFEST,
    provenance_path: Path = DEFAULT_PROVENANCE,
    download_dir: Path = DEFAULT_DOWNLOADS,
    *,
    require_metadata: bool = False,
    require_raw: bool = False,
    verify_raw_hash: bool = False,
) -> dict[str, Any]:
    if verify_raw_hash and not require_raw:
        raise AuditError("--verify-raw-hash requires --require-raw")
    manifest_path = require_under_experiment(manifest_path)
    provenance_path = require_under_experiment(provenance_path)
    download_dir = require_under_experiment(download_dir)
    try:
        manifest = load_manifest(manifest_path)
    except AcquisitionError as error:
        raise AuditError(str(error)) from error
    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AuditError(f"cannot read provenance: {error}") from error

    verified_fields: dict[str, Any] = {}
    for dotted_key, expected in manifest["provenance_expectations"].items():
        actual = nested_value(provenance, dotted_key)
        if actual != expected:
            raise AuditError(
                f"provenance mismatch at {dotted_key}: {actual!r} != {expected!r}"
            )
        verified_fields[dotted_key] = actual

    source = provenance.get("source", {})
    url_fields = {
        "raw_archive": "raw_archive_url",
        "guide_reference": "guide_reference_url",
        "file_inventory": "file_inventory_url",
        "supplementary_table": "supplementary_table_url",
    }
    for asset_name, source_key in url_fields.items():
        expected_url = manifest["assets"][asset_name]["url"]
        if source.get(source_key) != expected_url:
            raise AuditError(f"provenance URL mismatch for {asset_name}")

    recorded_metadata = source.get("downloaded_metadata")
    if not isinstance(recorded_metadata, dict):
        raise AuditError("provenance has no downloaded_metadata map")
    for asset in manifest["assets"].values():
        if asset["group"] != "metadata":
            continue
        filename = str(asset["filename"])
        recorded = recorded_metadata.get(filename)
        if not isinstance(recorded, dict):
            raise AuditError(f"provenance has no metadata record for {filename}")
        if recorded.get("size_bytes") != asset["size_bytes"]:
            raise AuditError(f"provenance byte-count mismatch for {filename}")
        if recorded.get("sha256") != asset["sha256"]:
            raise AuditError(f"provenance SHA-256 mismatch for {filename}")

    downloaded: dict[str, Any] = {}
    if require_metadata or require_raw:
        try:
            receipt = load_receipt(download_dir)
        except AcquisitionError as error:
            raise AuditError(str(error)) from error
        expected_manifest_hash = sha256_file(manifest_path)
        if receipt.get("dataset") != manifest["dataset"]:
            raise AuditError("download receipt dataset does not match the manifest")
        if receipt.get("manifest_sha256") != expected_manifest_hash:
            raise AuditError("download receipt is not bound to the current manifest")

        groups = {"metadata"} if require_metadata else set()
        if require_raw:
            groups.add("raw")
        for name, asset in manifest["assets"].items():
            if asset["group"] not in groups:
                continue
            path = download_dir / asset["filename"]
            receipt_entry = receipt.get("assets", {}).get(name)
            if not isinstance(receipt_entry, dict):
                raise AuditError(f"download receipt has no entry for {name}")
            if not receipt_entry_matches_binding(
                receipt_entry,
                name,
                asset,
                expected_manifest_hash,
            ):
                expected_binding = asset_binding(name, asset, expected_manifest_hash)
                raise AuditError(
                    f"download receipt binding mismatch for {name}: "
                    f"expected {expected_binding}"
                )
            try:
                local = verify_local_asset(path, asset)
            except AcquisitionError as error:
                raise AuditError(str(error)) from error
            for key in ("filename", "size_bytes", "sha256"):
                if receipt_entry.get(key) != local[key]:
                    raise AuditError(f"download receipt mismatch for {name}.{key}")
            local["sha256_recomputed_now"] = True
            downloaded[name] = local

    return {
        "status": "pass",
        "manifest_sha256": sha256_file(manifest_path),
        "provenance_sha256": sha256_file(provenance_path),
        "verified_fields": verified_fields,
        "downloaded_assets": downloaded,
        "raw_hash_interpretation": (
            "The raw archive has no source-pinned SHA-256. Its receipt hash is a "
            "locally observed transport hash, not an independent source checksum."
        ),
    }


def parse_file_inventory(
    path: Path,
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise AuditError(f"cannot read file inventory: {error}") from error
    if not lines or lines[0].split("\t") != [
        "#Archive/File",
        "Name",
        "Time",
        "Size",
        "Type",
    ]:
        raise AuditError("unexpected file inventory header")

    archive_record: dict[str, Any] | None = None
    members: dict[str, dict[str, Any]] = {}
    for line_number, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != 5:
            raise AuditError(
                f"file inventory line {line_number} has {len(fields)} fields"
            )
        kind, name, timestamp, raw_size, file_type = fields
        try:
            size = int(raw_size)
        except ValueError as error:
            raise AuditError(
                f"file inventory line {line_number} has invalid size"
            ) from error
        record = {
            "name": name,
            "timestamp": timestamp,
            "size_bytes": size,
            "type": file_type,
        }
        if kind == "Archive":
            if archive_record is not None:
                raise AuditError("file inventory contains multiple archive rows")
            archive_record = record
        elif kind == "File":
            if name in members:
                raise AuditError(f"duplicate file inventory member: {name}")
            members[name] = record
        else:
            raise AuditError(f"unknown file inventory record type: {kind}")

    if archive_record is None:
        raise AuditError("file inventory has no archive row")
    return archive_record, members


def normalized_member_name(name: str) -> str:
    while name.startswith("./"):
        name = name[2:]
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or len(path.parts) != 1:
        raise AuditError(f"unsafe or nested archive member path: {name}")
    return path.as_posix()


def expected_member_grid(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    contract = manifest["archive_contract"]
    roles = contract["roles"]
    expected: dict[str, dict[str, Any]] = {}
    for channel_record in contract["channels"]:
        channel = int(channel_record["channel"])
        accession = str(channel_record["accession"])
        for modality in contract["modalities"]:
            for role, suffix in roles.items():
                filename = f"{accession}_channel{channel}_{modality}_{role}{suffix}"
                expected[filename] = {
                    "channel": channel,
                    "accession": accession,
                    "modality": modality,
                    "role": role,
                }
    return expected


def inventory_archive(
    archive_path: Path,
    filelist_path: Path,
    manifest: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, tarfile.TarInfo]]:
    archive_record, official_members = parse_file_inventory(filelist_path)
    raw_asset = manifest["assets"]["raw_archive"]
    if archive_record["name"] != raw_asset["filename"]:
        raise AuditError("official inventory names a different raw archive")
    if archive_record["size_bytes"] != raw_asset["size_bytes"]:
        raise AuditError("official inventory raw archive size differs from manifest")
    if archive_path.stat().st_size != raw_asset["size_bytes"]:
        raise AuditError("local raw archive byte count differs from manifest")

    expected_grid = expected_member_grid(manifest)
    if set(official_members) != set(expected_grid):
        missing = sorted(set(expected_grid) - set(official_members))
        extra = sorted(set(official_members) - set(expected_grid))
        raise AuditError(
            f"official member grid mismatch; missing={missing[:3]}, extra={extra[:3]}"
        )
    expected_count = int(manifest["archive_contract"]["expected_member_count"])
    if len(official_members) != expected_count:
        raise AuditError(
            f"official inventory has {len(official_members)} members, "
            f"expected {expected_count}"
        )

    tar_members: dict[str, tarfile.TarInfo] = {}
    try:
        with tarfile.open(archive_path, mode="r:") as archive:
            for member in archive:
                name = normalized_member_name(member.name)
                if not member.isfile():
                    raise AuditError(f"non-regular archive member: {name}")
                if name in tar_members:
                    raise AuditError(f"duplicate archive member: {name}")
                tar_members[name] = member
    except (OSError, tarfile.TarError) as error:
        raise AuditError(f"cannot inventory raw archive: {error}") from error

    if set(tar_members) != set(official_members):
        missing = sorted(set(official_members) - set(tar_members))
        extra = sorted(set(tar_members) - set(official_members))
        raise AuditError(
            f"archive member mismatch; missing={missing[:3]}, extra={extra[:3]}"
        )
    for name, member in tar_members.items():
        expected_size = official_members[name]["size_bytes"]
        if member.size != expected_size:
            raise AuditError(
                f"{name}: tar member size {member.size} != official {expected_size}"
            )
        match = MEMBER_PATTERN.fullmatch(name)
        if match is None:
            raise AuditError(f"unexpected member filename grammar: {name}")
        grid = expected_grid[name]
        for key in ("accession", "modality", "role"):
            if match.group(key) != str(grid[key]):
                raise AuditError(f"{name}: filename {key} does not match manifest")
        if int(match.group("channel")) != int(grid["channel"]):
            raise AuditError(f"{name}: filename channel does not match manifest")

    report = {
        "status": "pass",
        "archive": {
            "filename": archive_path.name,
            "size_bytes": archive_path.stat().st_size,
        },
        "file_inventory_sha256": sha256_file(filelist_path),
        "member_count": len(tar_members),
        "member_total_compressed_bytes": sum(
            member.size for member in tar_members.values()
        ),
        "channels": [
            {
                "channel": record["channel"],
                "accession": record["accession"],
                "member_count": 9,
            }
            for record in manifest["archive_contract"]["channels"]
        ],
    }
    return report, tar_members


def extracted_member(
    archive: tarfile.TarFile,
    member: tarfile.TarInfo,
) -> BinaryIO:
    handle = archive.extractfile(member)
    if handle is None:
        raise AuditError(f"cannot read archive member: {member.name}")
    return handle


def gzip_text_stats(
    archive: tarfile.TarFile,
    member: tarfile.TarInfo,
    *,
    capture_rows: bool = False,
    sample_limit: int = 5,
    expected_fields: int | None = None,
) -> tuple[dict[str, Any], list[list[str]]]:
    digest = hashlib.sha256()
    line_count = 0
    sample: list[list[str]] = []
    rows: list[list[str]] = []
    try:
        with extracted_member(archive, member) as compressed:
            with gzip.GzipFile(fileobj=compressed, mode="rb") as uncompressed:
                for raw_line in uncompressed:
                    digest.update(raw_line)
                    line_count += 1
                    try:
                        decoded = raw_line.rstrip(b"\r\n").decode("utf-8")
                    except UnicodeDecodeError as error:
                        raise AuditError(
                            f"{member.name}: non-UTF-8 TSV line {line_count}"
                        ) from error
                    fields = decoded.split("\t")
                    if not fields or not fields[0]:
                        raise AuditError(
                            f"{member.name}: empty first field at line {line_count}"
                        )
                    if expected_fields is not None and len(fields) != expected_fields:
                        raise AuditError(
                            f"{member.name}: expected {expected_fields} fields at "
                            f"line {line_count}, found {len(fields)}"
                        )
                    if len(sample) < sample_limit:
                        sample.append(fields)
                    if capture_rows:
                        rows.append(fields)
    except (OSError, EOFError) as error:
        raise AuditError(f"{member.name}: invalid gzip stream: {error}") from error
    return (
        {
            "rows": line_count,
            "uncompressed_sha256": digest.hexdigest(),
            "sample": sample,
            "gzip_crc_verified": True,
        },
        rows,
    )


def read_matrix_header(
    archive: tarfile.TarFile,
    member: tarfile.TarInfo,
) -> dict[str, Any]:
    try:
        with extracted_member(archive, member) as compressed:
            with gzip.GzipFile(fileobj=compressed, mode="rb") as matrix:
                banner = matrix.readline().decode("ascii").strip().split()
                if len(banner) != 5 or banner[:3] != [
                    "%%MatrixMarket",
                    "matrix",
                    "coordinate",
                ]:
                    raise AuditError(f"{member.name}: unsupported Matrix Market banner")
                field = banner[3].lower()
                symmetry = banner[4].lower()
                if field not in {"integer", "real", "pattern"}:
                    raise AuditError(f"{member.name}: unsupported field {field}")
                if symmetry != "general":
                    raise AuditError(f"{member.name}: unsupported symmetry {symmetry}")
                while True:
                    line = matrix.readline()
                    if not line:
                        raise AuditError(f"{member.name}: missing dimension line")
                    if not line.startswith(b"%"):
                        break
                dimensions = line.split()
                if len(dimensions) != 3:
                    raise AuditError(f"{member.name}: invalid dimension line")
                rows, columns, nnz = (int(value) for value in dimensions)
    except (UnicodeDecodeError, ValueError, OSError, EOFError) as error:
        if isinstance(error, AuditError):
            raise
        raise AuditError(
            f"{member.name}: invalid Matrix Market header: {error}"
        ) from error
    if rows <= 0 or columns <= 0 or nnz < 0 or nnz > rows * columns:
        raise AuditError(f"{member.name}: impossible Matrix Market dimensions")
    return {
        "rows": rows,
        "columns": columns,
        "nnz": nnz,
        "field": field,
        "symmetry": symmetry,
        "gzip_crc_verified": False,
    }


def hash_compressed_member(
    archive: tarfile.TarFile,
    member: tarfile.TarInfo,
) -> str:
    digest = hashlib.sha256()
    with extracted_member(archive, member) as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()


def histogram_quantile(histogram: Counter[int], quantile: float) -> int:
    total = sum(histogram.values())
    if total == 0:
        return 0
    wanted = max(1, math.ceil(quantile * total))
    cumulative = 0
    for value, count in sorted(histogram.items()):
        cumulative += count
        if cumulative >= wanted:
            return value
    raise AssertionError("quantile traversal did not terminate")


def parse_matrix_prefix(
    matrix: BinaryIO, member_name: str
) -> tuple[str, str, int, int, int]:
    banner = matrix.readline().decode("ascii").strip().split()
    if len(banner) != 5 or banner[:3] != [
        "%%MatrixMarket",
        "matrix",
        "coordinate",
    ]:
        raise AuditError(f"{member_name}: unsupported Matrix Market banner")
    field = banner[3].lower()
    symmetry = banner[4].lower()
    if field not in {"integer", "real", "pattern"}:
        raise AuditError(f"{member_name}: unsupported field {field}")
    if symmetry != "general":
        raise AuditError(f"{member_name}: unsupported symmetry {symmetry}")
    while True:
        line = matrix.readline()
        if not line:
            raise AuditError(f"{member_name}: missing dimension line")
        if not line.startswith(b"%"):
            break
    parts = line.split()
    if len(parts) != 3:
        raise AuditError(f"{member_name}: invalid dimension line")
    rows, columns, nnz = (int(value) for value in parts)
    if rows <= 0 or columns <= 0 or nnz < 0 or nnz > rows * columns:
        raise AuditError(f"{member_name}: impossible Matrix Market dimensions")
    return field, symmetry, rows, columns, nnz


def scan_assignment_matrix(
    archive: tarfile.TarFile,
    member: tarfile.TarInfo,
    *,
    expected_rows: int,
    expected_columns: int,
    max_nnz: int,
) -> tuple[dict[str, Any], array]:
    try:
        with extracted_member(archive, member) as compressed:
            with gzip.GzipFile(fileobj=compressed, mode="rb") as matrix:
                field, symmetry, rows, columns, declared_nnz = parse_matrix_prefix(
                    matrix, member.name
                )
                if field != "integer" or symmetry != "general":
                    raise AuditError(
                        f"{member.name}: assignment matrices must be integer general"
                    )
                if rows != expected_rows or columns != expected_columns:
                    raise AuditError(
                        f"{member.name}: assignment dimensions changed between "
                        "header validation and scan"
                    )
                if declared_nnz > max_nnz:
                    raise AuditError(
                        f"{member.name}: declared nnz {declared_nnz} exceeds "
                        f"remaining safety budget {max_nnz}"
                    )
                row_support = array("Q", [0]) * rows
                column_support = array("I", [0]) * columns
                observed_nnz = 0
                previous_coordinate: tuple[int, int] | None = None
                row_major_possible = True
                column_major_possible = True
                for line_number, line in enumerate(matrix, start=1):
                    if observed_nnz >= declared_nnz:
                        raise AuditError(
                            f"{member.name}: more entries than declared at "
                            f"line {line_number}"
                        )
                    parts = line.split()
                    if len(parts) != 3:
                        raise AuditError(
                            f"{member.name}: malformed coordinate line {line_number}"
                        )
                    row = int(parts[0])
                    column = int(parts[1])
                    if not 1 <= row <= rows or not 1 <= column <= columns:
                        raise AuditError(
                            f"{member.name}: out-of-bounds coordinate at "
                            f"line {line_number}"
                        )
                    value = int(parts[2])
                    if value <= 0:
                        raise AuditError(
                            f"{member.name}: assignment value must be positive "
                            f"at line {line_number}"
                        )
                    coordinate = (row, column)
                    if previous_coordinate is not None:
                        if coordinate == previous_coordinate:
                            raise AuditError(
                                f"{member.name}: duplicate coordinate at "
                                f"line {line_number}"
                            )
                        row_major_possible = row_major_possible and (
                            coordinate > previous_coordinate
                        )
                        column_major_possible = column_major_possible and (
                            (column, row)
                            > (previous_coordinate[1], previous_coordinate[0])
                        )
                        if not row_major_possible and not column_major_possible:
                            raise AuditError(
                                f"{member.name}: coordinates are not in a canonical "
                                "row-major or column-major order"
                            )
                    previous_coordinate = coordinate
                    row_support[row - 1] += 1
                    column_support[column - 1] += 1
                    observed_nnz += 1
    except (UnicodeDecodeError, ValueError, OverflowError, OSError, EOFError) as error:
        if isinstance(error, AuditError):
            raise
        raise AuditError(
            f"{member.name}: invalid coordinate stream: {error}"
        ) from error

    if observed_nnz != declared_nnz:
        raise AuditError(
            f"{member.name}: observed {observed_nnz} entries, declared {declared_nnz}"
        )
    histogram = Counter(column_support)
    active_histogram = Counter(
        {support: count for support, count in histogram.items() if support > 0}
    )
    row_bytes = row_support.tobytes()
    return (
        {
            "declared_nnz": declared_nnz,
            "observed_nnz": observed_nnz,
            "coordinate_order": ("row_major" if row_major_possible else "column_major"),
            "supported_features": sum(value > 0 for value in row_support),
            "active_cells": columns - histogram.get(0, 0),
            "features_per_active_cell": {
                "min": min(active_histogram, default=0),
                "median": histogram_quantile(active_histogram, 0.5),
                "p95": histogram_quantile(active_histogram, 0.95),
                "max": max(active_histogram, default=0),
                "histogram": {
                    str(key): value for key, value in sorted(active_histogram.items())
                },
            },
            "row_support_sha256": hashlib.sha256(row_bytes).hexdigest(),
            "gzip_crc_verified": True,
        },
        row_support,
    )


def write_support_table(
    path: Path,
    feature_rows: list[list[str]],
    row_support: array,
) -> dict[str, Any]:
    if len(feature_rows) != len(row_support):
        raise AuditError(f"{path.name}: feature/support row mismatch")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(descriptor, "wb") as raw:
            with gzip.GzipFile(
                filename="",
                fileobj=raw,
                mode="wb",
                mtime=0,
            ) as compressed:
                with io.TextIOWrapper(
                    compressed,
                    encoding="utf-8",
                    newline="",
                ) as handle:
                    handle.write("row_index\tfeature_id\tfeature_name\tnnz_cells\n")
                    for index, (fields, count) in enumerate(
                        zip(feature_rows, row_support, strict=True), start=1
                    ):
                        feature_id = fields[0].replace("\t", " ")
                        feature_name = (
                            fields[1] if len(fields) > 1 else fields[0]
                        ).replace("\t", " ")
                        handle.write(
                            f"{index}\t{feature_id}\t{feature_name}\t{count}\n"
                        )
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise
    return {
        "path": str(path.relative_to(HERE)),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def create_staging_directory(final_directory: Path) -> Path:
    final_directory = require_under_experiment(final_directory)
    final_directory.parent.mkdir(parents=True, exist_ok=True)
    if final_directory.is_symlink():
        raise AuditError(f"support directory cannot be a symlink: {final_directory}")
    if final_directory.exists() and not final_directory.is_dir():
        raise AuditError(f"support path is not a directory: {final_directory}")
    return Path(
        tempfile.mkdtemp(
            dir=final_directory.parent,
            prefix=f".{final_directory.name}.stage-",
        )
    )


def publish_staged_directory(staged: Path, destination: Path) -> None:
    destination = require_under_experiment(destination)
    backup: Path | None = None
    if destination.exists():
        backup = Path(
            tempfile.mkdtemp(
                dir=destination.parent,
                prefix=f".{destination.name}.previous-",
            )
        )
        backup.rmdir()
        os.replace(destination, backup)
    try:
        os.replace(staged, destination)
    except BaseException:
        if backup is not None and backup.exists() and not destination.exists():
            os.replace(backup, destination)
        raise
    if backup is not None:
        shutil.rmtree(backup)


def schema_audit(
    archive_path: Path,
    filelist_path: Path,
    manifest: dict[str, Any],
    *,
    verify_barcode_hashes: bool = False,
    scan_assignments: bool = False,
    support_dir: Path | None = None,
    input_binding: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if scan_assignments and input_binding is None:
        input_binding = {
            "manifest_contract_sha256": hashlib.sha256(
                json.dumps(
                    manifest,
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode("utf-8")
            ).hexdigest(),
            "raw_archive_sha256": sha256_file(archive_path),
            "file_inventory_sha256": sha256_file(filelist_path),
            "raw_hash_recomputed_for_this_run": True,
        }
    inventory, _ = inventory_archive(archive_path, filelist_path, manifest)
    grid = expected_member_grid(manifest)
    feature_hashes: dict[str, set[str]] = {
        modality: set() for modality in manifest["archive_contract"]["modalities"]
    }
    channel_reports: list[dict[str, Any]] = []
    declared_assignment_nnz = 0

    try:
        with tarfile.open(archive_path, mode="r:") as archive:
            members = {
                normalized_member_name(member.name): member for member in archive
            }
            for channel_record in manifest["archive_contract"]["channels"]:
                channel = int(channel_record["channel"])
                accession = str(channel_record["accession"])
                channel_report: dict[str, Any] = {
                    "channel": channel,
                    "accession": accession,
                    "modalities": {},
                }
                barcode_hashes: dict[str, str] = {}
                barcode_stats_by_modality: dict[str, dict[str, Any]] = {}
                for modality, dimensions in manifest["archive_contract"][
                    "modalities"
                ].items():
                    prefix = f"{accession}_channel{channel}_{modality}_"
                    matrix_name = f"{prefix}matrix.mtx.gz"
                    features_name = f"{prefix}features.tsv.gz"
                    barcode_name = f"{prefix}barcode.tsv.gz"
                    for name in (matrix_name, features_name, barcode_name):
                        if name not in grid or name not in members:
                            raise AuditError(f"missing schema member: {name}")

                    header = read_matrix_header(archive, members[matrix_name])
                    if header["rows"] != int(dimensions["feature_rows"]):
                        raise AuditError(f"{matrix_name}: unexpected feature dimension")
                    if header["columns"] != int(dimensions["barcode_columns"]):
                        raise AuditError(f"{matrix_name}: unexpected barcode dimension")

                    feature_stats, _ = gzip_text_stats(archive, members[features_name])
                    if feature_stats["rows"] != header["rows"]:
                        raise AuditError(
                            f"{features_name}: row count does not match matrix"
                        )
                    feature_hashes[modality].add(feature_stats["uncompressed_sha256"])
                    modality_report: dict[str, Any] = {
                        "matrix": header,
                        "features": feature_stats,
                    }
                    if verify_barcode_hashes:
                        barcode_stats, _ = gzip_text_stats(
                            archive,
                            members[barcode_name],
                            sample_limit=5 if modality == "transcriptome" else 0,
                            expected_fields=1,
                        )
                        expected_barcodes = int(dimensions["barcode_columns"])
                        if barcode_stats["rows"] != expected_barcodes:
                            raise AuditError(
                                f"{barcode_name}: barcode row count "
                                f"{barcode_stats['rows']} != expected "
                                f"{expected_barcodes}"
                            )
                        barcode_hash = barcode_stats["uncompressed_sha256"]
                        barcode_hashes[modality] = barcode_hash
                        barcode_stats_by_modality[modality] = barcode_stats
                        modality_report["barcode"] = {
                            **barcode_stats,
                            "exact_axis_check": "uncompressed-byte identity",
                        }
                    if modality in {"guides", "labels"}:
                        declared_assignment_nnz += int(header["nnz"])
                    channel_report["modalities"][modality] = modality_report

                if verify_barcode_hashes and len(set(barcode_hashes.values())) != 1:
                    raise AuditError(
                        f"channel {channel}: barcode files are not byte-identical"
                    )
                if verify_barcode_hashes:
                    barcode_stats = barcode_stats_by_modality["transcriptome"]
                    expected_barcodes = int(
                        manifest["archive_contract"]["modalities"]["transcriptome"][
                            "barcode_columns"
                        ]
                    )
                    if barcode_stats["rows"] != expected_barcodes:
                        raise AuditError(
                            f"channel {channel}: barcode row count "
                            f"{barcode_stats['rows']} != expected {expected_barcodes}"
                        )
                    channel_report["barcode_axis_sha256"] = next(
                        iter(barcode_hashes.values())
                    )
                    channel_report["barcode_axis"] = {
                        **barcode_stats,
                        "expected_rows": expected_barcodes,
                        "uncompressed_byte_identity_across_modalities": True,
                    }
                channel_reports.append(channel_report)

            if any(len(values) != 1 for values in feature_hashes.values()):
                raise AuditError("feature axes differ between channels")

            max_assignment_nnz = int(manifest["safety"]["max_assignment_nnz_total"])
            if declared_assignment_nnz > max_assignment_nnz:
                raise AuditError(
                    f"declared assignment nnz {declared_assignment_nnz} exceeds "
                    f"safety ceiling {max_assignment_nnz}"
                )

            if scan_assignments:
                if support_dir is None:
                    raise AuditError(
                        "assignment scanning requires an explicit support directory"
                    )
                support_dir = require_under_experiment(support_dir)
                staged_support_dir = create_staging_directory(support_dir)
                support_artifacts: list[dict[str, Any]] = []
                observed_assignment_nnz = 0
                try:
                    for channel_report in channel_reports:
                        channel = int(channel_report["channel"])
                        accession = str(channel_report["accession"])
                        for modality in ("guides", "labels"):
                            prefix = f"{accession}_channel{channel}_{modality}_"
                            features_name = f"{prefix}features.tsv.gz"
                            matrix_name = f"{prefix}matrix.mtx.gz"
                            _, feature_rows = gzip_text_stats(
                                archive,
                                members[features_name],
                                capture_rows=True,
                                sample_limit=0,
                            )
                            dimensions = manifest["archive_contract"]["modalities"][
                                modality
                            ]
                            remaining_nnz = max_assignment_nnz - observed_assignment_nnz
                            scan_report, row_support = scan_assignment_matrix(
                                archive,
                                members[matrix_name],
                                expected_rows=int(dimensions["feature_rows"]),
                                expected_columns=int(dimensions["barcode_columns"]),
                                max_nnz=remaining_nnz,
                            )
                            observed_assignment_nnz += int(scan_report["observed_nnz"])
                            support_filename = (
                                f"channel{channel:02d}_{modality}_row_support.tsv.gz"
                            )
                            staged_support_path = staged_support_dir / support_filename
                            support_record = write_support_table(
                                staged_support_path,
                                feature_rows,
                                row_support,
                            )
                            support_record["path"] = str(
                                (support_dir / support_filename).relative_to(HERE)
                            )
                            support_artifacts.append(support_record)
                            scan_report["support_table"] = support_record
                            channel_report["modalities"][modality][
                                "assignment_scan"
                            ] = scan_report
                    expected_support_tables = (
                        len(manifest["archive_contract"]["channels"]) * 2
                    )
                    if len(support_artifacts) != expected_support_tables:
                        raise AuditError(
                            f"built {len(support_artifacts)} support tables, "
                            f"expected {expected_support_tables}"
                        )
                    atomic_write_json(
                        staged_support_dir / "COMPLETE.json",
                        {
                            "status": "complete",
                            "expected_support_tables": expected_support_tables,
                            "support_tables": support_artifacts,
                            "input_binding": input_binding,
                        },
                    )
                    publish_staged_directory(staged_support_dir, support_dir)
                except BaseException:
                    if staged_support_dir.exists():
                        shutil.rmtree(staged_support_dir)
                    raise
    except (OSError, tarfile.TarError) as error:
        if isinstance(error, AuditError):
            raise
        raise AuditError(f"schema audit failed: {error}") from error

    report = {
        "status": "pass",
        "inventory": inventory,
        "feature_axis_sha256": {
            modality: next(iter(values)) for modality, values in feature_hashes.items()
        },
        "barcode_axes_verified": verify_barcode_hashes,
        "assignment_coordinates_scanned": scan_assignments,
        "transcriptome_coordinates_scanned": False,
        "declared_assignment_nnz_total": declared_assignment_nnz,
        "channels": channel_reports,
    }
    if input_binding is not None:
        report["input_binding"] = input_binding
    if scan_assignments and support_dir is not None:
        report["support_completion_manifest"] = str(
            (support_dir / "COMPLETE.json").relative_to(HERE)
        )
    return report


def verify_all_gzip_members(
    archive_path: Path,
    filelist_path: Path,
    manifest: dict[str, Any],
) -> dict[str, Any]:
    inventory, _ = inventory_archive(archive_path, filelist_path, manifest)
    reports: list[dict[str, Any]] = []
    try:
        with tarfile.open(archive_path, mode="r:") as archive:
            for member in archive:
                name = normalized_member_name(member.name)
                digest = hashlib.sha256()
                uncompressed_bytes = 0
                with extracted_member(archive, member) as compressed:
                    with gzip.GzipFile(fileobj=compressed, mode="rb") as stream:
                        while block := stream.read(8 << 20):
                            digest.update(block)
                            uncompressed_bytes += len(block)
                reports.append(
                    {
                        "name": name,
                        "uncompressed_bytes": uncompressed_bytes,
                        "uncompressed_sha256": digest.hexdigest(),
                    }
                )
    except (OSError, EOFError, tarfile.TarError) as error:
        raise AuditError(f"gzip integrity audit failed: {error}") from error
    return {
        "status": "pass",
        "inventory": inventory,
        "member_count": len(reports),
        "all_members_read_to_gzip_eof": True,
        "members": reports,
    }


def emit_report(payload: dict[str, Any], report_path: Path) -> None:
    report_path = require_under_experiment(report_path)
    atomic_write_json(report_path, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


def invalidate_report(path: Path) -> Path:
    path = require_under_experiment(path)
    if path.is_symlink():
        raise AuditError(f"report path cannot be a symlink: {path}")
    if path.exists() and not path.is_file():
        raise AuditError(f"report path is not a regular file: {path}")
    path.unlink(missing_ok=True)
    return path


def build_input_binding(provenance_report: dict[str, Any]) -> dict[str, Any]:
    downloaded = provenance_report["downloaded_assets"]
    return {
        "manifest_sha256": provenance_report["manifest_sha256"],
        "provenance_sha256": provenance_report["provenance_sha256"],
        "raw_archive_sha256": downloaded["raw_archive"]["sha256"],
        "file_inventory_sha256": downloaded["file_inventory"]["sha256"],
        "raw_hash_recomputed_for_this_run": True,
    }


def common_paths(
    args: argparse.Namespace,
) -> tuple[dict[str, Any], Path, Path]:
    manifest_path = require_under_experiment(Path(args.manifest))
    download_dir = require_under_experiment(Path(args.download_dir))
    manifest = load_manifest(manifest_path)
    archive_path = download_dir / manifest["assets"]["raw_archive"]["filename"]
    filelist_path = download_dir / manifest["assets"]["file_inventory"]["filename"]
    return manifest, archive_path, filelist_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--download-dir", default=str(DEFAULT_DOWNLOADS))
    subparsers = parser.add_subparsers(dest="command", required=True)

    provenance = subparsers.add_parser("provenance")
    provenance.add_argument("--require-metadata", action="store_true")
    provenance.add_argument("--require-raw", action="store_true")
    provenance.add_argument("--verify-raw-hash", action="store_true")
    provenance.add_argument(
        "--report",
        default=str(DEFAULT_REPORT_ROOT / "provenance.json"),
    )

    inventory = subparsers.add_parser("inventory")
    inventory.add_argument(
        "--report",
        default=str(DEFAULT_REPORT_ROOT / "inventory.json"),
    )

    schema = subparsers.add_parser("schema")
    schema.add_argument("--verify-barcode-hashes", action="store_true")
    schema.add_argument("--scan-assignments", action="store_true")
    schema.add_argument(
        "--support-dir",
        default=str(DEFAULT_REPORT_ROOT / "support"),
    )
    schema.add_argument(
        "--report",
        default=str(DEFAULT_REPORT_ROOT / "schema.json"),
    )

    gzip_parser = subparsers.add_parser("verify-gzip")
    gzip_parser.add_argument("--confirm-full-stream", action="store_true")
    gzip_parser.add_argument(
        "--report",
        default=str(DEFAULT_REPORT_ROOT / "gzip-integrity.json"),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        report_path = invalidate_report(Path(args.report))
        if args.command == "provenance":
            manifest_path = require_under_experiment(Path(args.manifest))
            download_dir = require_under_experiment(Path(args.download_dir))
            payload = verify_provenance(
                manifest_path,
                DEFAULT_PROVENANCE,
                download_dir,
                require_metadata=args.require_metadata,
                require_raw=args.require_raw,
                verify_raw_hash=args.verify_raw_hash,
            )
        else:
            manifest, archive_path, filelist_path = common_paths(args)
            provenance_report = verify_provenance(
                Path(args.manifest),
                DEFAULT_PROVENANCE,
                Path(args.download_dir),
                require_metadata=True,
                require_raw=True,
                verify_raw_hash=True,
            )
            input_binding = build_input_binding(provenance_report)
            if args.command == "inventory":
                payload, _ = inventory_archive(archive_path, filelist_path, manifest)
            elif args.command == "schema":
                support_dir = require_under_experiment(Path(args.support_dir))
                payload = schema_audit(
                    archive_path,
                    filelist_path,
                    manifest,
                    verify_barcode_hashes=args.verify_barcode_hashes,
                    scan_assignments=args.scan_assignments,
                    support_dir=support_dir if args.scan_assignments else None,
                    input_binding=input_binding,
                )
            elif args.command == "verify-gzip":
                if not args.confirm_full_stream:
                    raise AuditError(
                        "full gzip validation requires --confirm-full-stream"
                    )
                payload = verify_all_gzip_members(archive_path, filelist_path, manifest)
            else:
                raise AssertionError(f"unhandled command: {args.command}")
            payload["input_binding"] = input_binding
        emit_report(payload, report_path)
    except (AcquisitionError, AuditError, FileNotFoundError) as error:
        print(f"audit failed: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
