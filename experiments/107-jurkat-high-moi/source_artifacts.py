#!/usr/bin/env python3
"""Fail-closed verification for Experiment 107 source-stage artifacts.

This module verifies provenance and derived cell-call artifacts without opening
or parsing transcriptome matrix members.  The raw archive is read only as an
opaque byte stream for its size and SHA-256.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import stat
from array import array
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any, BinaryIO, Iterator

import numpy as np
import scipy


HERE = Path(__file__).resolve().parent
SHA256_LENGTH = 64


class SourceArtifactError(RuntimeError):
    """A source artifact failed a frozen integrity or semantic check."""


@dataclass(frozen=True)
class SourceArtifactPaths:
    """All files consumed by the source-artifact verifier."""

    root: Path
    protocol: Path
    manifest: Path
    provenance: Path
    archive: Path
    file_inventory: Path
    schema_report: Path
    gzip_report: Path
    support_report: Path
    oracle_receipt: Path
    preflight_receipt: Path
    preflight_code: Path
    preprocess_code: Path
    r_oracle_code: Path
    source_stage_code: Path
    work: Path

    @classmethod
    def from_root(cls, root: Path) -> SourceArtifactPaths:
        """Build the canonical Experiment 107 path layout under ``root``."""

        root = Path(root)
        return cls(
            root=root,
            protocol=root / "preregistered_protocol.json",
            manifest=root / "ACQUISITION_MANIFEST.json",
            provenance=root / "PROVENANCE.json",
            archive=root / "downloads" / "GSE249595_RAW.tar",
            file_inventory=root / "downloads" / "filelist.txt",
            schema_report=root / "audit" / "schema.json",
            gzip_report=root / "audit" / "gzip-integrity.json",
            support_report=root / "audit" / "support" / "COMPLETE.json",
            oracle_receipt=root / "R_ORACLE_RECEIPT.json",
            preflight_receipt=root / "SOURCE_PREFLIGHT_RECEIPT.json",
            preflight_code=root / "source_preflight.py",
            preprocess_code=root / "preprocess.py",
            r_oracle_code=root / "r_oracle.R",
            source_stage_code=root / "source_stage.py",
            work=root / "work",
        )


@dataclass(frozen=True)
class VerificationContract:
    """Immutable values that cannot be supplied by mutable reports."""

    dataset: str
    protocol_sha256: str
    manifest_sha256: str
    provenance_sha256: str
    raw_archive_sha256: str
    raw_archive_size_bytes: int
    file_inventory_sha256: str
    channels: tuple[int, ...]
    accessions: tuple[str, ...]
    barcode_columns: int
    feature_rows: tuple[tuple[str, int], ...]
    production_iterations: int = 10_000
    lower: int = 100
    fdr_threshold: float = 0.01
    manifest_assignment_coordinate_order: str = "row_or_column_major_strict"
    coordinate_order: str = "column_major"
    required_oracle_checks: tuple[str, ...] = (
        "column_totals",
        "barcode_ranks_knee",
        "barcode_ranks_inflection",
        "safe_good_turing_ambient_profile",
        "multinomial_observed_data_log_probability",
        "benjamini_hochberg",
    )

    def rows_by_modality(self) -> dict[str, int]:
        """Return the frozen feature counts keyed by modality."""

        return dict(self.feature_rows)


PRODUCTION_CONTRACT = VerificationContract(
    dataset="GSE249595",
    protocol_sha256=(
        "68f6598277b73201885e94582fba0d15295b936a7184312e08356eb2604a3748"
    ),
    manifest_sha256=(
        "59e37105c31b68cbca84b3c77834d1fa53b6316b2b351c1d0029f51f0f837c74"
    ),
    provenance_sha256=(
        "0c3277a5ee1fa9b179d5f91683375c7f2a9b5a72dd38c2e00c35df4950a1c722"
    ),
    raw_archive_sha256=(
        "3c0c5fa8cab2983dfe0fa3b38b31662bb68e191e6dbaf3d03af02a67587a67a9"
    ),
    raw_archive_size_bytes=3_753_820_160,
    file_inventory_sha256=(
        "ce492b775742632b3b07c854c37712b62ef6bb3262050edb8c81950a96db8eca"
    ),
    channels=tuple(range(1, 17)),
    accessions=tuple(f"GSM79514{value:02d}" for value in range(13, 29)),
    barcode_columns=6_794_880,
    feature_rows=(
        ("transcriptome", 20_606),
        ("guides", 83_401),
        ("labels", 7),
    ),
)


PREPARE_ARRAY_KEYS = frozenset(
    {
        "totals",
        "detected_features",
        "mitochondrial_counts",
        "feature_totals",
        "ambient_counts",
        "ambient_probability",
        "tested_indices",
        "observed_data_log_probability",
    }
)
CALL_ARRAY_KEYS = frozenset(
    {
        "tested_indices",
        "p_value",
        "fdr",
        "limited",
        "accepted_indices",
    }
)


def _fail(message: str) -> None:
    raise SourceArtifactError(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _require_sha256(value: Any, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != SHA256_LENGTH
        or any(character not in "0123456789abcdef" for character in value)
    ):
        _fail(f"{label} is not a lowercase SHA-256")
    return value


def _require_int(value: Any, label: str, *, minimum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        _fail(f"{label} is not an integer")
    if minimum is not None and value < minimum:
        _fail(f"{label} is below {minimum}")
    return value


def _require_number(
    value: Any,
    label: str,
    *,
    minimum: float | None = None,
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        _fail(f"{label} is not numeric")
    number = float(value)
    if not math.isfinite(number):
        _fail(f"{label} is not finite")
    if minimum is not None and number < minimum:
        _fail(f"{label} is below {minimum}")
    return number


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON constant: {value}")


def _root(paths: SourceArtifactPaths) -> Path:
    try:
        resolved = paths.root.resolve(strict=True)
    except OSError as error:
        raise SourceArtifactError(
            f"cannot resolve experiment root {paths.root}: {error}"
        ) from error
    if paths.root.is_symlink() or not resolved.is_dir():
        _fail(f"experiment root is missing, unsafe, or not a directory: {paths.root}")
    return resolved


def _confined_path(
    paths: SourceArtifactPaths,
    path: Path,
    label: str,
    *,
    must_exist: bool = True,
) -> Path:
    root = _root(paths)
    try:
        resolved = Path(path).resolve(strict=must_exist)
    except OSError as error:
        raise SourceArtifactError(f"cannot resolve {label} {path}: {error}") from error
    try:
        resolved.relative_to(root)
    except ValueError:
        _fail(f"{label} escapes the experiment root: {path}")
    if must_exist and Path(path).is_symlink():
        _fail(f"{label} cannot be a symlink: {path}")
    return resolved


@contextmanager
def _open_regular_file(
    paths: SourceArtifactPaths,
    path: Path,
    label: str,
) -> Iterator[BinaryIO]:
    resolved = _confined_path(paths, path, label)
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(resolved, flags)
    except OSError as error:
        raise SourceArtifactError(f"cannot open {label} {resolved}: {error}") from error
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            _fail(f"{label} is not a regular file: {resolved}")
        with os.fdopen(descriptor, "rb", closefd=False) as handle:
            yield handle
    finally:
        os.close(descriptor)


def _hash_handle(handle: BinaryIO, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    handle.seek(0)
    for block in iter(lambda: handle.read(block_bytes), b""):
        digest.update(block)
    return digest.hexdigest()


def sha256_file(
    path: Path,
    *,
    paths: SourceArtifactPaths | None = None,
    label: str = "file",
) -> str:
    """Hash a regular file, rejecting symlinks and in-flight mutation."""

    if paths is None:
        root = Path(path).resolve(strict=True).parent
        paths = SourceArtifactPaths.from_root(root)
        paths = SourceArtifactPaths(
            **{
                **paths.__dict__,
                "root": root,
            }
        )
    with _open_regular_file(paths, path, label) as handle:
        before = os.fstat(handle.fileno())
        digest = _hash_handle(handle)
        after = os.fstat(handle.fileno())
    identity_before = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_ctime_ns,
    )
    identity_after = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    )
    if identity_before != identity_after:
        _fail(f"{label} changed while it was being hashed: {path}")
    return digest


def _load_json(
    paths: SourceArtifactPaths,
    path: Path,
    label: str,
) -> tuple[dict[str, Any], str]:
    with _open_regular_file(paths, path, label) as handle:
        payload = handle.read()
    digest = hashlib.sha256(payload).hexdigest()
    try:
        parsed = json.loads(
            payload,
            object_pairs_hook=_reject_duplicate_keys,
            parse_constant=_reject_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise SourceArtifactError(f"cannot parse {label} {path}: {error}") from error
    if not isinstance(parsed, dict):
        _fail(f"{label} is not a JSON object: {path}")
    return parsed, digest


def _relative_report_path(
    paths: SourceArtifactPaths,
    value: Any,
    label: str,
) -> Path:
    if not isinstance(value, str) or not value:
        _fail(f"{label} is not a non-empty relative path")
    pure = PurePosixPath(value)
    if pure.is_absolute() or ".." in pure.parts or "." in pure.parts:
        _fail(f"{label} is unsafe: {value!r}")
    candidate = paths.root.joinpath(*pure.parts)
    return _confined_path(paths, candidate, label)


def _expected_relative(paths: SourceArtifactPaths, path: Path) -> str:
    resolved = _confined_path(paths, path, "expected artifact")
    return resolved.relative_to(_root(paths)).as_posix()


def _expected_binding(
    *,
    manifest_sha256: str,
    provenance_sha256: str,
    archive_sha256: str,
    file_inventory_sha256: str,
) -> dict[str, Any]:
    return {
        "manifest_sha256": manifest_sha256,
        "provenance_sha256": provenance_sha256,
        "raw_archive_sha256": archive_sha256,
        "file_inventory_sha256": file_inventory_sha256,
        "raw_hash_recomputed_for_this_run": True,
    }


def _validate_binding(
    report: dict[str, Any],
    expected: dict[str, Any],
    label: str,
) -> None:
    binding = report.get("input_binding")
    if not isinstance(binding, dict):
        _fail(f"{label} has no input_binding object")
    for key, value in expected.items():
        if binding.get(key) != value:
            _fail(f"{label} input binding mismatch at {key}")


def _manifest_channels(
    manifest: dict[str, Any],
    contract: VerificationContract,
) -> dict[int, str]:
    try:
        records = manifest["archive_contract"]["channels"]
    except (KeyError, TypeError) as error:
        raise SourceArtifactError("manifest has no channel contract") from error
    if not isinstance(records, list):
        _fail("manifest channels are not a list")
    result: dict[int, str] = {}
    for record in records:
        if not isinstance(record, dict):
            _fail("manifest channel record is not an object")
        channel = _require_int(record.get("channel"), "manifest channel", minimum=1)
        accession = record.get("accession")
        if not isinstance(accession, str) or not accession:
            _fail(f"manifest channel {channel} has no accession")
        if channel in result:
            _fail(f"manifest repeats channel {channel}")
        result[channel] = accession
    expected = dict(zip(contract.channels, contract.accessions, strict=True))
    if result != expected:
        _fail("manifest channel/accession grid differs from the frozen contract")
    return result


def _validate_protocol_and_manifest(
    paths: SourceArtifactPaths,
    contract: VerificationContract,
) -> tuple[dict[str, Any], dict[str, Any], dict[int, str], str, str]:
    protocol, protocol_sha256 = _load_json(paths, paths.protocol, "protocol")
    manifest, manifest_sha256 = _load_json(paths, paths.manifest, "manifest")
    _, provenance_sha256 = _load_json(paths, paths.provenance, "provenance")

    if protocol_sha256 != contract.protocol_sha256:
        _fail("preregistered protocol SHA-256 changed")
    if manifest_sha256 != contract.manifest_sha256:
        _fail("acquisition manifest SHA-256 changed")
    if provenance_sha256 != contract.provenance_sha256:
        _fail("provenance SHA-256 changed")
    if protocol.get("experiment") != "107-jurkat-high-moi":
        _fail("protocol names a different experiment")
    expected_status = (
        "acquisition_complete_and_protocol_frozen_before_"
        "transcriptome_coordinate_values_are_parsed"
    )
    if protocol.get("status") != expected_status:
        _fail("protocol is not in the frozen pre-coordinate state")
    if manifest.get("dataset") != contract.dataset:
        _fail("manifest names a different dataset")

    immutable = protocol.get("immutable_inputs")
    if not isinstance(immutable, dict):
        _fail("protocol has no immutable_inputs object")
    expected_immutable = {
        "acquisition_manifest": manifest_sha256,
        "provenance": provenance_sha256,
    }
    for name, expected_sha256 in expected_immutable.items():
        record = immutable.get(name)
        if not isinstance(record, dict) or record.get("sha256") != expected_sha256:
            _fail(f"protocol immutable binding mismatch for {name}")

    try:
        archive_asset = manifest["assets"]["raw_archive"]
        inventory_asset = manifest["assets"]["file_inventory"]
        archive_contract = manifest["archive_contract"]
        modalities = archive_contract["modalities"]
    except (KeyError, TypeError) as error:
        raise SourceArtifactError("manifest is missing archive fields") from error
    if archive_asset.get("filename") != paths.archive.name:
        _fail("manifest raw archive filename differs from the configured path")
    if archive_asset.get("size_bytes") != contract.raw_archive_size_bytes:
        _fail("manifest raw archive size differs from the frozen contract")
    if inventory_asset.get("filename") != paths.file_inventory.name:
        _fail("manifest inventory filename differs from the configured path")
    if inventory_asset.get("sha256") != contract.file_inventory_sha256:
        _fail("manifest inventory SHA-256 differs from the frozen contract")
    if archive_contract.get("expected_member_count") != len(contract.channels) * 9:
        _fail("manifest member count differs from the frozen channel grid")

    expected_rows = contract.rows_by_modality()
    if set(modalities) != set(expected_rows):
        _fail("manifest modalities differ from the frozen contract")
    for modality, rows in expected_rows.items():
        record = modalities.get(modality)
        if not isinstance(record, dict):
            _fail(f"manifest modality {modality} is not an object")
        if record.get("feature_rows") != rows:
            _fail(f"manifest {modality} feature count changed")
        if record.get("barcode_columns") != contract.barcode_columns:
            _fail(f"manifest {modality} barcode count changed")
    if archive_contract.get("roles") != {
        "barcode": ".tsv.gz",
        "features": ".tsv.gz",
        "matrix": ".mtx.gz",
    }:
        _fail("manifest archive roles changed")
    safety = manifest.get("safety")
    if not isinstance(safety, dict):
        _fail("manifest has no safety object")
    if safety.get("never_scan_transcriptome_coordinates") is not True:
        _fail("manifest no longer forbids transcriptome coordinate scanning")
    if (
        safety.get("assignment_coordinate_order")
        != contract.manifest_assignment_coordinate_order
    ):
        _fail("manifest assignment coordinate order changed")

    channels = _manifest_channels(manifest, contract)
    return protocol, manifest, channels, manifest_sha256, provenance_sha256


def _validate_runtime_versions(protocol: dict[str, Any]) -> None:
    try:
        frozen_numpy = protocol["cell_and_assignment_freeze"][
            "transcriptome_cell_call"
        ]["rng_contract"]["numpy_version"]
        frozen_scipy = protocol["cell_and_assignment_freeze"]["guide_assignment"][
            "scipy_version"
        ]
    except (KeyError, TypeError) as error:
        raise SourceArtifactError(
            "protocol is missing frozen NumPy/SciPy runtime versions"
        ) from error
    if np.__version__ != frozen_numpy:
        _fail(
            f"NumPy version changed: current {np.__version__}, frozen {frozen_numpy}"
        )
    if scipy.__version__ != frozen_scipy:
        _fail(
            f"SciPy version changed: current {scipy.__version__}, frozen {frozen_scipy}"
        )


def _validate_source_preflight(
    paths: SourceArtifactPaths,
    protocol: dict[str, Any],
) -> tuple[str, str, str]:
    receipt, receipt_sha256 = _load_json(
        paths,
        paths.preflight_receipt,
        "source preflight receipt",
    )
    try:
        frozen_numpy = protocol["cell_and_assignment_freeze"][
            "transcriptome_cell_call"
        ]["rng_contract"]["numpy_version"]
        frozen_scipy = protocol["cell_and_assignment_freeze"]["guide_assignment"][
            "scipy_version"
        ]
    except (KeyError, TypeError) as error:
        raise SourceArtifactError(
            "protocol is missing frozen NumPy/SciPy runtime versions"
        ) from error
    preflight_sha256 = sha256_file(
        paths.preflight_code,
        paths=paths,
        label="source preflight code",
    )
    expected_scalars = {
        "status": "pass",
        "numpy_version": frozen_numpy,
        "scipy_version": frozen_scipy,
        "protocol_sha256": sha256_file(
            paths.protocol,
            paths=paths,
            label="protocol",
        ),
        "preflight_implementation_sha256": preflight_sha256,
        "held_outcomes_opened": False,
    }
    for field, expected in expected_scalars.items():
        if receipt.get(field) != expected:
            _fail(f"source preflight receipt mismatch at {field}")

    implementation_names = (
        "source_stage.py",
        "preprocess.py",
        "guide_calling.py",
        "demultiplex.py",
    )
    expected_implementations = {
        name: sha256_file(
            paths.root / name,
            paths=paths,
            label=f"source implementation {name}",
        )
        for name in implementation_names
    }
    if receipt.get("implementation_sha256") != expected_implementations:
        _fail("source preflight implementation hashes changed")
    test_names = (
        "test_preprocess.py",
        "test_source_stage.py",
        "test_guide_calling.py",
        "test_demultiplex.py",
    )
    expected_tests = {
        name: sha256_file(
            paths.root / name,
            paths=paths,
            label=f"source preflight test {name}",
        )
        for name in test_names
    }
    if receipt.get("test_module_sha256") != expected_tests:
        _fail("source preflight test hashes changed")
    if (
        not isinstance(receipt.get("tests_run"), int)
        or receipt["tests_run"] <= 0
        or receipt.get("failures") != 0
        or receipt.get("errors") != 0
        or receipt.get("skipped") != 0
    ):
        _fail("source preflight synthetic tests did not fully pass")
    _require_sha256(
        receipt.get("transcript_sha256"),
        "source preflight transcript hash",
    )
    return receipt_sha256, str(frozen_numpy), str(frozen_scipy)


def _validate_archive_and_inventory(
    paths: SourceArtifactPaths,
    contract: VerificationContract,
    manifest: dict[str, Any],
) -> tuple[str, str]:
    archive_path = _confined_path(paths, paths.archive, "raw archive")
    archive_size = archive_path.stat().st_size
    if archive_size != contract.raw_archive_size_bytes:
        _fail(
            "raw archive byte count changed: "
            f"{archive_size} != {contract.raw_archive_size_bytes}"
        )
    archive_sha256 = sha256_file(
        archive_path,
        paths=paths,
        label="raw archive",
    )
    if archive_sha256 != contract.raw_archive_sha256:
        _fail("raw archive SHA-256 changed")

    inventory_path = _confined_path(paths, paths.file_inventory, "file inventory")
    inventory_asset = manifest["assets"]["file_inventory"]
    if inventory_path.stat().st_size != inventory_asset.get("size_bytes"):
        _fail("file inventory byte count changed")
    inventory_sha256 = sha256_file(
        inventory_path,
        paths=paths,
        label="file inventory",
    )
    if inventory_sha256 != contract.file_inventory_sha256:
        _fail("file inventory SHA-256 changed")
    return archive_sha256, inventory_sha256


def _expected_member_names(
    manifest: dict[str, Any],
    channels: dict[int, str],
) -> set[str]:
    roles = manifest["archive_contract"]["roles"]
    modalities = manifest["archive_contract"]["modalities"]
    return {
        f"{accession}_channel{channel}_{modality}_{role}{suffix}"
        for channel, accession in channels.items()
        for modality in modalities
        for role, suffix in roles.items()
    }


def _validate_inventory_summary(
    inventory: Any,
    *,
    manifest: dict[str, Any],
    channels: dict[int, str],
    contract: VerificationContract,
    label: str,
) -> None:
    if not isinstance(inventory, dict):
        _fail(f"{label} inventory is not an object")
    if inventory.get("status") != "pass":
        _fail(f"{label} inventory did not pass")
    archive = inventory.get("archive")
    if not isinstance(archive, dict):
        _fail(f"{label} inventory has no archive object")
    if archive.get("filename") != manifest["assets"]["raw_archive"]["filename"]:
        _fail(f"{label} inventory names a different archive")
    if archive.get("size_bytes") != contract.raw_archive_size_bytes:
        _fail(f"{label} inventory archive size changed")
    expected_members = len(contract.channels) * 9
    if inventory.get("member_count") != expected_members:
        _fail(f"{label} inventory member count changed")
    if inventory.get("file_inventory_sha256") != contract.file_inventory_sha256:
        _fail(f"{label} inventory file-list hash changed")
    records = inventory.get("channels")
    if not isinstance(records, list):
        _fail(f"{label} inventory channels are not a list")
    observed: dict[int, tuple[str, int]] = {}
    for record in records:
        if not isinstance(record, dict):
            _fail(f"{label} inventory channel record is not an object")
        channel = _require_int(
            record.get("channel"),
            f"{label} inventory channel",
            minimum=1,
        )
        if channel in observed:
            _fail(f"{label} inventory repeats channel {channel}")
        observed[channel] = (
            record.get("accession"),
            record.get("member_count"),
        )
    expected = {channel: (accession, 9) for channel, accession in channels.items()}
    if observed != expected:
        _fail(f"{label} inventory channel grid changed")


@dataclass(frozen=True)
class _SupportTableSummary:
    rows: int
    total_nnz: int
    supported_features: int
    row_support_sha256: str


def _validate_support_table(
    paths: SourceArtifactPaths,
    path: Path,
    *,
    expected_size: int,
    expected_sha256: str,
    label: str,
) -> _SupportTableSummary:
    if path.stat().st_size != expected_size:
        _fail(f"{label} byte count changed")
    actual_sha256 = sha256_file(path, paths=paths, label=label)
    if actual_sha256 != expected_sha256:
        _fail(f"{label} SHA-256 changed")

    counts = array("Q")
    try:
        with _open_regular_file(paths, path, label) as raw:
            with gzip.GzipFile(fileobj=raw, mode="rb") as compressed:
                header = compressed.readline().decode("utf-8").rstrip("\r\n")
                if header != "row_index\tfeature_id\tfeature_name\tnnz_cells":
                    _fail(f"{label} has an unexpected header")
                for expected_index, raw_line in enumerate(compressed, start=1):
                    fields = raw_line.decode("utf-8").rstrip("\r\n").split("\t")
                    if len(fields) != 4:
                        _fail(f"{label} row {expected_index} is malformed")
                    try:
                        row_index = int(fields[0])
                        count = int(fields[3])
                    except ValueError as error:
                        raise SourceArtifactError(
                            f"{label} row {expected_index} is not numeric"
                        ) from error
                    if row_index != expected_index:
                        _fail(f"{label} row indices are not canonical")
                    if not fields[1] or not fields[2]:
                        _fail(f"{label} row {expected_index} has an empty feature")
                    if count < 0 or count > np.iinfo(np.uint64).max:
                        _fail(f"{label} row {expected_index} has an invalid count")
                    counts.append(count)
    except SourceArtifactError:
        raise
    except (EOFError, OSError, UnicodeDecodeError) as error:
        raise SourceArtifactError(
            f"{label} is not a valid complete gzip TSV"
        ) from error
    return _SupportTableSummary(
        rows=len(counts),
        total_nnz=sum(counts),
        supported_features=sum(count > 0 for count in counts),
        row_support_sha256=hashlib.sha256(counts.tobytes()).hexdigest(),
    )


def _validate_audit_reports(
    paths: SourceArtifactPaths,
    contract: VerificationContract,
    manifest: dict[str, Any],
    channels: dict[int, str],
    binding: dict[str, Any],
) -> dict[int, dict[str, Any]]:
    schema, _ = _load_json(paths, paths.schema_report, "schema report")
    gzip_report, _ = _load_json(paths, paths.gzip_report, "gzip report")
    support, _ = _load_json(paths, paths.support_report, "support report")

    _validate_binding(schema, binding, "schema report")
    _validate_binding(gzip_report, binding, "gzip report")
    _validate_binding(support, binding, "support report")

    if schema.get("status") != "pass":
        _fail("schema report did not pass")
    if schema.get("barcode_axes_verified") is not True:
        _fail("schema report did not verify barcode axes")
    if schema.get("assignment_coordinates_scanned") is not True:
        _fail("schema report did not scan assignment coordinates")
    if schema.get("transcriptome_coordinates_scanned") is not False:
        _fail("schema report claims transcriptome coordinates were scanned")
    if schema.get("support_completion_manifest") != _expected_relative(
        paths,
        paths.support_report,
    ):
        _fail("schema report points to a different support completion manifest")

    _validate_inventory_summary(
        schema.get("inventory"),
        manifest=manifest,
        channels=channels,
        contract=contract,
        label="schema report",
    )
    _validate_inventory_summary(
        gzip_report.get("inventory"),
        manifest=manifest,
        channels=channels,
        contract=contract,
        label="gzip report",
    )
    if schema.get("inventory") != gzip_report.get("inventory"):
        _fail("schema and gzip reports embed different archive inventories")

    if gzip_report.get("status") != "pass":
        _fail("gzip report did not pass")
    if gzip_report.get("all_members_read_to_gzip_eof") is not True:
        _fail("gzip report did not read every member to gzip EOF")
    expected_names = _expected_member_names(manifest, channels)
    if gzip_report.get("member_count") != len(expected_names):
        _fail("gzip report member count changed")
    gzip_members = gzip_report.get("members")
    if not isinstance(gzip_members, list) or len(gzip_members) != len(expected_names):
        _fail("gzip report member records are incomplete")
    gzip_by_name: dict[str, dict[str, Any]] = {}
    for record in gzip_members:
        if not isinstance(record, dict):
            _fail("gzip member record is not an object")
        name = record.get("name")
        if not isinstance(name, str) or name in gzip_by_name:
            _fail(f"gzip member name is invalid or duplicated: {name!r}")
        _require_int(
            record.get("uncompressed_bytes"),
            f"gzip member {name} uncompressed bytes",
            minimum=0,
        )
        _require_sha256(
            record.get("uncompressed_sha256"),
            f"gzip member {name} hash",
        )
        gzip_by_name[name] = record
    if set(gzip_by_name) != expected_names:
        _fail("gzip report member names differ from the frozen manifest grid")

    if support.get("status") != "complete":
        _fail("support report is incomplete")
    expected_support_count = len(contract.channels) * 2
    if support.get("expected_support_tables") != expected_support_count:
        _fail("support report expected-table count changed")
    support_records = support.get("support_tables")
    if not isinstance(support_records, list) or len(support_records) != (
        expected_support_count
    ):
        _fail("support report table records are incomplete")
    support_by_path: dict[str, dict[str, Any]] = {}
    support_summaries: dict[str, _SupportTableSummary] = {}
    for record in support_records:
        if not isinstance(record, dict):
            _fail("support table record is not an object")
        relative = record.get("path")
        if not isinstance(relative, str) or relative in support_by_path:
            _fail(f"support table path is invalid or duplicated: {relative!r}")
        actual_path = _relative_report_path(paths, relative, "support table")
        expected_size = _require_int(
            record.get("size_bytes"),
            f"support table {relative} size",
            minimum=0,
        )
        expected_sha256 = _require_sha256(
            record.get("sha256"),
            f"support table {relative} hash",
        )
        support_summaries[relative] = _validate_support_table(
            paths,
            actual_path,
            expected_size=expected_size,
            expected_sha256=expected_sha256,
            label=f"support table {relative}",
        )
        support_by_path[relative] = record

    rows_by_modality = contract.rows_by_modality()
    feature_axis = schema.get("feature_axis_sha256")
    if not isinstance(feature_axis, dict) or set(feature_axis) != set(rows_by_modality):
        _fail("schema feature-axis hash map is incomplete")
    for modality, value in feature_axis.items():
        _require_sha256(value, f"schema {modality} feature-axis hash")

    schema_channels = schema.get("channels")
    if not isinstance(schema_channels, list):
        _fail("schema channels are not a list")
    schema_by_channel: dict[int, dict[str, Any]] = {}
    schema_support_paths: set[str] = set()
    declared_assignment_nnz = 0
    for record in schema_channels:
        if not isinstance(record, dict):
            _fail("schema channel record is not an object")
        channel = _require_int(record.get("channel"), "schema channel", minimum=1)
        if channel in schema_by_channel:
            _fail(f"schema repeats channel {channel}")
        if record.get("accession") != channels.get(channel):
            _fail(f"schema accession changed for channel {channel}")
        barcode_axis = record.get("barcode_axis")
        if not isinstance(barcode_axis, dict):
            _fail(f"channel {channel} has no verified barcode axis")
        if barcode_axis.get("uncompressed_byte_identity_across_modalities") is not True:
            _fail(f"channel {channel} barcode axes are not byte-identical")
        if barcode_axis.get("gzip_crc_verified") is not True:
            _fail(f"channel {channel} barcode-axis gzip CRC was not verified")
        if barcode_axis.get("rows") != contract.barcode_columns:
            _fail(f"channel {channel} barcode-axis row count changed")
        if barcode_axis.get("expected_rows") != contract.barcode_columns:
            _fail(f"channel {channel} expected barcode row count changed")
        barcode_sha256 = _require_sha256(
            record.get("barcode_axis_sha256"),
            f"channel {channel} barcode-axis hash",
        )
        if barcode_axis.get("uncompressed_sha256") != barcode_sha256:
            _fail(f"channel {channel} barcode-axis hashes disagree")

        modalities = record.get("modalities")
        if not isinstance(modalities, dict) or set(modalities) != set(rows_by_modality):
            _fail(f"channel {channel} modality schema is incomplete")
        for modality, expected_rows in rows_by_modality.items():
            modality_record = modalities.get(modality)
            if not isinstance(modality_record, dict):
                _fail(f"channel {channel} {modality} schema is not an object")
            matrix = modality_record.get("matrix")
            features = modality_record.get("features")
            barcode = modality_record.get("barcode")
            if not all(
                isinstance(value, dict) for value in (matrix, features, barcode)
            ):
                _fail(f"channel {channel} {modality} schema is incomplete")
            if matrix.get("rows") != expected_rows:
                _fail(f"channel {channel} {modality} matrix rows changed")
            if matrix.get("columns") != contract.barcode_columns:
                _fail(f"channel {channel} {modality} matrix columns changed")
            if matrix.get("field") != "integer" or matrix.get("symmetry") != "general":
                _fail(f"channel {channel} {modality} matrix type changed")
            _require_int(
                matrix.get("nnz"),
                f"channel {channel} {modality} matrix nnz",
                minimum=0,
            )
            if features.get("rows") != expected_rows:
                _fail(f"channel {channel} {modality} feature rows changed")
            if features.get("gzip_crc_verified") is not True:
                _fail(f"channel {channel} {modality} feature gzip CRC failed")
            if features.get("uncompressed_sha256") != feature_axis[modality]:
                _fail(f"channel {channel} {modality} feature-axis hash changed")
            if barcode.get("rows") != contract.barcode_columns:
                _fail(f"channel {channel} {modality} barcode rows changed")
            if barcode.get("gzip_crc_verified") is not True:
                _fail(f"channel {channel} {modality} barcode gzip CRC failed")
            if barcode.get("uncompressed_sha256") != barcode_sha256:
                _fail(f"channel {channel} {modality} barcode hash changed")
            if barcode.get("exact_axis_check") != "uncompressed-byte identity":
                _fail(f"channel {channel} {modality} barcode check changed")

            prefix = f"{channels[channel]}_channel{channel}_{modality}_"
            for role, suffix in manifest["archive_contract"]["roles"].items():
                name = f"{prefix}{role}{suffix}"
                gzip_record = gzip_by_name[name]
                if role == "features":
                    expected_sha256 = features["uncompressed_sha256"]
                    if gzip_record["uncompressed_sha256"] != expected_sha256:
                        _fail(f"gzip feature hash disagrees for {name}")
                elif role == "barcode":
                    if gzip_record["uncompressed_sha256"] != barcode_sha256:
                        _fail(f"gzip barcode hash disagrees for {name}")

            if modality == "transcriptome":
                if "assignment_scan" in modality_record:
                    _fail(f"channel {channel} transcriptome coordinates were scanned")
                continue
            scan = modality_record.get("assignment_scan")
            if not isinstance(scan, dict):
                _fail(f"channel {channel} {modality} assignment scan is missing")
            if scan.get("gzip_crc_verified") is not True:
                _fail(f"channel {channel} {modality} assignment gzip CRC failed")
            if scan.get("coordinate_order") != contract.coordinate_order:
                _fail(f"channel {channel} {modality} coordinate order changed")
            matrix_nnz = matrix["nnz"]
            if scan.get("declared_nnz") != matrix_nnz:
                _fail(f"channel {channel} {modality} declared nnz changed")
            if scan.get("observed_nnz") != matrix_nnz:
                _fail(f"channel {channel} {modality} observed nnz changed")
            declared_assignment_nnz += matrix_nnz
            support_record = scan.get("support_table")
            if not isinstance(support_record, dict):
                _fail(f"channel {channel} {modality} support record is missing")
            support_path = support_record.get("path")
            expected_path = (
                f"audit/support/channel{channel:02d}_{modality}_row_support.tsv.gz"
            )
            if support_path != expected_path:
                _fail(f"channel {channel} {modality} support path changed")
            if support_record != support_by_path.get(expected_path):
                _fail(
                    f"channel {channel} {modality} support bindings disagree "
                    "between reports"
                )
            support_summary = support_summaries[expected_path]
            if support_summary.rows != expected_rows:
                _fail(f"channel {channel} {modality} support row count changed")
            if support_summary.total_nnz != matrix_nnz:
                _fail(f"channel {channel} {modality} support nnz changed")
            if support_summary.supported_features != scan.get("supported_features"):
                _fail(f"channel {channel} {modality} supported-feature count changed")
            if support_summary.row_support_sha256 != scan.get("row_support_sha256"):
                _fail(f"channel {channel} {modality} row-support hash changed")
            schema_support_paths.add(expected_path)
        schema_by_channel[channel] = record

    if set(schema_by_channel) != set(contract.channels):
        _fail("schema report does not contain exactly the frozen channels")
    if schema_support_paths != set(support_by_path):
        _fail("schema and support reports name different support tables")
    if schema.get("declared_assignment_nnz_total") != declared_assignment_nnz:
        _fail("schema assignment nnz total is inconsistent")
    maximum_nnz = manifest.get("safety", {}).get("max_assignment_nnz_total")
    if not isinstance(maximum_nnz, int) or declared_assignment_nnz > maximum_nnz:
        _fail("schema assignment nnz exceeds the frozen safety limit")
    return schema_by_channel


def _validate_oracle(
    paths: SourceArtifactPaths,
    protocol: dict[str, Any],
    contract: VerificationContract,
) -> tuple[str, str, str]:
    receipt, receipt_sha256 = _load_json(
        paths,
        paths.oracle_receipt,
        "R oracle receipt",
    )
    if receipt.get("status") != "pass":
        _fail("R oracle receipt did not pass")
    checks = receipt.get("checks")
    if (
        not isinstance(checks, list)
        or len(checks) != len(set(checks))
        or set(checks) != set(contract.required_oracle_checks)
    ):
        _fail("R oracle receipt check set changed or is incomplete")

    try:
        semantic = protocol["immutable_inputs"]["emptydrops_semantic_reference"]
        good_turing = semantic["good_turing_reference"]
    except (KeyError, TypeError) as error:
        raise SourceArtifactError(
            "protocol has no pinned EmptyDrops semantic reference"
        ) from error
    expected_fields = {
        "DropletUtils_version": semantic.get("package_version"),
        "DropletUtils_commit": semantic.get("commit"),
        "DropletUtils_empty_drops_sha256": semantic.get("empty_drops_source_sha256"),
        "DropletUtils_barcode_ranks_sha256": semantic.get(
            "barcode_ranks_source_sha256"
        ),
        "DropletUtils_ambient_profile_sha256": semantic.get(
            "ambient_profile_source_sha256"
        ),
        "edgeR_version": good_turing.get("package_version"),
        "edgeR_commit": good_turing.get("commit"),
        "edgeR_good_turing_r_sha256": good_turing.get("r_source_sha256"),
        "edgeR_good_turing_c_sha256": good_turing.get("c_source_sha256"),
    }
    for field, expected in expected_fields.items():
        if receipt.get(field) != expected:
            _fail(f"R oracle receipt mismatch at {field}")
    for field in expected_fields:
        if field.endswith("_sha256"):
            _require_sha256(receipt[field], f"R oracle receipt {field}")

    preprocess_sha256 = sha256_file(
        paths.preprocess_code,
        paths=paths,
        label="preprocess code",
    )
    r_oracle_sha256 = sha256_file(
        paths.r_oracle_code,
        paths=paths,
        label="R oracle code",
    )
    if receipt.get("python_port_sha256") != preprocess_sha256:
        _fail("R oracle receipt is stale for preprocess.py")
    if receipt.get("r_oracle_sha256") != r_oracle_sha256:
        _fail("R oracle receipt is stale for r_oracle.R")
    _require_int(receipt.get("synthetic_nnz"), "R oracle synthetic nnz", minimum=1)
    shape = receipt.get("synthetic_shape")
    if (
        not isinstance(shape, list)
        or len(shape) != 2
        or any(
            isinstance(value, bool) or not isinstance(value, int) or value <= 0
            for value in shape
        )
    ):
        _fail("R oracle synthetic shape is invalid")
    return receipt_sha256, preprocess_sha256, r_oracle_sha256


def _load_npz(
    paths: SourceArtifactPaths,
    path: Path,
    expected_sha256: Any,
    label: str,
) -> dict[str, np.ndarray]:
    expected = _require_sha256(expected_sha256, f"{label} expected hash")
    try:
        with _open_regular_file(paths, path, label) as handle:
            actual = _hash_handle(handle)
            if actual != expected:
                _fail(f"{label} SHA-256 changed")
            handle.seek(0)
            with np.load(handle, allow_pickle=False) as payload:
                return {name: payload[name] for name in payload.files}
    except SourceArtifactError:
        raise
    except (OSError, ValueError) as error:
        raise SourceArtifactError(f"cannot load {label}: {error}") from error


def _require_array_keys(
    arrays: dict[str, np.ndarray],
    expected: frozenset[str],
    label: str,
) -> None:
    if set(arrays) != expected or len(arrays) != len(expected):
        _fail(f"{label} array keys changed: {sorted(arrays)}")


def _require_dtype(
    array: np.ndarray,
    dtype: str,
    label: str,
) -> None:
    if array.dtype != np.dtype(dtype):
        _fail(f"{label} dtype changed: {array.dtype} != {np.dtype(dtype)}")


def _require_vector(
    array: np.ndarray,
    length: int,
    label: str,
) -> None:
    if array.ndim != 1 or array.shape != (length,):
        _fail(f"{label} shape changed: {array.shape} != {(length,)}")


def _require_sorted_indices(
    array: np.ndarray,
    upper_bound: int,
    label: str,
) -> None:
    if array.ndim != 1 or array.dtype != np.dtype("int32"):
        _fail(f"{label} must be a one-dimensional int32 vector")
    if array.size:
        if int(array[0]) < 0 or int(array[-1]) >= upper_bound:
            _fail(f"{label} contains an out-of-range index")
        if np.any(np.diff(array.astype(np.int64, copy=False)) <= 0):
            _fail(f"{label} is not strictly increasing and unique")


@dataclass(frozen=True)
class _PreparedArrays:
    totals: np.ndarray
    tested_indices: np.ndarray


def _validate_prepare_arrays(
    paths: SourceArtifactPaths,
    contract: VerificationContract,
    arrays_path: Path,
    report: dict[str, Any],
) -> _PreparedArrays:
    arrays = _load_npz(
        paths,
        arrays_path,
        report.get("arrays_sha256"),
        "prepared arrays",
    )
    _require_array_keys(arrays, PREPARE_ARRAY_KEYS, "prepared arrays")
    columns = contract.barcode_columns
    rows = contract.rows_by_modality()["transcriptome"]

    exact_specs = {
        "totals": (columns, "int64"),
        "detected_features": (columns, "uint16"),
        "mitochondrial_counts": (columns, "int64"),
        "feature_totals": (rows, "int64"),
        "ambient_counts": (rows, "int64"),
        "ambient_probability": (rows, "float64"),
        "observed_data_log_probability": (
            _require_int(
                report.get("tested_barcodes"),
                "prepare report tested_barcodes",
                minimum=0,
            ),
            "float64",
        ),
    }
    for name, (length, dtype) in exact_specs.items():
        _require_vector(arrays[name], length, f"prepared {name}")
        _require_dtype(arrays[name], dtype, f"prepared {name}")
    tested = arrays["tested_indices"]
    _require_sorted_indices(tested, columns, "prepared tested_indices")
    if tested.size != report["tested_barcodes"]:
        _fail("prepare report tested count disagrees with tested_indices")

    totals = arrays["totals"]
    detected = arrays["detected_features"]
    mitochondrial = arrays["mitochondrial_counts"]
    feature_totals = arrays["feature_totals"]
    ambient_counts = arrays["ambient_counts"]
    ambient_probability = arrays["ambient_probability"]
    observed = arrays["observed_data_log_probability"]
    if np.any(totals < 0):
        _fail("prepared totals contain negative values")
    if np.any(detected > rows):
        _fail("prepared detected_features exceeds the feature count")
    if np.any(mitochondrial < 0) or np.any(mitochondrial > totals):
        _fail("prepared mitochondrial counts are invalid")
    if np.any(feature_totals < 0) or np.any(ambient_counts < 0):
        _fail("prepared feature totals contain negative values")
    if np.any(ambient_counts > feature_totals):
        _fail("prepared ambient counts exceed feature totals")
    if int(feature_totals.sum(dtype=np.int64)) != int(totals.sum(dtype=np.int64)):
        _fail("prepared row and column totals disagree")
    expected_ambient_total = int(totals[totals <= contract.lower].sum(dtype=np.int64))
    if int(ambient_counts.sum(dtype=np.int64)) != expected_ambient_total:
        _fail("prepared ambient counts disagree with the lower-threshold barcodes")
    if (
        np.any(~np.isfinite(ambient_probability))
        or np.any(ambient_probability < 0)
        or not np.isclose(
            ambient_probability.sum(dtype=np.float64),
            1.0,
            rtol=0,
            atol=1e-12,
        )
    ):
        _fail("prepared ambient probability is invalid")
    if np.any(~np.isfinite(observed)):
        _fail("prepared observed log probabilities are not finite")
    expected_tested = np.flatnonzero(totals > contract.lower).astype(np.int32)
    if not np.array_equal(tested, expected_tested):
        _fail("prepared tested_indices disagree with the frozen lower threshold")

    expected_report_values = {
        "nonzero_barcodes": int(np.count_nonzero(totals)),
        "maximum_total": int(totals.max(initial=0)),
        "ambient_barcodes": int(np.count_nonzero(totals <= contract.lower)),
        "expressed_features": int(np.count_nonzero(feature_totals)),
    }
    for field, expected in expected_report_values.items():
        if report.get(field) != expected:
            _fail(f"prepare report {field} disagrees with prepared arrays")
    mitochondrial_features = _require_int(
        report.get("mitochondrial_features"),
        "prepare report mitochondrial_features",
        minimum=0,
    )
    if mitochondrial_features > rows:
        _fail("prepare report mitochondrial feature count is impossible")
    return _PreparedArrays(
        totals=totals,
        tested_indices=tested,
    )


def _benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    count = p_values.size
    order = np.argsort(p_values, kind="stable")
    ranked = p_values[order] * count / np.arange(1, count + 1)
    adjusted_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(p_values)
    adjusted[order] = np.minimum(adjusted_sorted, 1.0)
    return adjusted


def _validate_call_arrays(
    paths: SourceArtifactPaths,
    contract: VerificationContract,
    output_path: Path,
    report: dict[str, Any],
    prepared: _PreparedArrays,
    knee: float,
) -> tuple[int, int]:
    arrays = _load_npz(
        paths,
        output_path,
        report.get("output_sha256"),
        "production call output",
    )
    _require_array_keys(arrays, CALL_ARRAY_KEYS, "production call output")
    tested = arrays["tested_indices"]
    accepted = arrays["accepted_indices"]
    _require_sorted_indices(tested, contract.barcode_columns, "call tested_indices")
    _require_sorted_indices(accepted, contract.barcode_columns, "accepted_indices")
    if not np.array_equal(tested, prepared.tested_indices):
        _fail("production call tested_indices differ from prepared tested_indices")
    tested_count = tested.size
    for name in ("p_value", "fdr", "limited"):
        _require_vector(arrays[name], tested_count, f"production {name}")
    _require_dtype(arrays["p_value"], "float64", "production p_value")
    _require_dtype(arrays["fdr"], "float64", "production fdr")
    _require_dtype(arrays["limited"], "bool", "production limited")

    p_value = arrays["p_value"]
    fdr = arrays["fdr"]
    limited = arrays["limited"]
    if np.any(~np.isfinite(p_value)) or np.any((p_value <= 0) | (p_value > 1)):
        _fail("production p-values are outside (0, 1]")
    if np.any(~np.isfinite(fdr)) or np.any((fdr < 0) | (fdr > 1)):
        _fail("production FDR values are outside [0, 1]")
    denominator = contract.production_iterations + 1
    numerators = p_value * denominator
    if np.any(~np.isclose(numerators, np.rint(numerators), rtol=0, atol=1e-10)):
        _fail("production p-values are not on the frozen Monte Carlo grid")
    if not np.array_equal(limited, np.rint(numerators).astype(np.int64) == 1):
        _fail("production limited flags disagree with p-values")

    tested_totals = prepared.totals[tested]
    temporary = p_value.copy()
    forced = tested_totals >= knee
    temporary[forced] = 0.0
    expected_fdr = _benjamini_hochberg(temporary)
    if not np.array_equal(fdr, expected_fdr):
        _fail("production FDR values do not reproduce frozen BH behavior")
    expected_accepted = tested[fdr <= contract.fdr_threshold]
    if not np.array_equal(accepted, expected_accepted):
        _fail("accepted_indices disagree with the frozen FDR threshold")

    expected_report_values = {
        "tested_barcodes": tested_count,
        "accepted_barcodes": accepted.size,
        "forced_retained_barcodes": int(np.count_nonzero(forced)),
        "limited_barcodes": int(np.count_nonzero(limited)),
    }
    for field, expected in expected_report_values.items():
        if report.get(field) != expected:
            _fail(f"production report {field} disagrees with output arrays")
    return tested_count, accepted.size


def _validate_channels(
    paths: SourceArtifactPaths,
    contract: VerificationContract,
    schema_by_channel: dict[int, dict[str, Any]],
    *,
    archive_sha256: str,
    oracle_receipt_sha256: str,
    preprocess_sha256: str,
    preflight_receipt_sha256: str,
    numpy_version: str,
    scipy_version: str,
) -> list[dict[str, int]]:
    source_stage_sha256 = sha256_file(
        paths.source_stage_code,
        paths=paths,
        label="source stage code",
    )
    summaries: list[dict[str, int]] = []
    for channel in contract.channels:
        directory = paths.work / f"channel{channel:02d}"
        prepare_path = directory / "prepare-report.json"
        arrays_path = directory / "cell-call-inputs.npz"
        call_path = directory / "cell-calls-production.json"
        output_path = directory / "cell-calls-production.npz"

        prepare, prepare_sha256 = _load_json(
            paths,
            prepare_path,
            f"channel {channel} prepare report",
        )
        schema = schema_by_channel[channel]
        if prepare.get("status") != "complete":
            _fail(f"channel {channel} prepare report is incomplete")
        if prepare.get("channel") != channel:
            _fail(f"channel {channel} prepare report channel changed")
        if prepare.get("accession") != schema.get("accession"):
            _fail(f"channel {channel} prepare accession changed")
        expected_matrix_member = (
            f"{schema['accession']}_channel{channel}_transcriptome_matrix.mtx.gz"
        )
        if prepare.get("matrix_member") != expected_matrix_member:
            _fail(f"channel {channel} prepare matrix member changed")
        if prepare.get("raw_archive_sha256") != archive_sha256:
            _fail(f"channel {channel} prepare used a different raw archive")
        expected_code_hashes = {
            "source_stage_sha256": source_stage_sha256,
            "preprocess_sha256": preprocess_sha256,
            "r_oracle_receipt_sha256": oracle_receipt_sha256,
            "preflight_receipt_sha256": preflight_receipt_sha256,
            "preflight_implementation_sha256": sha256_file(
                paths.preflight_code,
                paths=paths,
                label="source preflight code",
            ),
            "numpy_version": numpy_version,
            "scipy_version": scipy_version,
        }
        for field, expected in expected_code_hashes.items():
            if prepare.get(field) != expected:
                _fail(f"channel {channel} prepare report is stale at {field}")
        if prepare.get("arrays_path") != _expected_relative(paths, arrays_path):
            _fail(f"channel {channel} prepare arrays path changed")
        matrix = schema["modalities"]["transcriptome"]["matrix"]
        if prepare.get("shape") != [matrix["rows"], matrix["columns"]]:
            _fail(f"channel {channel} prepare shape changed")
        if prepare.get("nnz") != matrix["nnz"]:
            _fail(f"channel {channel} prepare nnz changed")
        if prepare.get("coordinate_order") != contract.coordinate_order:
            _fail(f"channel {channel} prepare coordinate order changed")
        knee = _require_number(
            prepare.get("knee"),
            f"channel {channel} prepare knee",
            minimum=float(contract.lower),
        )
        _require_number(
            prepare.get("inflection"),
            f"channel {channel} prepare inflection",
            minimum=0.0,
        )
        _require_number(
            prepare.get("elapsed_seconds"),
            f"channel {channel} prepare elapsed_seconds",
            minimum=0.0,
        )
        prepared = _validate_prepare_arrays(
            paths,
            contract,
            arrays_path,
            prepare,
        )

        call, _ = _load_json(
            paths,
            call_path,
            f"channel {channel} production call report",
        )
        if call.get("status") != "complete" or call.get("stage") != "production":
            _fail(f"channel {channel} production call report is incomplete")
        if call.get("channel") != channel:
            _fail(f"channel {channel} production call report channel changed")
        if call.get("iterations") != contract.production_iterations:
            _fail(
                f"channel {channel} did not use exactly "
                f"{contract.production_iterations} production iterations"
            )
        _require_int(
            call.get("iteration_block_size"),
            f"channel {channel} iteration block size",
            minimum=1,
        )
        if call.get("prepare_report_sha256") != prepare_sha256:
            _fail(
                f"channel {channel} production call is stale after prepare report "
                "replacement"
            )
        arrays_sha256 = _require_sha256(
            prepare.get("arrays_sha256"),
            f"channel {channel} prepared arrays hash",
        )
        if call.get("prepared_arrays_sha256") != arrays_sha256:
            _fail(
                f"channel {channel} production call is stale after prepared arrays "
                "replacement"
            )
        for field, expected in expected_code_hashes.items():
            if call.get(field) != expected:
                _fail(f"channel {channel} production call is stale at {field}")
        if call.get("output_path") != _expected_relative(paths, output_path):
            _fail(f"channel {channel} production output path changed")
        call_knee = _require_number(
            call.get("knee"),
            f"channel {channel} production knee",
            minimum=float(contract.lower),
        )
        if call_knee != knee:
            _fail(f"channel {channel} production knee differs from prepare")
        _require_number(
            call.get("elapsed_seconds"),
            f"channel {channel} call elapsed_seconds",
            minimum=0.0,
        )
        tested, accepted = _validate_call_arrays(
            paths,
            contract,
            output_path,
            call,
            prepared,
            knee,
        )
        summaries.append(
            {
                "channel": channel,
                "tested_barcodes": tested,
                "accepted_barcodes": accepted,
            }
        )
    return summaries


def verify_source_artifacts(
    paths: SourceArtifactPaths | None = None,
    contract: VerificationContract = PRODUCTION_CONTRACT,
) -> dict[str, Any]:
    """Verify every source artifact required before assignment or modeling."""

    paths = paths or SourceArtifactPaths.from_root(HERE)
    if len(contract.channels) != len(contract.accessions):
        _fail("verification contract channel/accession lengths differ")
    for field, value in {
        "protocol": contract.protocol_sha256,
        "manifest": contract.manifest_sha256,
        "provenance": contract.provenance_sha256,
        "raw archive": contract.raw_archive_sha256,
        "file inventory": contract.file_inventory_sha256,
    }.items():
        _require_sha256(value, f"frozen {field} hash")

    (
        protocol,
        manifest,
        channels,
        manifest_sha256,
        provenance_sha256,
    ) = _validate_protocol_and_manifest(paths, contract)
    _validate_runtime_versions(protocol)
    (
        preflight_receipt_sha256,
        numpy_version,
        scipy_version,
    ) = _validate_source_preflight(paths, protocol)
    archive_sha256, inventory_sha256 = _validate_archive_and_inventory(
        paths,
        contract,
        manifest,
    )
    binding = _expected_binding(
        manifest_sha256=manifest_sha256,
        provenance_sha256=provenance_sha256,
        archive_sha256=archive_sha256,
        file_inventory_sha256=inventory_sha256,
    )
    schema_by_channel = _validate_audit_reports(
        paths,
        contract,
        manifest,
        channels,
        binding,
    )
    (
        oracle_receipt_sha256,
        preprocess_sha256,
        r_oracle_sha256,
    ) = _validate_oracle(paths, protocol, contract)
    channel_summaries = _validate_channels(
        paths,
        contract,
        schema_by_channel,
        archive_sha256=archive_sha256,
        oracle_receipt_sha256=oracle_receipt_sha256,
        preprocess_sha256=preprocess_sha256,
        preflight_receipt_sha256=preflight_receipt_sha256,
        numpy_version=numpy_version,
        scipy_version=scipy_version,
    )
    return {
        "status": "pass",
        "dataset": contract.dataset,
        "raw_archive_size_bytes": contract.raw_archive_size_bytes,
        "raw_archive_sha256": archive_sha256,
        "manifest_sha256": manifest_sha256,
        "provenance_sha256": provenance_sha256,
        "file_inventory_sha256": inventory_sha256,
        "oracle_receipt_sha256": oracle_receipt_sha256,
        "preflight_receipt_sha256": preflight_receipt_sha256,
        "preprocess_sha256": preprocess_sha256,
        "r_oracle_sha256": r_oracle_sha256,
        "numpy_version": numpy_version,
        "scipy_version": scipy_version,
        "channels": channel_summaries,
        "total_tested_barcodes": sum(
            record["tested_barcodes"] for record in channel_summaries
        ),
        "total_accepted_barcodes": sum(
            record["accepted_barcodes"] for record in channel_summaries
        ),
        "transcriptome_matrix_values_read": False,
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--root",
        type=Path,
        default=HERE,
        help="Experiment 107 root; production constants remain mandatory.",
    )
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        result = verify_source_artifacts(SourceArtifactPaths.from_root(args.root))
    except SourceArtifactError as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
