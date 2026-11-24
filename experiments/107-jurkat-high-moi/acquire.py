#!/usr/bin/env python3
"""Fail-closed, resumable acquisition for the frozen GSE249595 assets."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable, Iterator


HERE = Path(__file__).resolve().parent
DEFAULT_MANIFEST = HERE / "ACQUISITION_MANIFEST.json"
DEFAULT_DOWNLOADS = HERE / "downloads"
RECEIPT_NAME = "download_receipt.json"
RECEIPT_SCHEMA_VERSION = 2
RESUME_SCHEMA_VERSION = 1
EXPECTED_ASSET_GROUPS = {
    "file_inventory": "metadata",
    "guide_reference": "metadata",
    "supplementary_table": "metadata",
    "raw_archive": "raw",
}
EXPECTED_SAFETY_VALUES = {
    "download_root": "downloads",
    "report_root": "audit",
    "assignment_coordinate_order": "row_or_column_major_strict",
    "never_scan_transcriptome_coordinates": True,
    "never_extract_archive": True,
}


class AcquisitionError(RuntimeError):
    """An acquisition contract or local artifact failed validation."""


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    file_descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(file_descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def require_under_experiment(path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_relative_to(HERE):
        raise AcquisitionError(f"path escapes experiment directory: {resolved}")
    return resolved


def load_manifest(path: Path = DEFAULT_MANIFEST) -> dict[str, Any]:
    require_under_experiment(path)
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AcquisitionError(f"cannot read manifest {path}: {error}") from error

    if manifest.get("schema_version") != 1:
        raise AcquisitionError("unsupported acquisition manifest schema")
    if manifest.get("dataset") != "GSE249595":
        raise AcquisitionError("manifest dataset must be GSE249595")

    assets = manifest.get("assets")
    if not isinstance(assets, dict) or set(assets) != set(EXPECTED_ASSET_GROUPS):
        raise AcquisitionError("manifest asset set is incomplete or unexpected")

    filenames: set[str] = set()
    for name, asset in assets.items():
        expected_group = EXPECTED_ASSET_GROUPS[name]
        if asset.get("group") != expected_group:
            raise AcquisitionError(f"{name}: asset group must be {expected_group!r}")
        filename = asset.get("filename")
        if not isinstance(filename, str) or Path(filename).name != filename:
            raise AcquisitionError(f"{name}: filename must be a plain basename")
        if filename in filenames:
            raise AcquisitionError(f"duplicate destination filename: {filename}")
        filenames.add(filename)
        url = asset.get("url")
        if not isinstance(url, str) or not url.startswith("https://"):
            raise AcquisitionError(f"{name}: only HTTPS source URLs are allowed")
        size = asset.get("size_bytes")
        if not isinstance(size, int) or size <= 0:
            raise AcquisitionError(f"{name}: invalid expected byte count")
        expected_hash = asset.get("sha256")
        if expected_hash is not None and (
            not isinstance(expected_hash, str)
            or len(expected_hash) != 64
            or any(character not in "0123456789abcdef" for character in expected_hash)
        ):
            raise AcquisitionError(f"{name}: invalid SHA-256")
        if asset["group"] == "metadata" and expected_hash is None:
            raise AcquisitionError(
                f"{name}: metadata must have a source-pinned SHA-256"
            )
        if name == "raw_archive":
            if expected_hash is not None:
                raise AcquisitionError(
                    "raw_archive SHA-256 must remain null until source-pinned"
                )
            if asset.get("sha256_status") != "not_source_pinned":
                raise AcquisitionError(
                    "raw_archive must declare sha256_status=not_source_pinned"
                )

    contract = manifest.get("archive_contract")
    if not isinstance(contract, dict):
        raise AcquisitionError("manifest has no archive contract")
    channels = contract.get("channels")
    if not isinstance(channels, list) or not channels:
        raise AcquisitionError("archive contract has no channels")
    channel_numbers = [record.get("channel") for record in channels]
    accessions = [record.get("accession") for record in channels]
    if channel_numbers != list(range(1, len(channels) + 1)):
        raise AcquisitionError("archive channels must be unique and sequential from 1")
    if any(not isinstance(accession, str) for accession in accessions) or len(
        set(accessions)
    ) != len(accessions):
        raise AcquisitionError("archive accessions must be unique strings")

    modalities = contract.get("modalities")
    if not isinstance(modalities, dict) or set(modalities) != {
        "transcriptome",
        "guides",
        "labels",
    }:
        raise AcquisitionError("archive modality contract is incomplete")
    barcode_columns: set[int] = set()
    for modality, dimensions in modalities.items():
        if not isinstance(dimensions, dict):
            raise AcquisitionError(f"{modality}: invalid dimension contract")
        rows = dimensions.get("feature_rows")
        columns = dimensions.get("barcode_columns")
        if (
            not isinstance(rows, int)
            or rows <= 0
            or not isinstance(columns, int)
            or columns <= 0
        ):
            raise AcquisitionError(f"{modality}: dimensions must be positive integers")
        barcode_columns.add(columns)
    if len(barcode_columns) != 1:
        raise AcquisitionError("all modalities must share one barcode dimension")

    roles = contract.get("roles")
    if roles != {
        "barcode": ".tsv.gz",
        "features": ".tsv.gz",
        "matrix": ".mtx.gz",
    }:
        raise AcquisitionError("archive role/suffix contract is invalid")
    expected_member_count = len(channels) * len(modalities) * len(roles)
    if contract.get("expected_member_count") != expected_member_count:
        raise AcquisitionError(
            "expected_member_count does not match channels x modalities x roles"
        )

    safety = manifest.get("safety")
    if not isinstance(safety, dict):
        raise AcquisitionError("manifest has no safety contract")
    maximum_nnz = safety.get("max_assignment_nnz_total")
    if not isinstance(maximum_nnz, int) or maximum_nnz <= 0:
        raise AcquisitionError("invalid assignment-nnz safety ceiling")
    for key, expected in EXPECTED_SAFETY_VALUES.items():
        if safety.get(key) != expected:
            raise AcquisitionError(
                f"manifest safety field {key!r} must be {expected!r}"
            )
    if not isinstance(manifest.get("provenance_expectations"), dict):
        raise AcquisitionError("manifest has no provenance expectations")

    return manifest


def verify_local_asset(path: Path, asset: dict[str, Any]) -> dict[str, Any]:
    if path.is_symlink():
        raise AcquisitionError(f"refusing symlinked downloaded asset: {path}")
    if not path.is_file():
        raise AcquisitionError(f"missing downloaded asset: {path}")
    actual_size = path.stat().st_size
    expected_size = int(asset["size_bytes"])
    if actual_size != expected_size:
        raise AcquisitionError(
            f"{path.name}: byte count {actual_size} != expected {expected_size}"
        )
    actual_hash = sha256_file(path)
    expected_hash = asset.get("sha256")
    if expected_hash is not None and actual_hash != expected_hash:
        raise AcquisitionError(
            f"{path.name}: SHA-256 {actual_hash} != expected {expected_hash}"
        )
    return {
        "filename": path.name,
        "size_bytes": actual_size,
        "sha256": actual_hash,
        "source_sha256_pinned": expected_hash is not None,
    }


def load_receipt(download_dir: Path) -> dict[str, Any]:
    download_dir = require_under_experiment(download_dir)
    receipt_path = download_dir / RECEIPT_NAME
    if receipt_path.is_symlink():
        raise AcquisitionError(f"download receipt cannot be a symlink: {receipt_path}")
    if not receipt_path.exists():
        return {"schema_version": RECEIPT_SCHEMA_VERSION, "assets": {}}
    if not receipt_path.is_file():
        raise AcquisitionError(
            f"download receipt is not a regular file: {receipt_path}"
        )
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AcquisitionError(f"invalid download receipt: {error}") from error
    if receipt.get("schema_version") != RECEIPT_SCHEMA_VERSION or not isinstance(
        receipt.get("assets"), dict
    ):
        raise AcquisitionError("unsupported or malformed download receipt")
    return receipt


def asset_binding(
    name: str,
    asset: dict[str, Any],
    manifest_sha256: str,
) -> dict[str, Any]:
    """Return the immutable source contract carried by one receipt entry."""
    return {
        "asset_name": name,
        "manifest_sha256": manifest_sha256,
        "group": asset["group"],
        "source_url": asset["url"],
        "expected_filename": asset["filename"],
        "expected_size_bytes": asset["size_bytes"],
        "expected_sha256": asset.get("sha256"),
        "sha256_status": asset.get(
            "sha256_status",
            "source_pinned" if asset.get("sha256") is not None else None,
        ),
    }


def receipt_entry_matches_binding(
    entry: Any,
    name: str,
    asset: dict[str, Any],
    manifest_sha256: str,
) -> bool:
    return isinstance(entry, dict) and entry.get("binding") == asset_binding(
        name, asset, manifest_sha256
    )


def select_assets(
    manifest: dict[str, Any],
    group: str,
    confirm_large_download: bool,
) -> list[tuple[str, dict[str, Any]]]:
    if group not in {"metadata", "raw", "all"}:
        raise AcquisitionError(f"unsupported download group: {group}")
    selected = [
        (name, asset)
        for name, asset in manifest["assets"].items()
        if group == "all" or asset["group"] == group
    ]
    includes_raw = any(name == "raw_archive" for name, _ in selected)
    if includes_raw and not confirm_large_download:
        raise AcquisitionError(
            "raw acquisition requires --confirm-large-download "
            "(expected 3,753,820,160 bytes)"
        )
    return sorted(
        selected,
        key=lambda item: (item[1]["group"] == "raw", item[0]),
    )


def curl_base_command() -> list[str]:
    return [
        "curl",
        "--disable",
        "--fail",
        "--location",
        "--show-error",
        "--proto",
        "=https",
        "--proto-redir",
        "=https",
        "--retry",
        "5",
        "--retry-all-errors",
        "--connect-timeout",
        "30",
    ]


def curl_probe_command(asset: dict[str, Any], headers_path: Path) -> list[str]:
    return [
        *curl_base_command(),
        "--silent",
        "--head",
        "--dump-header",
        str(headers_path),
        "--output",
        "/dev/null",
        str(asset["url"]),
    ]


def curl_command(
    asset: dict[str, Any],
    partial_path: Path,
    headers_path: Path,
    *,
    if_range: str | None = None,
) -> list[str]:
    command = [
        *curl_base_command(),
        "--max-filesize",
        str(asset["size_bytes"]),
        "--dump-header",
        str(headers_path),
        "--continue-at",
        "-",
        "--output",
        str(partial_path),
    ]
    if if_range is not None:
        command.extend(["--header", f"If-Range: {if_range}"])
    command.append(str(asset["url"]))
    return command


def parse_final_http_headers(path: Path) -> dict[str, Any]:
    try:
        payload = path.read_text(encoding="iso-8859-1")
    except OSError as error:
        raise AcquisitionError(f"cannot read curl response headers: {error}") from error
    blocks = [
        block
        for block in re.split(r"\r?\n\r?\n", payload)
        if block.lstrip().startswith("HTTP/")
    ]
    if not blocks:
        raise AcquisitionError("curl response did not contain an HTTP header block")
    lines = blocks[-1].splitlines()
    status_fields = lines[0].split()
    if len(status_fields) < 2 or not status_fields[1].isdigit():
        raise AcquisitionError("curl response contained an invalid HTTP status line")
    status_code = int(status_fields[1])
    if not 200 <= status_code < 300:
        raise AcquisitionError(f"unexpected HTTP status {status_code}")
    headers: dict[str, str] = {}
    for line in lines[1:]:
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        headers[key.strip().lower()] = value.strip()
    return {
        "status_code": status_code,
        "etag": headers.get("etag"),
        "last_modified": headers.get("last-modified"),
        "content_length": headers.get("content-length"),
    }


def validate_remote_metadata(
    metadata: dict[str, Any],
    asset: dict[str, Any],
) -> dict[str, Any]:
    raw_content_length = metadata.get("content_length")
    content_length: int | None = None
    if raw_content_length is not None:
        try:
            content_length = int(raw_content_length)
        except ValueError as error:
            raise AcquisitionError("remote Content-Length is not an integer") from error
        if content_length != int(asset["size_bytes"]):
            raise AcquisitionError(
                f"remote Content-Length {content_length} != expected "
                f"{asset['size_bytes']}"
            )
    return {
        "etag": metadata.get("etag"),
        "last_modified": metadata.get("last_modified"),
        "content_length": content_length,
    }


def probe_remote_asset(
    asset: dict[str, Any],
    download_dir: Path,
) -> dict[str, Any]:
    descriptor, temporary_name = tempfile.mkstemp(
        dir=download_dir,
        prefix=".remote-headers.",
        suffix=".tmp",
    )
    os.close(descriptor)
    headers_path = Path(temporary_name)
    try:
        subprocess.run(curl_probe_command(asset, headers_path), check=True)
        return validate_remote_metadata(parse_final_http_headers(headers_path), asset)
    finally:
        headers_path.unlink(missing_ok=True)


def resume_state_path(partial_path: Path) -> Path:
    return partial_path.with_name(f"{partial_path.name}.resume.json")


def load_resume_state(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AcquisitionError(f"invalid resume state {path.name}: {error}") from error
    if payload.get("schema_version") != RESUME_SCHEMA_VERSION:
        raise AcquisitionError(f"unsupported resume state schema in {path.name}")
    return payload


def remote_validator(metadata: dict[str, Any]) -> str | None:
    etag = metadata.get("etag")
    if isinstance(etag, str) and etag and not etag.startswith("W/"):
        return etag
    last_modified = metadata.get("last_modified")
    if isinstance(last_modified, str) and last_modified:
        return last_modified
    return None


@contextmanager
def asset_lock(download_dir: Path, filename: str) -> Iterator[None]:
    lock_path = download_dir / f".{filename}.lock"
    with lock_path.open("a+", encoding="utf-8") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise AcquisitionError(
                f"{filename}: another acquisition is active"
            ) from error
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def reject_symlink(path: Path) -> None:
    if path.is_symlink():
        raise AcquisitionError(f"refusing symlinked acquisition path: {path}")


def download_asset(
    name: str,
    asset: dict[str, Any],
    download_dir: Path,
    *,
    manifest_sha256: str,
    dry_run: bool,
    probe: Callable[[dict[str, Any], Path], dict[str, Any]] = probe_remote_asset,
    prior_receipt_entry: dict[str, Any] | None = None,
) -> dict[str, Any]:
    destination = download_dir / str(asset["filename"])
    partial = destination.with_name(f"{destination.name}.part")
    resume_path = resume_state_path(partial)
    response_headers = destination.with_name(
        f".{destination.name}.response-headers.tmp"
    )
    expected_size = int(asset["size_bytes"])
    binding = asset_binding(name, asset, manifest_sha256)

    for path in (destination, partial, resume_path, response_headers):
        reject_symlink(path)
    if destination.exists():
        if partial.exists() or resume_path.exists():
            raise AcquisitionError(
                f"{name}: final file coexists with partial/resume state; "
                "refusing ambiguity"
            )
        verified = verify_local_asset(destination, asset)
        if asset.get("sha256") is None:
            if not receipt_entry_matches_binding(
                prior_receipt_entry,
                name,
                asset,
                manifest_sha256,
            ):
                raise AcquisitionError(
                    f"{name}: source-unpinned final has no receipt bound to the "
                    "current manifest; refusing reuse"
                )
            if prior_receipt_entry.get("sha256") != verified["sha256"]:
                raise AcquisitionError(
                    f"{name}: source-unpinned final differs from its receipt"
                )
        verified["status"] = (
            "reused_source_verified"
            if asset.get("sha256") is not None
            else "reused_local_hash_observed"
        )
        return verified

    if resume_path.exists() and not partial.exists():
        raise AcquisitionError(
            f"{name}: resume state exists without a partial file; refusing ambiguity"
        )
    if partial.exists() and partial.stat().st_size > expected_size:
        raise AcquisitionError(
            f"{name}: partial file exceeds expected size; remove it manually"
        )

    if dry_run:
        if_range = None
        if partial.exists():
            if not resume_path.is_file():
                raise AcquisitionError(
                    f"{name}: partial file has no authenticated resume state"
                )
            state = load_resume_state(resume_path)
            if state.get("binding") != binding:
                raise AcquisitionError(
                    f"{name}: partial resume state does not match current manifest"
                )
            if_range = remote_validator(state.get("remote", {}))
            if if_range is None:
                raise AcquisitionError(
                    f"{name}: partial file has no stable remote validator"
                )
        command = curl_command(
            asset,
            partial,
            response_headers,
            if_range=if_range,
        )
        return {
            "filename": destination.name,
            "partial_filename": partial.name,
            "resume_state_filename": resume_path.name,
            "expected_size_bytes": expected_size,
            "source_sha256_pinned": asset.get("sha256") is not None,
            "status": "dry_run",
            "command": command,
        }

    if shutil.which("curl") is None:
        raise AcquisitionError("curl is required but was not found")
    download_dir.mkdir(parents=True, exist_ok=True)
    with asset_lock(download_dir, destination.name):
        for path in (destination, partial, resume_path, response_headers):
            reject_symlink(path)
        if destination.exists():
            if partial.exists() or resume_path.exists():
                raise AcquisitionError(
                    f"{name}: final file coexists with partial/resume state; "
                    "refusing ambiguity"
                )
            verified = verify_local_asset(destination, asset)
            if asset.get("sha256") is None:
                if not receipt_entry_matches_binding(
                    prior_receipt_entry,
                    name,
                    asset,
                    manifest_sha256,
                ):
                    raise AcquisitionError(
                        f"{name}: source-unpinned final has no receipt bound to "
                        "the current manifest; refusing reuse"
                    )
                if prior_receipt_entry.get("sha256") != verified["sha256"]:
                    raise AcquisitionError(
                        f"{name}: source-unpinned final differs from its receipt"
                    )
            verified["status"] = (
                "reused_source_verified"
                if asset.get("sha256") is not None
                else "reused_local_hash_observed"
            )
            return verified
        if resume_path.exists() and not partial.exists():
            raise AcquisitionError(
                f"{name}: resume state exists without a partial file; "
                "refusing ambiguity"
            )
        if partial.exists() and partial.stat().st_size > expected_size:
            raise AcquisitionError(
                f"{name}: partial file exceeds expected size; remove it manually"
            )

        current_remote = probe(asset, download_dir)
        current_validator = remote_validator(current_remote)
        if partial.exists():
            if not resume_path.is_file():
                raise AcquisitionError(
                    f"{name}: partial file has no authenticated resume state"
                )
            resume_state = load_resume_state(resume_path)
            if resume_state.get("binding") != binding:
                raise AcquisitionError(
                    f"{name}: partial resume state does not match current manifest"
                )
            stored_remote = resume_state.get("remote")
            if not isinstance(stored_remote, dict):
                raise AcquisitionError(f"{name}: malformed remote resume metadata")
            stored_validator = remote_validator(stored_remote)
            if stored_validator is None or current_validator is None:
                raise AcquisitionError(
                    f"{name}: source exposes no stable validator; restart transfer"
                )
            if stored_remote != current_remote:
                raise AcquisitionError(
                    f"{name}: remote validators changed; restart transfer"
                )
        else:
            if resume_path.exists():
                raise AcquisitionError(
                    f"{name}: resume state exists without a partial file"
                )
            atomic_write_json(
                resume_path,
                {
                    "schema_version": RESUME_SCHEMA_VERSION,
                    "binding": binding,
                    "remote": current_remote,
                    "recorded_at_utc": datetime.now(UTC).isoformat(),
                },
            )

        if partial.exists() and partial.stat().st_size == expected_size:
            actual_hash = sha256_file(partial)
            expected_hash = asset.get("sha256")
            if expected_hash is not None and actual_hash != expected_hash:
                raise AcquisitionError(
                    f"{name}: completed partial has a SHA-256 mismatch"
                )
            os.replace(partial, destination)
            resume_path.unlink()
            return {
                "filename": destination.name,
                "size_bytes": expected_size,
                "sha256": actual_hash,
                "source_sha256_pinned": expected_hash is not None,
                "status": (
                    "recovered_source_verified_partial"
                    if expected_hash is not None
                    else "recovered_local_hash_observed_partial"
                ),
            }

        command = curl_command(
            asset,
            partial,
            response_headers,
            if_range=current_validator if partial.exists() else None,
        )
        try:
            subprocess.run(command, check=True)
            response_remote = parse_final_http_headers(response_headers)
            response_validator = remote_validator(response_remote)
            if (
                current_validator is not None
                and response_validator != current_validator
            ):
                raise AcquisitionError(
                    f"{name}: remote validator changed during transfer"
                )
        except BaseException:
            if not partial.exists():
                resume_path.unlink(missing_ok=True)
            raise
        finally:
            response_headers.unlink(missing_ok=True)

        actual_size = partial.stat().st_size
        if actual_size != expected_size:
            raise AcquisitionError(
                f"{name}: completed transfer has {actual_size} bytes; "
                f"expected {expected_size}; partial file was preserved"
            )

        actual_hash = sha256_file(partial)
        expected_hash = asset.get("sha256")
        if expected_hash is not None and actual_hash != expected_hash:
            raise AcquisitionError(
                f"{name}: SHA-256 mismatch; partial file was preserved"
            )

        os.replace(partial, destination)
        resume_path.unlink()
        return {
            "filename": destination.name,
            "size_bytes": actual_size,
            "sha256": actual_hash,
            "source_sha256_pinned": expected_hash is not None,
            "status": (
                "downloaded_source_verified"
                if expected_hash is not None
                else "downloaded_local_hash_observed"
            ),
        }


def run_download(
    group: str,
    *,
    manifest_path: Path = DEFAULT_MANIFEST,
    download_dir: Path = DEFAULT_DOWNLOADS,
    confirm_large_download: bool = False,
    dry_run: bool = False,
) -> dict[str, Any]:
    manifest_path = require_under_experiment(manifest_path)
    download_dir = require_under_experiment(download_dir)
    manifest = load_manifest(manifest_path)
    selected = select_assets(manifest, group, confirm_large_download)
    manifest_sha256 = sha256_file(manifest_path)
    if not dry_run:
        download_dir.mkdir(parents=True, exist_ok=True)

    def acquire_selected(receipt: dict[str, Any]) -> dict[str, Any]:
        results: dict[str, Any] = {}
        for name, asset in selected:
            prior_entry = receipt["assets"].get(name)
            result = download_asset(
                name,
                asset,
                download_dir,
                manifest_sha256=manifest_sha256,
                dry_run=dry_run,
                prior_receipt_entry=prior_entry,
            )
            result["binding"] = asset_binding(name, asset, manifest_sha256)
            results[name] = result
        return results

    if dry_run:
        receipt = (
            load_receipt(download_dir)
            if download_dir.exists()
            else {"schema_version": RECEIPT_SCHEMA_VERSION, "assets": {}}
        )
        results = acquire_selected(receipt)
    else:
        with asset_lock(download_dir, "download-receipt"):
            receipt = load_receipt(download_dir)
            results = acquire_selected(receipt)
            recorded_at = datetime.now(UTC).isoformat()
            retained_assets = {
                name: entry
                for name, entry in receipt["assets"].items()
                if name in manifest["assets"]
                and receipt_entry_matches_binding(
                    entry,
                    name,
                    manifest["assets"][name],
                    manifest_sha256,
                )
            }
            updated_receipt = {
                "schema_version": RECEIPT_SCHEMA_VERSION,
                "dataset": manifest["dataset"],
                "manifest_sha256": manifest_sha256,
                "updated_at_utc": recorded_at,
                "assets": retained_assets,
            }
            updated_receipt["assets"].update(results)
            atomic_write_json(download_dir / RECEIPT_NAME, updated_receipt)

    payload = {
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "dataset": manifest["dataset"],
        "manifest_sha256": manifest_sha256,
        "recorded_at_utc": datetime.now(UTC).isoformat(),
        "requested_group": group,
        "dry_run": dry_run,
        "assets": results,
    }
    if not dry_run:
        payload["recorded_at_utc"] = updated_receipt["updated_at_utc"]
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("group", choices=("metadata", "raw", "all"))
    parser.add_argument(
        "--confirm-large-download",
        action="store_true",
        help="required for any command that includes the 3.75 GB raw archive",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="validate the plan and print curl commands without writing or downloading",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = run_download(
            args.group,
            confirm_large_download=args.confirm_large_download,
            dry_run=args.dry_run,
        )
    except (AcquisitionError, subprocess.CalledProcessError) as error:
        print(f"acquisition failed: {error}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
