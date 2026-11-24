#!/usr/bin/env python3
"""Fail-closed two-transfer acquisition for Experiment 116 Low1/Low2 LFCs."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import tempfile
import urllib.request
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
AMENDMENT_PATH = HERE / "ACQUISITION_AMENDMENT.v1.json"
FREEZE_PATH = HERE / "PREVALUE_FREEZE.v1.json"
RECEIPT_PATH = HERE / "ACQUISITION_RECEIPT.v1.json"
CONFIRMATION = "I_CONFIRM_EXPERIMENT_116_TWO_TRANSFER_ACQUISITION_V1"
ASSET_NAMES = ("low1_lfc_matrix", "low2_lfc_matrix")
BLOCK_BYTES = 8 << 20
USER_AGENT = "vcc-experiment-116-acquisition/1"


class AcquisitionError(RuntimeError):
    """The authorization, remote binding, or independent verification failed."""


def sha256_file(path: Path, block_bytes: int = BLOCK_BYTES) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
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
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def load_amendment(path: Path = AMENDMENT_PATH) -> dict[str, Any]:
    amendment = json.loads(path.read_text(encoding="utf-8"))
    if amendment.get("schema_version") != 1:
        raise AcquisitionError("unsupported acquisition amendment schema")
    if amendment.get("experiment") != "116-dld1-lowmoi":
        raise AcquisitionError("wrong experiment in acquisition amendment")
    if amendment.get("status") != (
        "authorized_pending_independent_two_transfer_verification"
    ):
        raise AcquisitionError("acquisition amendment status changed")

    base = amendment.get("base_source_manifest", {})
    base_path = ROOT / base.get("path", "")
    if not base_path.is_file() or sha256_file(base_path) != base.get("sha256"):
        raise AcquisitionError("base SOURCE_MANIFEST.json binding failed")

    authorization = amendment.get("authorization", {})
    if authorization.get("authorized_asset_names") != list(ASSET_NAMES):
        raise AcquisitionError("only the ordered Low1/Low2 allowlist is authorized")
    for forbidden in (
        "download_geo_raw_archive",
        "download_any_other_processed_asset",
        "open_vcc_held_perturbation_truth",
        "open_existing_outcomes",
        "model_fitting",
        "held_context_evaluation",
        "cell_eval2",
        "count_generation",
        "challenge_submission",
        "commit",
        "push",
    ):
        if authorization.get(forbidden) is not False:
            raise AcquisitionError(f"{forbidden} must remain false")

    destination_root = (ROOT / authorization["authorized_destination_root"]).resolve()
    expected_root = (ROOT / "data/external/dld1-lowmoi").resolve()
    if destination_root != expected_root:
        raise AcquisitionError("acquisition destination root changed")

    assets = amendment.get("assets", {})
    if tuple(assets) != ASSET_NAMES:
        raise AcquisitionError("asset set or order changed")
    for name, asset in assets.items():
        if set(asset) != {
            "filename",
            "destination",
            "url",
            "size_bytes",
            "last_modified_http",
            "content_type",
            "accept_ranges",
            "etag",
            "bounded_header_probe_bytes",
            "header_line_bytes",
            "header_sha256",
            "target_label_count",
            "target_labels_sha256",
        }:
            raise AcquisitionError(f"{name}: unexpected asset fields")
        destination = (ROOT / asset["destination"]).resolve()
        if destination.parent != expected_root:
            raise AcquisitionError(f"{name}: destination escapes authorized root")
        if destination.name != asset["filename"]:
            raise AcquisitionError(f"{name}: filename/destination mismatch")
        if not asset["url"].startswith(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE337nnn/GSE337988/suppl/"
        ):
            raise AcquisitionError(f"{name}: URL is outside the frozen GEO root")
        if not isinstance(asset["size_bytes"], int) or asset["size_bytes"] <= 0:
            raise AcquisitionError(f"{name}: invalid byte count")
        if asset["etag"] is not None:
            raise AcquisitionError(f"{name}: unexpected ETag binding")
        for digest_key in ("header_sha256", "target_labels_sha256"):
            if re.fullmatch(r"[0-9a-f]{64}", asset[digest_key]) is None:
                raise AcquisitionError(f"{name}: invalid {digest_key}")

    contract = amendment.get("acceptance_contract", {})
    required = {
        "required_independent_full_transfers": 2,
        "cross_transfer_byte_reuse": False,
        "require_exact_size_each_transfer": True,
        "require_identical_sha256": True,
        "require_byte_for_byte_equality": True,
        "require_valid_complete_gzip_stream_each_transfer": True,
        "require_live_head_binding_before_each_asset": True,
        "publish_sha256_only_after_all_checks_pass": True,
        "atomic_finalization": True,
        "mismatch_policy": (
            "fail_closed_without_publishing_or_renaming_a_final_asset"
        ),
    }
    for key, value in required.items():
        if contract.get(key) != value:
            raise AcquisitionError(f"acceptance contract drift for {key}")
    return amendment


def verify_prevalue_freeze(path: Path = FREEZE_PATH) -> dict[str, Any]:
    if not path.is_file():
        raise AcquisitionError("PREVALUE_FREEZE.v1.json is required")
    freeze = json.loads(path.read_text(encoding="utf-8"))
    if freeze.get("status") != "frozen_before_matrix_value_access":
        raise AcquisitionError("pre-value freeze status changed")
    if freeze.get("matrix_values_accessed") is not False:
        raise AcquisitionError("pre-value freeze must record no value access")
    for relative, digest in freeze.get("bound_files", {}).items():
        candidate = (ROOT / relative).resolve()
        if not candidate.is_relative_to(HERE.resolve()):
            raise AcquisitionError("pre-value freeze binds an out-of-scope file")
        if not candidate.is_file() or sha256_file(candidate) != digest:
            raise AcquisitionError(f"pre-value bound file drift: {relative}")
    return freeze


def head_binding(asset: dict[str, Any]) -> dict[str, Any]:
    request = urllib.request.Request(
        asset["url"],
        method="HEAD",
        headers={
            "User-Agent": USER_AGENT,
            "Accept-Encoding": "identity",
            "Cache-Control": "no-cache",
        },
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        headers = response.headers
        result = {
            "status": response.status,
            "content_length": int(headers["Content-Length"]),
            "last_modified": headers.get("Last-Modified"),
            "content_type": headers.get_content_type(),
            "accept_ranges": headers.get("Accept-Ranges"),
            "etag": headers.get("ETag"),
        }
    expected = {
        "status": 200,
        "content_length": asset["size_bytes"],
        "last_modified": asset["last_modified_http"],
        "content_type": asset["content_type"],
        "accept_ranges": asset["accept_ranges"],
        "etag": asset["etag"],
    }
    if result != expected:
        raise AcquisitionError(
            f"{asset['filename']}: live HEAD binding drift: {result!r}"
        )
    return result


def _validate_get_headers(
    response: Any,
    *,
    asset: dict[str, Any],
    offset: int,
) -> dict[str, Any]:
    remaining = asset["size_bytes"] - offset
    status = response.status
    headers = response.headers
    if offset == 0:
        if status != 200:
            raise AcquisitionError(
                f"{asset['filename']}: initial GET returned HTTP {status}"
            )
        content_range = None
    else:
        if status != 206:
            raise AcquisitionError(
                f"{asset['filename']}: resumed GET returned HTTP {status}"
            )
        content_range = headers.get("Content-Range")
        expected_range = (
            f"bytes {offset}-{asset['size_bytes'] - 1}/{asset['size_bytes']}"
        )
        if content_range != expected_range:
            raise AcquisitionError(
                f"{asset['filename']}: Content-Range drift: {content_range!r}"
            )
    declared = headers.get("Content-Length")
    if declared is None or int(declared) != remaining:
        raise AcquisitionError(f"{asset['filename']}: GET byte count drift")
    if headers.get("Last-Modified") != asset["last_modified_http"]:
        raise AcquisitionError(f"{asset['filename']}: GET Last-Modified drift")
    return {
        "status": status,
        "offset": offset,
        "content_length": int(declared),
        "content_range": content_range,
        "last_modified": headers.get("Last-Modified"),
        "content_type": headers.get_content_type(),
        "cache_control": headers.get("Cache-Control"),
    }


def download_resumable(
    asset: dict[str, Any],
    partial_path: Path,
    *,
    transfer_name: str,
) -> dict[str, Any]:
    partial_path.parent.mkdir(parents=True, exist_ok=True)
    if partial_path.is_symlink():
        raise AcquisitionError(f"{partial_path}: symlink partial is forbidden")
    offset = partial_path.stat().st_size if partial_path.exists() else 0
    if offset > asset["size_bytes"]:
        raise AcquisitionError(f"{partial_path}: partial exceeds expected size")
    if offset == asset["size_bytes"]:
        return {
            "transfer_name": transfer_name,
            "resumed_from_bytes": offset,
            "network_bytes_this_invocation": 0,
            "get": None,
        }

    headers = {
        "User-Agent": f"{USER_AGENT}-{transfer_name}",
        "Accept-Encoding": "identity",
        "Cache-Control": "no-cache",
    }
    if offset:
        headers["Range"] = f"bytes={offset}-"
    request = urllib.request.Request(asset["url"], headers=headers)
    network_bytes = 0
    with urllib.request.urlopen(request, timeout=120) as response:
        get_metadata = _validate_get_headers(
            response,
            asset=asset,
            offset=offset,
        )
        mode = "ab" if offset else "wb"
        with partial_path.open(mode) as handle:
            while block := response.read(BLOCK_BYTES):
                handle.write(block)
                network_bytes += len(block)
            handle.flush()
            os.fsync(handle.fileno())
    if partial_path.stat().st_size != asset["size_bytes"]:
        raise AcquisitionError(f"{partial_path}: completed byte count mismatch")
    return {
        "transfer_name": transfer_name,
        "resumed_from_bytes": offset,
        "network_bytes_this_invocation": network_bytes,
        "get": get_metadata,
    }


def verify_complete_gzip(path: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    uncompressed_bytes = 0
    try:
        with gzip.open(path, "rb") as handle:
            while block := handle.read(BLOCK_BYTES):
                digest.update(block)
                uncompressed_bytes += len(block)
    except (EOFError, OSError) as error:
        raise AcquisitionError(f"{path.name}: invalid complete gzip stream") from error
    return {
        "uncompressed_size_bytes": uncompressed_bytes,
        "uncompressed_sha256": digest.hexdigest(),
    }


def files_identical(left: Path, right: Path) -> bool:
    if left.stat().st_size != right.stat().st_size:
        return False
    with left.open("rb") as left_handle, right.open("rb") as right_handle:
        while True:
            left_block = left_handle.read(BLOCK_BYTES)
            right_block = right_handle.read(BLOCK_BYTES)
            if left_block != right_block:
                return False
            if not left_block:
                return True


def _existing_receipt_is_valid(
    receipt: dict[str, Any],
    amendment: dict[str, Any],
) -> bool:
    if receipt.get("status") != "accepted_after_two_independent_transfers":
        return False
    if receipt.get("acquisition_amendment_sha256") != sha256_file(AMENDMENT_PATH):
        return False
    recorded = receipt.get("assets", {})
    if tuple(recorded) != ASSET_NAMES:
        return False
    for name in ASSET_NAMES:
        asset = amendment["assets"][name]
        destination = ROOT / asset["destination"]
        item = recorded[name]
        if (
            not destination.is_file()
            or destination.stat().st_size != asset["size_bytes"]
            or sha256_file(destination) != item.get("sha256")
        ):
            return False
    return True


def acquire_one(name: str, asset: dict[str, Any]) -> dict[str, Any]:
    if name not in ASSET_NAMES:
        raise AcquisitionError(f"asset is not authorized: {name}")
    destination = (ROOT / asset["destination"]).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    primary = destination.parent / f".{destination.name}.primary.part"
    verification = destination.parent / f".{destination.name}.verification.part"
    if destination.exists():
        raise AcquisitionError(
            f"{destination}: final exists without a valid complete receipt"
        )

    live_head = head_binding(asset)
    primary_transfer = download_resumable(
        asset,
        primary,
        transfer_name="primary",
    )
    verification_transfer = download_resumable(
        asset,
        verification,
        transfer_name="verification",
    )

    primary_sha = sha256_file(primary)
    verification_sha = sha256_file(verification)
    primary_gzip = verify_complete_gzip(primary)
    verification_gzip = verify_complete_gzip(verification)
    if primary_sha != verification_sha:
        raise AcquisitionError(f"{asset['filename']}: transfer SHA-256 mismatch")
    if primary_gzip != verification_gzip:
        raise AcquisitionError(
            f"{asset['filename']}: uncompressed stream verification mismatch"
        )
    if not files_identical(primary, verification):
        raise AcquisitionError(f"{asset['filename']}: byte comparison mismatch")

    os.replace(primary, destination)
    directory_fd = os.open(destination.parent, os.O_RDONLY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
    verification.unlink()
    return {
        "filename": asset["filename"],
        "path": asset["destination"],
        "size_bytes": asset["size_bytes"],
        "sha256": primary_sha,
        "gzip": primary_gzip,
        "live_head": live_head,
        "primary_transfer": primary_transfer,
        "verification_transfer": verification_transfer,
        "independent_full_transfers": 2,
        "byte_for_byte_equal": True,
        "accepted": True,
    }


def acquire() -> dict[str, Any]:
    amendment = load_amendment()
    freeze = verify_prevalue_freeze()
    if RECEIPT_PATH.exists():
        receipt = json.loads(RECEIPT_PATH.read_text(encoding="utf-8"))
        if _existing_receipt_is_valid(receipt, amendment):
            return receipt
        raise AcquisitionError("existing acquisition receipt or final files drifted")

    accepted: dict[str, Any] = {}
    for name in ASSET_NAMES:
        accepted[name] = acquire_one(name, amendment["assets"][name])

    receipt = {
        "schema_version": 1,
        "experiment": amendment["experiment"],
        "dataset": amendment["dataset"],
        "accepted_at": "2026-08-30",
        "status": "accepted_after_two_independent_transfers",
        "acquisition_amendment_sha256": sha256_file(AMENDMENT_PATH),
        "prevalue_freeze_sha256": sha256_file(FREEZE_PATH),
        "prevalue_bound_file_count": len(freeze["bound_files"]),
        "assets": accepted,
        "raw_archive_downloaded": False,
        "other_processed_asset_downloaded": False,
        "matrix_values_parsed": False,
    }
    atomic_json(RECEIPT_PATH, receipt)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("authenticate")
    acquire_parser = subparsers.add_parser("acquire")
    acquire_parser.add_argument("--confirm", required=True)
    args = parser.parse_args()

    if args.command == "authenticate":
        amendment = load_amendment()
        freeze = verify_prevalue_freeze()
        print(
            json.dumps(
                {
                    "status": "authorized",
                    "assets": list(amendment["assets"]),
                    "prevalue_freeze_sha256": sha256_file(FREEZE_PATH),
                    "bound_file_count": len(freeze["bound_files"]),
                },
                indent=2,
                sort_keys=True,
            )
        )
        return
    if args.confirm != CONFIRMATION:
        raise AcquisitionError("explicit acquisition confirmation token is required")
    print(json.dumps(acquire(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
