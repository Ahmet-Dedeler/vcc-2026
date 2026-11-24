#!/usr/bin/env python3
"""Verify and safely extract the official STATE VCC support archive."""

from __future__ import annotations

import argparse
import base64
import hashlib
import json
import os
import shutil
import stat
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from zipfile import BadZipFile, ZipFile

import google_crc32c


EXPECTED_SIZE = 8_716_992_349
EXPECTED_CRC32C_BASE64 = "MBP7uA=="
EXPECTED_GENERATION = "1763329752965344"


def digest_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> tuple[str, str, int]:
    sha256 = hashlib.sha256()
    crc32c = google_crc32c.Checksum()
    size = 0
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            size += len(chunk)
            sha256.update(chunk)
            crc32c.update(chunk)
    return sha256.hexdigest(), base64.b64encode(crc32c.digest()).decode("ascii"), size


def validate_members(archive: ZipFile, destination: Path) -> list[dict[str, object]]:
    root = destination.resolve()
    members: list[dict[str, object]] = []
    for info in archive.infolist():
        member_path = (root / info.filename).resolve()
        if os.path.commonpath((str(root), str(member_path))) != str(root):
            raise ValueError(f"unsafe ZIP member path: {info.filename!r}")

        unix_mode = info.external_attr >> 16
        if stat.S_ISLNK(unix_mode):
            raise ValueError(f"symbolic links are not allowed in support archive: {info.filename!r}")

        members.append(
            {
                "name": info.filename,
                "is_dir": info.is_dir(),
                "compressed_bytes": info.compress_size,
                "uncompressed_bytes": info.file_size,
                "crc32": f"{info.CRC:08x}",
            }
        )
    return members


def extracted_payload_root(staging: Path, destination: Path) -> Path:
    """Return the directory that should become destination after extraction.

    The official archive is wrapped in a directory named
    ``competition_support_set``. Avoid preserving that wrapper twice when the
    requested destination has the same name.
    """
    entries = list(staging.iterdir())
    if len(entries) == 1 and entries[0].is_dir() and entries[0].name == destination.name:
        return entries[0]
    return staging


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()

    archive_path = args.archive.resolve()
    destination = args.destination.resolve()
    manifest_path = args.manifest.resolve()

    sha256, crc32c_b64, size = digest_file(archive_path)
    if size != EXPECTED_SIZE:
        raise SystemExit(f"size mismatch: observed {size}, expected {EXPECTED_SIZE}")
    if crc32c_b64 != EXPECTED_CRC32C_BASE64:
        raise SystemExit(
            f"CRC32C mismatch: observed {crc32c_b64}, expected {EXPECTED_CRC32C_BASE64}"
        )

    try:
        with ZipFile(archive_path) as archive:
            members = validate_members(archive, destination)
    except BadZipFile as exc:
        raise SystemExit(f"invalid ZIP archive: {exc}") from exc

    uncompressed_bytes = sum(int(member["uncompressed_bytes"]) for member in members)
    free_bytes = shutil.disk_usage(destination.parent).free
    required_free_bytes = int(uncompressed_bytes * 1.20)
    if not args.verify_only and free_bytes < required_free_bytes:
        raise SystemExit(
            f"insufficient disk space: {free_bytes} free, "
            f"{required_free_bytes} required for extraction headroom"
        )

    manifest = {
        "verified_at": datetime.now(timezone.utc).isoformat(),
        "archive": {
            "path": str(archive_path),
            "size_bytes": size,
            "sha256": sha256,
            "crc32c_base64": crc32c_b64,
            "gcs_generation": EXPECTED_GENERATION,
        },
        "zip": {
            "member_count": len(members),
            "compressed_bytes_total": sum(int(member["compressed_bytes"]) for member in members),
            "uncompressed_bytes_total": uncompressed_bytes,
            "members": members,
        },
        "extraction": {
            "requested": not args.verify_only,
            "destination": str(destination),
            "available_bytes_before": free_bytes,
            "required_headroom_bytes": required_free_bytes,
        },
    }

    if not args.verify_only:
        if destination.exists():
            raise SystemExit(f"destination already exists: {destination}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(
            tempfile.mkdtemp(prefix=f".{destination.name}.extracting-", dir=destination.parent)
        )
        try:
            with ZipFile(archive_path) as archive:
                archive.extractall(staging)
            payload_root = extracted_payload_root(staging, destination)
            payload_root.rename(destination)
            if staging.exists():
                staging.rmdir()
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        manifest["extraction"]["completed"] = True

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest["archive"], sort_keys=True))
    print(
        json.dumps(
            {
                "member_count": len(members),
                "uncompressed_bytes_total": uncompressed_bytes,
                "extracted": not args.verify_only,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
