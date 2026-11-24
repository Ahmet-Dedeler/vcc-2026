#!/usr/bin/env python3
"""Fail-closed acquisition helper for GSE295214.

Metadata is the default and safe path. Multi-gigabyte assets require the exact
confirmation token frozen in SOURCE_MANIFEST.json.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import sys
import tempfile
import urllib.request
import urllib.error
import urllib.parse
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MANIFEST_PATH = HERE / "SOURCE_MANIFEST.json"
PSEUDOBULK_AMENDMENT_PATH = HERE / "PSEUDOBULK_ACQUISITION_AMENDMENT.json"
DATA_ROOT = ROOT / "data" / "external" / "gse295214"
DOWNLOAD_RECEIPT_PATH = HERE / "work" / "pseudobulk-download-receipt.json"


class AcquisitionError(RuntimeError):
    """The source contract or download authorization failed."""


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, request, fp, code, msg, headers, newurl):
        return None


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def load_manifest(path: Path = MANIFEST_PATH) -> dict[str, Any]:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1:
        raise AcquisitionError("unsupported manifest schema")
    if manifest.get("dataset") != "GSE295214":
        raise AcquisitionError("manifest dataset must be GSE295214")
    if manifest.get("scope") != "metadata_only":
        raise AcquisitionError("manifest scope must remain metadata_only")
    assets = manifest.get("assets")
    if not isinstance(assets, dict) or not assets:
        raise AcquisitionError("manifest assets are missing")
    for name, asset in assets.items():
        if asset.get("group") not in {"metadata", "derived_metadata", "large"}:
            raise AcquisitionError(f"{name}: invalid asset group")
        path_value = asset.get("local_path")
        if not isinstance(path_value, str):
            raise AcquisitionError(f"{name}: local_path is missing")
        destination = (ROOT / path_value).resolve()
        if not destination.is_relative_to(DATA_ROOT.resolve()):
            raise AcquisitionError(f"{name}: destination escapes the data root")
        size = asset.get("size_bytes")
        if not isinstance(size, int) or size <= 0:
            raise AcquisitionError(f"{name}: invalid size")
        expected_hash = asset.get("sha256")
        if asset["group"] != "large" and (
            not isinstance(expected_hash, str) or len(expected_hash) != 64
        ):
            raise AcquisitionError(f"{name}: metadata must have a SHA-256")
        if asset["group"] == "large" and asset.get("downloaded") is not False:
            raise AcquisitionError(f"{name}: large asset must remain marked unacquired")
    return manifest


def load_pseudobulk_amendment(
    manifest: dict[str, Any],
    path: Path = PSEUDOBULK_AMENDMENT_PATH,
) -> dict[str, Any]:
    amendment = json.loads(path.read_text(encoding="utf-8"))
    if amendment.get("schema_version") != 1:
        raise AcquisitionError("unsupported pseudobulk amendment schema")
    if amendment.get("experiment") != "115-gse295214-h1":
        raise AcquisitionError("pseudobulk amendment experiment changed")
    if amendment.get("dataset") != manifest["dataset"]:
        raise AcquisitionError("pseudobulk amendment dataset changed")
    parent = amendment.get("parent_manifest")
    if not isinstance(parent, dict):
        raise AcquisitionError("pseudobulk amendment parent binding is missing")
    if parent.get("path") != str(MANIFEST_PATH.relative_to(ROOT)):
        raise AcquisitionError("pseudobulk amendment parent path changed")
    if parent.get("sha256") != sha256_file(MANIFEST_PATH):
        raise AcquisitionError("pseudobulk amendment parent hash changed")
    authorization = amendment.get("authorization")
    if not isinstance(authorization, dict):
        raise AcquisitionError("pseudobulk amendment authorization is missing")
    if authorization.get("asset_name") != "public_pseudobulk_database":
        raise AcquisitionError("pseudobulk amendment asset changed")
    if authorization.get("download_public_pseudobulk_database") is not True:
        raise AcquisitionError("pseudobulk download is not authorized")
    if authorization.get("download_geo_raw_archive") is not False:
        raise AcquisitionError("raw GEO archive must remain unauthorized")
    token = authorization.get("confirmation_token")
    if not isinstance(token, str) or not token:
        raise AcquisitionError("pseudobulk confirmation token is invalid")

    asset = manifest["assets"]["public_pseudobulk_database"]
    binding = amendment.get("asset_binding")
    expected_binding = {
        "repository": "xingjiepan/SCMG_data",
        "repository_revision": manifest["official_sources"]["public_derived_data"][
            "revision"
        ],
        "path": "data/pseudo_bulk_perturbation_database.h5ad",
        "huggingface_file_oid": "5b49230eb742fb17a278f660ca97e6e4e3c9a71a",
        "huggingface_lfs_oid_sha256": asset["sha256"],
        "huggingface_file_introducing_commit": (
            "6c8ede815030c7672e960be42a36781b20bc830c"
        ),
        "size_bytes": asset["size_bytes"],
    }
    if binding != expected_binding:
        raise AcquisitionError("pseudobulk asset binding changed")
    if authorization.get("max_download_bytes") != asset["size_bytes"]:
        raise AcquisitionError("pseudobulk byte ceiling changed")
    return amendment


def authorize_asset(
    manifest: dict[str, Any],
    asset_name: str,
    confirmation: str | None,
    amendment: dict[str, Any] | None = None,
) -> dict[str, Any]:
    try:
        asset = manifest["assets"][asset_name]
    except KeyError as error:
        raise AcquisitionError(f"unknown asset: {asset_name}") from error
    group = asset["group"]
    if group == "derived_metadata":
        raise AcquisitionError("derived metadata is produced, not downloaded")
    if group == "large":
        if asset_name != "public_pseudobulk_database":
            raise AcquisitionError(f"{asset_name} is not authorized by the amendment")
        if amendment is None:
            amendment = load_pseudobulk_amendment(manifest)
        expected = amendment["authorization"]["confirmation_token"]
        if confirmation != expected:
            raise AcquisitionError(
                f"{asset_name} is a large download; pass the exact confirmation token"
            )
        if (
            amendment["authorization"]["download_public_pseudobulk_database"]
            is not True
        ):
            raise AcquisitionError(
                f"{asset_name} remains disabled by the acquisition amendment"
            )
    return asset


def verify_asset(path: Path, asset: dict[str, Any]) -> None:
    if path.is_symlink() or not path.is_file():
        raise AcquisitionError(f"missing or symlinked asset: {path}")
    if path.stat().st_size != asset["size_bytes"]:
        raise AcquisitionError(f"{path.name}: byte count mismatch")
    expected_hash = asset.get("sha256")
    if expected_hash and sha256_file(path) != expected_hash:
        raise AcquisitionError(f"{path.name}: SHA-256 mismatch")


def download_asset(asset: dict[str, Any]) -> Path:
    destination = (ROOT / asset["local_path"]).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        verify_asset(destination, asset)
        return destination
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".part",
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with (
            urllib.request.urlopen(asset["url"]) as source,
            temporary.open("wb") as sink,
        ):
            while block := source.read(8 << 20):
                sink.write(block)
        verify_asset(temporary, asset)
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return destination


def resolve_pinned_pseudobulk_url(
    asset: dict[str, Any],
    amendment: dict[str, Any],
) -> tuple[str, dict[str, str]]:
    opener = urllib.request.build_opener(NoRedirect)
    request = urllib.request.Request(asset["url"], method="HEAD")
    try:
        opener.open(request, timeout=60)
    except urllib.error.HTTPError as error:
        if error.code not in {301, 302, 303, 307, 308}:
            raise
        headers = {key.lower(): value for key, value in error.headers.items()}
        location = urllib.parse.urljoin(asset["url"], headers["location"])
    else:
        raise AcquisitionError("expected a pinned Hugging Face redirect")
    binding = amendment["asset_binding"]
    if headers.get("x-repo-commit") != binding["repository_revision"]:
        raise AcquisitionError("pseudobulk repository revision drift")
    if int(headers.get("x-linked-size", "0")) != binding["size_bytes"]:
        raise AcquisitionError("pseudobulk remote byte count drift")
    if (
        headers.get("x-linked-etag", "").strip('"')
        != binding["huggingface_lfs_oid_sha256"]
    ):
        raise AcquisitionError("pseudobulk remote content hash drift")
    return location, headers


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
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def write_download_receipt(
    destination: Path,
    asset: dict[str, Any],
    amendment: dict[str, Any],
    headers: dict[str, str],
) -> None:
    receipt = {
        "schema_version": 1,
        "experiment": "115-gse295214-h1",
        "dataset": "GSE295214",
        "status": "verified_public_pseudobulk_download",
        "completed_at_utc": datetime.datetime.now(datetime.UTC)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z"),
        "source_manifest_sha256": sha256_file(MANIFEST_PATH),
        "acquisition_amendment_sha256": sha256_file(PSEUDOBULK_AMENDMENT_PATH),
        "asset": {
            "path": str(destination.relative_to(ROOT)),
            "size_bytes": destination.stat().st_size,
            "sha256": sha256_file(destination),
            "repository_revision": headers["x-repo-commit"],
            "remote_lfs_sha256": headers["x-linked-etag"].strip('"'),
        },
        "authorization": amendment["authorization"],
    }
    atomic_json(DOWNLOAD_RECEIPT_PATH, receipt)


def download_pseudobulk(
    asset: dict[str, Any],
    amendment: dict[str, Any],
) -> Path:
    destination = (ROOT / asset["local_path"]).resolve()
    if destination.parent != DATA_ROOT.resolve():
        raise AcquisitionError("pseudobulk destination escapes the data root")
    destination.parent.mkdir(parents=True, exist_ok=True)
    signed_url, headers = resolve_pinned_pseudobulk_url(asset, amendment)
    if destination.exists():
        verify_asset(destination, asset)
        write_download_receipt(destination, asset, amendment, headers)
        return destination

    partial = destination.with_name(f".{destination.name}.part")
    lock = destination.with_name(f".{destination.name}.download.lock")
    if partial.is_symlink() or lock.is_symlink():
        raise AcquisitionError("unsafe pseudobulk partial or lock path")
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as error:
        raise AcquisitionError("pseudobulk download lock already exists") from error
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(f"pid={os.getpid()}\n")
        offset = partial.stat().st_size if partial.exists() else 0
        maximum = amendment["authorization"]["max_download_bytes"]
        if offset > maximum:
            raise AcquisitionError("pseudobulk partial exceeds the byte ceiling")
        if offset < asset["size_bytes"]:
            request = urllib.request.Request(
                signed_url,
                headers={"Range": f"bytes={offset}-"},
            )
            with (
                urllib.request.urlopen(request, timeout=120) as source,
                partial.open("ab") as sink,
            ):
                if source.status != 206:
                    raise AcquisitionError(
                        "pseudobulk server did not honor resume range"
                    )
                content_range = source.headers.get("Content-Range", "")
                if not content_range.startswith(f"bytes {offset}-"):
                    raise AcquisitionError("pseudobulk resume range changed")
                total = offset
                next_report = ((total // (256 << 20)) + 1) * (256 << 20)
                while block := source.read(8 << 20):
                    total += len(block)
                    if total > maximum:
                        raise AcquisitionError(
                            "pseudobulk download exceeded byte ceiling"
                        )
                    sink.write(block)
                    if total >= next_report:
                        print(
                            f"downloaded {total / (1 << 30):.2f} GiB",
                            file=sys.stderr,
                            flush=True,
                        )
                        next_report += 256 << 20
                sink.flush()
                os.fsync(sink.fileno())
        verify_asset(partial, asset)
        os.replace(partial, destination)
        destination.chmod(0o444)
        write_download_receipt(destination, asset, amendment, headers)
        return destination
    finally:
        lock.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--asset",
        required=True,
        help="manifest asset name; metadata assets are safe by default",
    )
    parser.add_argument("--confirm-large-download")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    manifest = load_manifest()
    amendment = None
    if args.asset == "public_pseudobulk_database":
        amendment = load_pseudobulk_amendment(manifest)
    asset = authorize_asset(
        manifest,
        args.asset,
        args.confirm_large_download,
        amendment=amendment,
    )
    if args.dry_run:
        print(json.dumps({"authorized": True, "asset": args.asset}, sort_keys=True))
        return
    if asset["group"] == "large":
        path = download_pseudobulk(asset, amendment)
    else:
        path = download_asset(asset)
    print(path)


if __name__ == "__main__":
    main()
