#!/usr/bin/env python3
"""Extract only the Pan2025_hESC target axis from public H5AD obs metadata.

The remote file is multi-gigabyte, but this reader uses bounded HTTP ranges and
never indexes X or layers. It exists to make the 2,354-row QC axis reproducible
without downloading the full public pseudobulk database.
"""

from __future__ import annotations

import collections
import hashlib
import io
import json
import os
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from acquire import DATA_ROOT, HERE, MANIFEST_PATH, ROOT, load_manifest, sha256_file


RECEIPT_PATH = HERE / "POST_QC_METADATA_RECEIPT.json"


class MetadataExtractionError(RuntimeError):
    """The bounded metadata extraction contract failed."""


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, request, fp, code, msg, headers, newurl):
        return None


def redirect_metadata(url: str) -> tuple[dict[str, str], str]:
    opener = urllib.request.build_opener(NoRedirect)
    request = urllib.request.Request(url, method="HEAD")
    try:
        opener.open(request)
    except urllib.error.HTTPError as error:
        if error.code not in {301, 302, 303, 307, 308}:
            raise
        headers = {key.lower(): value for key, value in error.headers.items()}
        location = urllib.parse.urljoin(url, headers["location"])
        return headers, location
    raise MetadataExtractionError("expected a pinned Hugging Face redirect")


class BoundedRangeFile(io.RawIOBase):
    def __init__(self, url: str, size: int, byte_cap: int, block_size: int = 1 << 20):
        self.url = url
        self.size = size
        self.byte_cap = byte_cap
        self.block_size = block_size
        self.position = 0
        self.bytes_fetched = 0
        self.cache: collections.OrderedDict[int, bytes] = collections.OrderedDict()

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.position

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        if whence == io.SEEK_SET:
            position = offset
        elif whence == io.SEEK_CUR:
            position = self.position + offset
        elif whence == io.SEEK_END:
            position = self.size + offset
        else:
            raise ValueError("invalid whence")
        if position < 0:
            raise ValueError("negative seek")
        self.position = position
        return position

    def _block(self, index: int) -> bytes:
        if index in self.cache:
            self.cache.move_to_end(index)
            return self.cache[index]
        start = index * self.block_size
        end = min(self.size, start + self.block_size) - 1
        requested = end - start + 1
        if self.bytes_fetched + requested > self.byte_cap:
            raise MetadataExtractionError("remote metadata byte cap exceeded")
        request = urllib.request.Request(
            self.url,
            headers={"Range": f"bytes={start}-{end}"},
        )
        with urllib.request.urlopen(request) as response:
            payload = response.read()
            status = response.status
            content_range = response.headers.get("Content-Range", "")
        if status != 206 or not content_range.startswith(f"bytes {start}-{end}/"):
            raise MetadataExtractionError("remote server did not honor the exact range")
        if len(payload) != requested:
            raise MetadataExtractionError("remote range length mismatch")
        self.bytes_fetched += len(payload)
        self.cache[index] = payload
        while len(self.cache) > 64:
            self.cache.popitem(last=False)
        return payload

    def readinto(self, buffer: bytearray) -> int:
        if self.position >= self.size:
            return 0
        requested = min(len(buffer), self.size - self.position)
        output = memoryview(buffer)
        written = 0
        while written < requested:
            index = self.position // self.block_size
            payload = self._block(index)
            offset = self.position % self.block_size
            take = min(requested - written, len(payload) - offset)
            output[written : written + take] = payload[offset : offset + take]
            written += take
            self.position += take
        return written


def decode_categorical(group: h5py.Group) -> np.ndarray:
    if set(group.keys()) != {"categories", "codes"}:
        raise MetadataExtractionError("unexpected categorical encoding")
    categories = np.array(
        [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in group["categories"][...]
        ],
        dtype=object,
    )
    codes = group["codes"][...]
    if np.any(codes < 0) or np.any(codes >= len(categories)):
        raise MetadataExtractionError("invalid categorical code")
    return categories[codes]


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def extract() -> dict[str, Any]:
    manifest = load_manifest()
    asset = manifest["assets"]["public_pseudobulk_database"]
    headers, signed_url = redirect_metadata(asset["url"])
    if (
        headers.get("x-repo-commit")
        != manifest["official_sources"]["public_derived_data"]["revision"]
    ):
        raise MetadataExtractionError("public data repository revision drift")
    if int(headers.get("x-linked-size", "0")) != asset["size_bytes"]:
        raise MetadataExtractionError("public pseudobulk byte count drift")
    linked_etag = headers.get("x-linked-etag", "").strip('"')
    if linked_etag != asset["sha256"]:
        raise MetadataExtractionError("public pseudobulk content hash drift")

    remote = BoundedRangeFile(
        signed_url,
        size=asset["size_bytes"],
        byte_cap=asset["range_byte_cap"],
    )
    with h5py.File(remote, "r") as handle:
        if "obs" not in handle:
            raise MetadataExtractionError("remote H5AD has no obs group")
        obs = handle["obs"]
        if "dataset" not in obs or "perturbed_gene_name" not in obs:
            raise MetadataExtractionError("required obs metadata is missing")
        datasets = decode_categorical(obs["dataset"]).astype(str)
        target_names = decode_categorical(obs["perturbed_gene_name"]).astype(str)
        if datasets.shape != target_names.shape:
            raise MetadataExtractionError("obs metadata axes disagree")
        targets = sorted(set(target_names[datasets == "Pan2025_hESC"]))

    expected = manifest["assets"]["post_qc_target_axis"]
    if len(targets) != expected["row_count"]:
        raise MetadataExtractionError("post-QC target count drift")
    payload = ("\n".join(targets) + "\n").encode()
    digest = hashlib.sha256(payload).hexdigest()
    if digest != expected["sha256"] or len(payload) != expected["size_bytes"]:
        raise MetadataExtractionError("post-QC target-axis content drift")
    output = DATA_ROOT / expected["filename"]
    if output.resolve().parent != DATA_ROOT.resolve():
        raise MetadataExtractionError("target output escapes the metadata directory")
    atomic_write(output, payload)

    receipt = {
        "schema_version": 1,
        "dataset": "GSE295214",
        "manifest_sha256": sha256_file(MANIFEST_PATH),
        "source": {
            "repository_revision": headers["x-repo-commit"],
            "file_size_bytes": asset["size_bytes"],
            "file_sha256": linked_etag,
            "access_mode": "bounded_http_range_h5ad_obs_only",
            "allowed_h5_paths": ["/obs/dataset", "/obs/perturbed_gene_name"],
            "forbidden_h5_paths": ["/X", "/layers"],
            "byte_cap": asset["range_byte_cap"],
            "bytes_fetched": remote.bytes_fetched,
        },
        "output": {
            "path": str(output.relative_to(ROOT)),
            "row_count": len(targets),
            "size_bytes": len(payload),
            "sha256": digest,
        },
    }
    atomic_write(
        RECEIPT_PATH,
        (json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode(),
    )
    return receipt


if __name__ == "__main__":
    print(json.dumps(extract(), indent=2, sort_keys=True))
