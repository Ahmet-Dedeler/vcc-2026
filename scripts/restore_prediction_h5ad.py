#!/usr/bin/env python3
"""Recover prediction content from a .vcc, not the original HDF5 byte layout.

Default: restore cell_0..cell_N-1 row names. --raw preserves the payload unchanged.
Existing outputs are never overwritten, and failures never publish a partial output.
The optional digest uses the September 8 storage manifest's content-hash algorithm
(CSR shape/counts/indices, target/context labels, gene order; row names validated
modulo cell_). It is NOT the novelty ledger's differently defined content digest.
"""
from __future__ import annotations

import argparse
from contextlib import suppress
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import tempfile

import h5py
import numpy as np


def strings(group: h5py.Group, name: str) -> list[str]:
    """Read the three string encodings used by our historical prediction files."""
    obj = group[name]
    if isinstance(obj, h5py.Group):
        if "values" in obj:
            if "mask" in obj and np.any(obj["mask"][...]):
                raise ValueError(f"{obj.name}: missing labels are not recoverable")
            values = obj["values"][...]
        elif "categories" in obj and "codes" in obj:
            categories = strings(obj, "categories")
            codes = obj["codes"][...]
            if np.any(codes < 0) or np.any(codes >= len(categories)):
                raise ValueError(f"{obj.name}: invalid categorical codes")
            return [categories[int(code)] for code in codes]
        else:
            raise ValueError(f"{obj.name}: unsupported string encoding")
    else:
        values = obj[...]
    return [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values]


def validate_index(f: h5py.File) -> list[str]:
    labels = strings(f["obs"], "_index")
    n = int(f["X"].attrs["shape"][0])
    if len(labels) != n:
        raise ValueError("obs index length does not match X")
    # Check EVERY row, not merely the first 1,000, and do not silently erase names.
    for i, label in enumerate(labels):
        if label not in (str(i), f"cell_{i}"):
            raise ValueError(f"unexpected obs index at row {i}: {label!r}")
    return labels


def storage_content_sha256(path: Path) -> str:
    """Reproduce lossless.tsv's original algorithm, recovered from the sweep code.

    The np.int64 tuple representation in the shape prefix is intentional: the
    original NumPy 2.x implementation hashed str(tuple(X.attrs['shape'])).
    """
    h = hashlib.sha256()
    with h5py.File(path, "r") as f:
        validate_index(f)
        shape = tuple(np.int64(v) for v in f["X"].attrs["shape"])
        h.update(str(shape).encode())
        for name in ("data", "indices", "indptr"):
            ds = f[f"X/{name}"]
            dtype = np.float32 if name == "data" else np.int64
            for start in range(0, len(ds), 4_000_000):
                h.update(np.asarray(ds[start:start + 4_000_000], dtype=dtype).tobytes())
        for name in ("target_gene", "context"):
            h.update("\0".join(strings(f["obs"], name)).encode())
        h.update("\0".join(strings(f["var"], "_index")).encode())
    return h.hexdigest()


def restore(package: Path, output: Path, *, raw: bool = False,
            expected_content: str | None = None) -> dict:
    if expected_content is not None and not re.fullmatch(r"[0-9a-f]{64}", expected_content):
        raise ValueError("expected content digest must be 64 lowercase hex characters")
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"refusing to overwrite {output}")
    if not output.parent.is_dir():
        raise FileNotFoundError(f"output directory does not exist: {output.parent}")
    if shutil.which("zstd") is None:
        raise RuntimeError("zstd is required")
    fd, name = tempfile.mkstemp(prefix=".restore-", suffix=".h5ad", dir=output.parent)
    tmp = Path(name)
    try:
        with os.fdopen(fd, "wb") as out, tarfile.open(package, "r:") as archive:
            members = [m for m in archive if m.name == "pred.h5ad.zst"]
            if len(members) != 1 or not members[0].isfile():
                raise ValueError("package must contain exactly one regular pred.h5ad.zst")
            with archive.extractfile(members[0]) as source:
                # No shell interpolation; quotes and metacharacters in paths are literal.
                proc = subprocess.Popen(["zstd", "-d", "-q", "-c"],
                                        stdin=subprocess.PIPE, stdout=out)
                try:
                    shutil.copyfileobj(source, proc.stdin, length=1 << 20)
                    proc.stdin.close()
                    if proc.wait() != 0:
                        raise RuntimeError("zstd failed; output was not published")
                except BaseException:
                    # A decompressor that exited early can break the pipe on close
                    # too; cleanup must still reap our child and discard scratch.
                    with suppress(OSError):
                        proc.stdin.close()
                    if proc.poll() is None:
                        proc.kill()
                    proc.wait()
                    raise
            out.flush()
            os.fsync(out.fileno())
        with h5py.File(tmp, "r" if raw else "r+") as f:
            labels = validate_index(f)
            if not raw and any(label != f"cell_{i}" for i, label in enumerate(labels)):
                del f["obs"]["_index"]
                ds = f["obs"].create_dataset(
                    "_index", data=[f"cell_{i}" for i in range(len(labels))],
                    dtype=h5py.string_dtype())
                ds.attrs["encoding-type"] = "string-array"
                ds.attrs["encoding-version"] = "0.2.0"
        digest = storage_content_sha256(tmp) if expected_content is not None else None
        if expected_content is not None and digest != expected_content:
            raise ValueError(f"content mismatch: expected {expected_content}, got {digest}")
        # Same-filesystem hard link is atomic and fails if another writer created output.
        os.link(tmp, output)
        return {"package": str(package), "output": str(output), "raw": raw,
                "storage_content_sha256": digest,
                "original_file_bytes_reproduced": "not claimed"}
    finally:
        tmp.unlink(missing_ok=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("vcc", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--raw", action="store_true")
    ap.add_argument("--expect-content-sha256")
    args = ap.parse_args()
    try:
        result = restore(args.vcc, args.out, raw=args.raw,
                         expected_content=args.expect_content_sha256)
    except (OSError, ValueError, RuntimeError, tarfile.TarError) as exc:
        ap.exit(1, f"restore failed: {exc}\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
