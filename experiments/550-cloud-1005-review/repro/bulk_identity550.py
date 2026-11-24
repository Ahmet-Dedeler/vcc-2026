#!/usr/bin/env python3
"""Compare per-(context, target_gene) summed counts of two prediction files.

Streams the CSR (or a small dense matrix) and never builds the dense genes-by-cells matrix.
Accepts an .h5ad or a .vcc (a tar whose only member is pred.h5ad.zst).

Exit 0 when every group sum matches. Exit 1 on the first mismatch, or 2 on a bad file.
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import h5py
import numpy as np


def _decode(x) -> str:
    if isinstance(x, bytes):
        return x.decode()
    return str(x)


def _obs_column(obs: h5py.Group, name: str) -> list[str]:
    node = obs[name]
    if isinstance(node, h5py.Group) and "codes" in node and "categories" in node:
        codes = node["codes"][:]
        cats = [_decode(c) for c in node["categories"][:]]
        return [cats[int(c)] if int(c) >= 0 else "" for c in codes]
    arr = node.asstr()[:] if getattr(node, "dtype", None) is not None and h5py.check_string_dtype(node.dtype) else node[:]
    return [_decode(x) for x in arr]


def _open_h5ad(path: Path):
    """Return (h5py file, cleanup path or None). A .vcc is decompressed to a temp h5ad."""
    if path.suffix == ".vcc":
        with tarfile.open(path) as archive:
            members = [m for m in archive.getmembers() if m.isfile() and Path(m.name).name == "pred.h5ad.zst"]
            if len(members) != 1:
                raise SystemExit(f"{path}: expected exactly one pred.h5ad.zst, found {len(members)}")
            raw = archive.extractfile(members[0]).read()
        if shutil.which("zstd") is None:
            raise SystemExit("reading a .vcc needs the zstd binary")
        proc = subprocess.run(["zstd", "-d", "-q", "-c"], input=raw, capture_output=True)
        if proc.returncode != 0 or not proc.stdout:
            raise SystemExit(f"zstd failed to decompress {path}")
        tmp = tempfile.NamedTemporaryFile(prefix="bulk550-", suffix=".h5ad", delete=False)
        tmp.write(proc.stdout)
        tmp.close()
        return h5py.File(tmp.name, "r"), Path(tmp.name)
    return h5py.File(path, "r"), None


def group_sums(path: Path) -> dict[tuple[str, str], np.ndarray]:
    handle, cleanup = _open_h5ad(path)
    try:
        ctx = _obs_column(handle["obs"], "context")
        tgt = _obs_column(handle["obs"], "target_gene")
        if len(ctx) != len(tgt):
            raise SystemExit(f"{path}: context/target_gene length mismatch")
        n = len(ctx)
        acc: dict[tuple[str, str], np.ndarray] = {}
        x = handle["X"]
        if isinstance(x, h5py.Group) and "indptr" in x:
            indptr = x["indptr"]
            indices = x["indices"]
            data = x["data"]
            shape = x.attrs.get("shape")
            n_genes = int(shape[1]) if shape is not None else 0
            # one row at a time; the accumulator is one gene vector per (context, target)
            for i in range(n):
                key = (ctx[i], tgt[i])
                start = int(indptr[i])
                stop = int(indptr[i + 1])
                if key not in acc:
                    acc[key] = np.zeros(n_genes, dtype=np.float64)
                if stop > start:
                    cols = indices[start:stop]
                    vals = data[start:stop]
                    vec = acc[key]
                    if vec.size <= int(cols.max()):
                        grown = np.zeros(int(cols.max()) + 1, dtype=np.float64)
                        grown[: vec.size] = vec
                        acc[key] = grown
                        vec = grown
                    np.add.at(vec, cols.astype(np.int64, copy=False), vals.astype(np.float64, copy=False))
        else:
            # dense, used by the synthetic test when a writer stored a small matrix
            matrix = x
            g = int(matrix.shape[1])
            for i in range(n):
                key = (ctx[i], tgt[i])
                if key not in acc:
                    acc[key] = np.zeros(g, dtype=np.float64)
                acc[key] += np.asarray(matrix[i], dtype=np.float64)
        return acc
    finally:
        handle.close()
        if cleanup is not None:
            cleanup.unlink(missing_ok=True)


def compare(left: Path, right: Path) -> int:
    a = group_sums(left)
    b = group_sums(right)
    keys = sorted(set(a) | set(b))
    bad = 0
    for key in keys:
        if key not in a or key not in b:
            print(f"MISSING {key} left={key in a} right={key in b}")
            bad += 1
            continue
        va, vb = a[key], b[key]
        if va.shape != vb.shape or not np.array_equal(va, vb):
            delta = None
            if va.shape == vb.shape:
                delta = float(np.max(np.abs(va - vb)))
            print(f"DIFF {key[0]} {key[1]} shape {va.shape} vs {vb.shape} max_abs {delta}")
            bad += 1
    if bad:
        print(f"FAIL {bad} group(s)")
        return 1
    digest = hashlib.sha256()
    for key in keys:
        digest.update(repr(key).encode())
        digest.update(a[key].tobytes())
    print(f"PASS {len(keys)} groups  sum-sha256 {digest.hexdigest()}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("left", type=Path)
    ap.add_argument("right", type=Path)
    a = ap.parse_args(argv)
    try:
        return compare(a.left, a.right)
    except SystemExit:
        raise
    except Exception as exc:
        print(f"ERROR {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
