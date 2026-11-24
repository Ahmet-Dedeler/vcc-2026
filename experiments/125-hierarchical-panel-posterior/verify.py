#!/usr/bin/env python3
"""Independently verify Experiment 125 immutable result artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "work" / "source-only-evaluation-v1"


class VerificationError(RuntimeError):
    """An immutable result verification failed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_array(array: np.ndarray) -> str:
    value = np.ascontiguousarray(array)
    digest = hashlib.sha256()
    digest.update(str(value.dtype).encode())
    digest.update(b"\0")
    digest.update(json.dumps(list(value.shape), separators=(",", ":")).encode())
    digest.update(b"\0")
    digest.update(value.tobytes(order="C"))
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"JSON root must be object: {path}")
    return value


def main() -> None:
    manifest = read_json(OUTPUT / "artifact-manifest.json")
    receipt = read_json(OUTPUT / "result-receipt.json")
    results = read_json(OUTPUT / "results.json")
    require(
        receipt["artifact_manifest_sha256"]
        == sha256_file(OUTPUT / "artifact-manifest.json"),
        "receipt manifest hash mismatch",
    )
    require(
        receipt["results_sha256"] == sha256_file(OUTPUT / "results.json"),
        "receipt results hash mismatch",
    )
    for name, record in manifest["files"].items():
        path = OUTPUT / name
        require(path.is_file() and not path.is_symlink(), f"missing artifact: {name}")
        require(path.stat().st_size == record["size_bytes"], f"size mismatch: {name}")
        require(sha256_file(path) == record["sha256"], f"file hash mismatch: {name}")
        if record["logical_array_sha256"] is not None:
            array = np.load(path, allow_pickle=False)
            require(
                sha256_array(array) == record["logical_array_sha256"],
                f"logical array hash mismatch: {name}",
            )
    candidate = np.load(OUTPUT / "observed-candidate-effect.npy", allow_pickle=False)
    baseline = np.load(OUTPUT / "observed-baseline-effect.npy", allow_pickle=False)
    selected = np.load(OUTPUT / "selected-candidate-index.npy", allow_pickle=False)
    require(candidate.shape == (3, 300, 300), "candidate shape changed")
    require(baseline.shape == candidate.shape, "baseline shape changed")
    require(selected.shape == (3, 5), "selection shape changed")
    diagonal = np.arange(300)
    require(np.all(candidate[:, diagonal, diagonal] == 0.0), "candidate self effect")
    require(np.all(baseline[:, diagonal, diagonal] == 0.0), "baseline self effect")
    require(np.all(np.isfinite(candidate)), "candidate nonfinite")
    require(
        results["invariants"]["transport_rows_passed"]
        == results["invariants"]["transport_rows_total"],
        "transport invariant count mismatch",
    )
    require(results["invariants"]["expanded_pds"]["all_pass"], "PDS invariant failed")
    require(receipt["decision"] == results["decision"], "decision mismatch")
    summary = {
        "verified": True,
        "decision": results["decision"],
        "artifact_manifest_sha256": sha256_file(
            OUTPUT / "artifact-manifest.json"
        ),
        "result_receipt_sha256": sha256_file(OUTPUT / "result-receipt.json"),
        "results_sha256": sha256_file(OUTPUT / "results.json"),
        "files_verified": len(manifest["files"]),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
