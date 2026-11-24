#!/usr/bin/env python3
"""Freeze Experiment 125 implementation before source effect loading."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "IMPLEMENTATION_FREEZE.json"
FILES = (
    "README.md",
    "preregistered_protocol.json",
    "SOURCE_ONLY_AUTHORIZATION.json",
    "PREREGISTRATION.sha256",
    "freeze_plans.py",
    "run.py",
    "test_run.py",
    "verify.py",
    "freeze_implementation.py",
    "IMPLEMENTATION_FREEZE.v1.json",
    "IMPLEMENTATION_FAILURE_RECOVERY.v1.json",
    "IMPLEMENTATION_FREEZE.v2.json",
    "IMPLEMENTATION_FAILURE_RECOVERY.v2.json",
)


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUTPUT.exists():
        raise RuntimeError(f"refusing to overwrite {OUTPUT}")
    files = {}
    for name in FILES:
        path = HERE / name
        if not path.is_file() or path.is_symlink():
            raise RuntimeError(f"missing implementation file: {path}")
        files[str(path.relative_to(ROOT))] = sha256_file(path)
    manifest = {
        "schema_version": 1,
        "experiment": "125-hierarchical-panel-posterior",
        "phase": "source_only_nested_evaluation",
        "created_before_source_effect_loading": False,
        "scientific_spec_frozen_before_source_effect_loading": True,
        "implementation_recovery": {
            "path": "experiments/125-hierarchical-panel-posterior/IMPLEMENTATION_FAILURE_RECOVERY.v1.json",
            "sha256": sha256_file(HERE / "IMPLEMENTATION_FAILURE_RECOVERY.v1.json"),
        },
        "plan_freeze_manifest_sha256": sha256_file(
            HERE / "work/preregistered-freeze-v1/freeze-manifest.json"
        ),
        "files": files,
    }
    OUTPUT.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.sync()
    OUTPUT.chmod(0o444)
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
