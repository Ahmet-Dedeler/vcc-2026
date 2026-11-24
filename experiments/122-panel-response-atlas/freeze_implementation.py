#!/usr/bin/env python3
"""Freeze the Experiment 122 implementation before source effect loading."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "IMPLEMENTATION_FREEZE.json"
FILES = (
    "experiments/122-panel-response-atlas/freeze_implementation.py",
    "experiments/122-panel-response-atlas/freeze_plans.py",
    "experiments/122-panel-response-atlas/IMPLEMENTATION_AMENDMENT_001.json",
    "experiments/122-panel-response-atlas/run.py",
    "experiments/122-panel-response-atlas/test_run.py",
)


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite {OUTPUT}")
    files = {relative: sha256_file(ROOT / relative) for relative in FILES}
    payload = {
        "schema_version": 1,
        "experiment": "122-panel-response-atlas",
        "phase": "phase_1_source_only",
        "created_before_effect_matrix_loading": False,
        "scientific_spec_frozen_before_effect_matrix_loading": True,
        "implementation_amendment": {
            "path": "experiments/122-panel-response-atlas/IMPLEMENTATION_AMENDMENT_001.json",
            "sha256": sha256_file(HERE / "IMPLEMENTATION_AMENDMENT_001.json")
        },
        "preregistration_sha256": sha256_file(HERE / "PREREGISTRATION.sha256"),
        "plan_freeze_manifest_sha256": sha256_file(
            HERE / "work/preregistered-freeze-v1/freeze-manifest.json"
        ),
        "files": files,
    }
    temporary = OUTPUT.with_name(f".{OUTPUT.name}.tmp")
    with temporary.open("xb") as handle:
        handle.write((json.dumps(payload, indent=2, sort_keys=True) + "\n").encode())
        handle.flush()
        os.fsync(handle.fileno())
    temporary.rename(OUTPUT)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
