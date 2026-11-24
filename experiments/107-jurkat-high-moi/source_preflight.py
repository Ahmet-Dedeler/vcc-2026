#!/usr/bin/env python3
"""Generate the lightweight synthetic preflight receipt for Experiment 107."""

from __future__ import annotations

import hashlib
import io
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any

import numpy as np
import scipy


HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "preregistered_protocol.json"
DEFAULT_RECEIPT = HERE / "SOURCE_PREFLIGHT_RECEIPT.json"
TEST_MODULES = (
    "test_preprocess.py",
    "test_source_stage.py",
    "test_guide_calling.py",
    "test_demultiplex.py",
)
BOUND_IMPLEMENTATIONS = (
    "source_stage.py",
    "preprocess.py",
    "guide_calling.py",
    "demultiplex.py",
)


class PreflightError(RuntimeError):
    """The source preflight could not establish its frozen prerequisites."""


def sha256_file(path: Path, block_bytes: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_bytes), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
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


def frozen_runtime_versions() -> tuple[str, str]:
    try:
        protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
        numpy_version = protocol["cell_and_assignment_freeze"][
            "transcriptome_cell_call"
        ]["rng_contract"]["numpy_version"]
        scipy_version = protocol["cell_and_assignment_freeze"]["guide_assignment"][
            "scipy_version"
        ]
    except (OSError, json.JSONDecodeError, KeyError, TypeError) as error:
        raise PreflightError("cannot read frozen NumPy/SciPy versions") from error
    return str(numpy_version), str(scipy_version)


def build_suite() -> unittest.TestSuite:
    loader = unittest.defaultTestLoader
    suite = unittest.TestSuite()
    for filename in TEST_MODULES:
        suite.addTests(
            loader.discover(
                start_dir=str(HERE),
                pattern=filename,
                top_level_dir=str(HERE),
            )
        )
    return suite


def generate_receipt(output_path: Path = DEFAULT_RECEIPT) -> dict[str, Any]:
    frozen_numpy, frozen_scipy = frozen_runtime_versions()
    if np.__version__ != frozen_numpy:
        raise PreflightError(
            f"NumPy version changed: current {np.__version__}, frozen {frozen_numpy}"
        )
    if scipy.__version__ != frozen_scipy:
        raise PreflightError(
            f"SciPy version changed: current {scipy.__version__}, frozen {frozen_scipy}"
        )

    transcript = io.StringIO()
    result = unittest.TextTestRunner(stream=transcript, verbosity=2).run(
        build_suite()
    )
    transcript_bytes = transcript.getvalue().encode("utf-8")
    payload = {
        "schema_version": 1,
        "status": "pass" if result.wasSuccessful() else "fail",
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "python_version": sys.version.split()[0],
        "protocol_sha256": sha256_file(PROTOCOL),
        "preflight_implementation_sha256": sha256_file(Path(__file__)),
        "implementation_sha256": {
            filename: sha256_file(HERE / filename)
            for filename in BOUND_IMPLEMENTATIONS
        },
        "test_module_sha256": {
            filename: sha256_file(HERE / filename) for filename in TEST_MODULES
        },
        "tests_run": result.testsRun,
        "failures": len(result.failures),
        "errors": len(result.errors),
        "skipped": len(result.skipped),
        "transcript_sha256": hashlib.sha256(transcript_bytes).hexdigest(),
        "held_outcomes_opened": False,
    }
    if not result.wasSuccessful():
        raise PreflightError(transcript.getvalue())
    atomic_write_json(output_path, payload)
    return payload


def main() -> int:
    try:
        receipt = generate_receipt()
    except PreflightError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
