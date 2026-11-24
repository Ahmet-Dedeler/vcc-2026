#!/usr/bin/env python3
"""Create Experiment 130's one-way pre-value implementation freeze."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "IMPLEMENTATION_FREEZE.json"
FILES = (
    HERE / "README.md",
    HERE / "RESULTS.md",
    HERE / "preregistered_protocol.json",
    HERE / "CANDIDATE_BINDING_AMENDMENT.schema.json",
    HERE / "VALUE_ACCESS_BOUNDARY.v1.json",
    HERE / "SOURCE_METADATA_RECEIPT.v1.json",
    HERE / "run.py",
    HERE / "test_run.py",
    HERE / "freeze.py",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", type=int, required=True)
    parser.add_argument("--ruff-passed", action="store_true")
    parser.add_argument("--compile-passed", action="store_true")
    args = parser.parse_args()
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite {OUTPUT}")
    if args.tests <= 0 or not args.ruff_passed or not args.compile_passed:
        raise ValueError("freeze requires tests, Ruff, and compilation")
    if (HERE / "CANDIDATE_BINDING_AMENDMENT.json").exists():
        raise ValueError("candidate binding must be absent at pre-value freeze")
    record = {
        "schema_version": 1,
        "experiment": "130-telohaec-panel-confirmation",
        "created_at": "2026-08-30",
        "status": "frozen_before_any_telohaec_expression_or_candidate_value_access",
        "candidate_binding_present": False,
        "biological_value_access": False,
        "files": {str(path.relative_to(ROOT)): sha256_file(path) for path in FILES},
        "preflight": {
            "synthetic_tests": f"{args.tests}/{args.tests} passed",
            "ruff": "passed",
            "python_compile": "passed"
        },
        "authorization": {
            "metadata_only": True,
            "telohaec_X_or_layers": False,
            "candidate_arrays": False,
            "hcasmc": False,
            "dld1": False,
            "pan": False,
            "kolf": False,
            "h1_or_orion_held_outcomes": False,
            "challenge_truth_or_context_values": False,
            "cell_eval2": False,
            "counts_or_submission": False,
            "commit_or_push": False
        }
    }
    OUTPUT.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    OUTPUT.chmod(0o444)
    print(sha256_file(OUTPUT))


if __name__ == "__main__":
    main()
