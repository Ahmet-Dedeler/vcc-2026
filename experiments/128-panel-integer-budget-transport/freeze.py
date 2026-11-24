#!/usr/bin/env python3
"""Create Experiment 128's one-way implementation freeze."""

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
    HERE / "preregistered_protocol.json",
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
        raise ValueError("freeze requires passing tests, Ruff, and compilation")
    record = {
        "schema_version": 1,
        "experiment": "128-panel-integer-budget-transport",
        "created_at_utc": "2026-08-30",
        "status": "frozen_before_experiment_128_released_control_value_access",
        "files": {
            str(path.relative_to(ROOT)): sha256_file(path)
            for path in FILES
        },
        "preflight": {
            "synthetic_tests": f"{args.tests}/{args.tests} passed",
            "ruff": "passed",
            "python_compile": "passed",
        },
        "authorization": {
            "released_control_values_after_this_freeze": True,
            "experiment_122_frozen_source_result": True,
            "staged_a15_effects": True,
            "held_outcomes": False,
            "challenge_perturbation_truth": False,
            "cell_eval2": False,
            "count_h5ad_or_vcc_build": False,
            "submission": False,
            "commit": False,
            "push": False,
        },
    }
    OUTPUT.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    OUTPUT.chmod(0o444)
    print(sha256_file(OUTPUT))


if __name__ == "__main__":
    main()
