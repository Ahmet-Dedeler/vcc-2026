#!/usr/bin/env python3
"""Contract tests for held-context pruned whitening operators."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "prepare_held_pruned_operators",
    HERE / "prepare_held_pruned_operators.py",
)
if SPEC is None or SPEC.loader is None:
    raise ImportError(HERE / "prepare_held_pruned_operators.py")
PREPARE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = PREPARE
SPEC.loader.exec_module(PREPARE)


def test_operator_contract_requires_prune_and_reference_axis() -> None:
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        operator = root / "operator.npz"
        reference = root / "reference.json"
        np.savez(
            operator,
            evals=np.ones(2, dtype=np.float64),
            evecs=np.eye(2, dtype=np.float32),
        )
        reference.write_text(
            json.dumps(
                {
                    "n_cols": 17_942,
                    "gamma": 0.65,
                    "lambda": None,
                    "wcols_sha256": "axis",
                }
            )
        )
        metadata = {
            "axis": "union",
            "donors": ["k562_gwps", "orion_hek293t"],
            "n_cols": 17_942,
            "gamma": 0.65,
            "lambda": None,
            "wcols_sha256": "axis",
            "basis_prune": "essential",
            "dropped_in_panel": [],
            "sha256": PREPARE.sha256_file(operator),
        }
        operator.with_suffix(".json").write_text(json.dumps(metadata))
        original_axis = PREPARE.expected_operator_axis
        PREPARE.expected_operator_axis = lambda donors: np.asarray(
            [3, 7], dtype=np.int64
        )
        axis_hash = PREPARE.hashlib.sha256(
            PREPARE.expected_operator_axis(()).tobytes()
        ).hexdigest()
        metadata["n_cols"] = 2
        metadata["wcols_sha256"] = axis_hash
        reference.write_text(
            json.dumps(
                {
                    "n_cols": 2,
                    "gamma": 0.65,
                    "lambda": None,
                    "wcols_sha256": axis_hash,
                }
            )
        )
        operator.with_suffix(".json").write_text(json.dumps(metadata))
        try:
            PREPARE.validate_operator(
                operator,
                donors=("k562_gwps", "orion_hek293t"),
                reference_metadata=reference,
            )
        finally:
            PREPARE.expected_operator_axis = original_axis
        metadata["wcols_sha256"] = "wrong"
        operator.with_suffix(".json").write_text(json.dumps(metadata))
        PREPARE.expected_operator_axis = lambda donors: np.asarray(
            [3, 7], dtype=np.int64
        )
        try:
            PREPARE.validate_operator(
                operator,
                donors=("k562_gwps", "orion_hek293t"),
                reference_metadata=reference,
            )
        except ValueError:
            pass
        else:
            raise AssertionError("operator with the wrong axis hash was accepted")
        finally:
            PREPARE.expected_operator_axis = original_axis


def main() -> None:
    tests = [
        value
        for name, value in sorted(globals().items())
        if name.startswith("test_") and callable(value)
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"{len(tests)} tests passed")


if __name__ == "__main__":
    main()
