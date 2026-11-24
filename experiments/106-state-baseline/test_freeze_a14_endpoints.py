#!/usr/bin/env python3
"""Static contract tests for the a14 held-context exporter."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "freeze_a14_endpoints", HERE / "freeze_a14_endpoints.py"
)
if SPEC is None or SPEC.loader is None:
    raise ImportError(HERE / "freeze_a14_endpoints.py")
FREEZE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = FREEZE
SPEC.loader.exec_module(FREEZE)


def test_a14_endpoint_constants_are_exact() -> None:
    assert FREEZE.A14_ALPHA == 0.982
    assert FREEZE.A14_CD4_FRACTION == 0.05


def test_exporter_has_no_truth_argument() -> None:
    import inspect

    parameters = inspect.signature(FREEZE.freeze).parameters
    assert {"truth", "reference", "real_effect"}.isdisjoint(parameters)
    assert "controls_h5ad" in parameters


def test_all_transitive_helper_sources_are_pinned() -> None:
    expected = {
        "experiments/104-context-heldout-credit-adapter/run.py",
        "experiments/100-sign-confidence/run.py",
        "experiments/062-cd4-post-whiten/run.py",
        "experiments/064-pm-cd4-composition/run.py",
    }
    assert set(FREEZE.HELPER_SOURCE_SHA256) == expected
    for relative, digest in FREEZE.HELPER_SOURCE_SHA256.items():
        assert FREEZE.sha256_file(FREEZE.ROOT / relative) == digest


def test_a14_biological_inputs_cover_held_lodo_dependencies() -> None:
    paths = FREEZE.a14_biological_inputs("hct116")
    assert "donors/k562_gwps/delta.npy" in paths
    assert "donors/orion_hek293t/delta.npy" in paths
    assert "donors/orion_hct116/delta.npy" not in paths
    assert "perturbmatch/delta.npy" in paths
    for condition in ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr"):
        assert f"cd4/{condition}/delta.npy" in paths
        assert f"cd4/{condition}/b_ctrl.npy" in paths
    assert "bundle/gene_names.csv" in paths
    assert "bundle/pert_counts.csv" in paths


def test_held_operator_requires_essential_prune() -> None:
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        path = root / "operator.npz"
        np.savez(
            path,
            evals=np.ones(2, dtype=np.float64),
            evecs=np.eye(2, dtype=np.float64),
        )
        metadata = {
            "axis": "union",
            "donors": ["k562_gwps", "orion_hek293t"],
            "n_cols": 2,
            "gamma": 0.65,
            "basis_prune": "essential",
            "dropped_in_panel": [],
            "sha256": FREEZE.sha256_file(path),
        }
        path.with_suffix(".json").write_text(json.dumps(metadata))
        observed = FREEZE.validate_held_operator("hct116", path)
        assert observed["basis_prune"] == "essential"
        metadata["basis_prune"] = "none"
        path.with_suffix(".json").write_text(json.dumps(metadata))
        try:
            FREEZE.validate_held_operator("hct116", path)
        except ValueError:
            pass
        else:
            raise AssertionError("unpruned held operator was accepted")


def test_chunked_control_sum_supports_backed_sparse_proxy() -> None:
    values = np.asarray(
        [[1, 0, 2], [0, 3, 0], [4, 0, 5], [0, 6, 0]],
        dtype=np.float32,
    )

    class BackedSparse:
        def __getitem__(self, key: object) -> sp.csr_matrix:
            return sp.csr_matrix(values[key])

    observed = FREEZE.sum_matrix_rows(BackedSparse(), 4, 3, chunk_rows=2)
    np.testing.assert_array_equal(observed, values.sum(axis=0))


def test_a14_reference_contract_rejects_biological_input_drift() -> None:
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        helper = root / "helper.py"
        biological = root / "gene_names.csv"
        helper.write_text("VALUE = 1\n")
        biological.write_text("gene_name\nA\n")
        reference = root / "reference.json"
        reference.write_text(
            json.dumps(
                {
                    "hashes": {
                        "inputs": {
                            "bundle/gene_names.csv": FREEZE.sha256_file(biological)
                        }
                    }
                }
            )
        )
        saved = (
            FREEZE.ROOT,
            FREEZE.A14_REFERENCE_META,
            FREEZE.A14_REFERENCE_META_SHA256,
            FREEZE.HELPER_SOURCE_SHA256,
            FREEZE.a14_biological_inputs,
            FREEZE.HELD_HARNESSES,
            FREEZE.HELD_PRUNED_OPERATORS,
        )
        FREEZE.ROOT = root
        FREEZE.A14_REFERENCE_META = reference
        FREEZE.A14_REFERENCE_META_SHA256 = FREEZE.sha256_file(reference)
        FREEZE.HELPER_SOURCE_SHA256 = {"helper.py": FREEZE.sha256_file(helper)}
        FREEZE.a14_biological_inputs = lambda surface: {
            "bundle/gene_names.csv": biological
        }
        FREEZE.HELD_HARNESSES = {"hct116": biological}
        FREEZE.HELD_PRUNED_OPERATORS = {"hct116": biological}
        biological.with_suffix(".h5ad.meta.json").write_text("metadata")
        biological.with_suffix(".json").write_text("metadata")
        (root / "scripts").mkdir()
        (root / "scripts/build_union_whitener.py").write_text("builder")
        (root / "data/donors/k562_essential").mkdir(parents=True)
        (root / "data/donors/k562_essential/perts.txt").write_text("A\n")
        try:
            FREEZE.validate_a14_reference("hct116")
            biological.write_text("gene_name\nB\n")
            try:
                FREEZE.validate_a14_reference("hct116")
            except ValueError:
                pass
            else:
                raise AssertionError("drifted biological input was accepted")
        finally:
            (
                FREEZE.ROOT,
                FREEZE.A14_REFERENCE_META,
                FREEZE.A14_REFERENCE_META_SHA256,
                FREEZE.HELPER_SOURCE_SHA256,
                FREEZE.a14_biological_inputs,
                FREEZE.HELD_HARNESSES,
                FREEZE.HELD_PRUNED_OPERATORS,
            ) = saved


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
