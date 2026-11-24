#!/usr/bin/env python3
"""Synthetic tests for Experiment 107 strict Nadig validation."""

from __future__ import annotations

import ast
import hashlib
import importlib.util
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


HERE = Path(__file__).resolve().parent


def load_module():
    spec = importlib.util.spec_from_file_location(
        "exp107_strict_nadig_validation",
        HERE / "strict_nadig_validation.py",
    )
    if spec is None or spec.loader is None:
        raise ImportError(HERE / "strict_nadig_validation.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


STRICT = load_module()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_role_bundle(
    root: Path,
    role: str,
    arrays: dict[str, np.ndarray],
    *,
    parent_roles: list[str],
    extra_transitive: tuple[Path, ...] = (),
) -> tuple[Path, Path]:
    transitive = root / f"{role.replace('.', '-')}-source.txt"
    transitive.write_text("source\n", encoding="utf-8")
    npz_path = root / f"{role.replace('.', '-')}.npz"
    np.savez(npz_path, **arrays)
    manifest_path = root / f"{role.replace('.', '-')}.json"
    manifest_path.write_text(
        json.dumps(
            {
                "schema_version": STRICT.ROLE_BUNDLE_SCHEMA_VERSION,
                "experiment": "107-jurkat-high-moi",
                "status": STRICT.ROLE_BUNDLE_STATUS,
                "role": role,
                "kind": role.split(".", 1)[0],
                "npz": STRICT.file_binding(npz_path, root=root),
                "arrays": {
                    name: STRICT.array_binding(value) for name, value in arrays.items()
                },
                "parent_roles": parent_roles,
                "transitive_inputs": [
                    STRICT.file_binding(path, root=root)
                    for path in (transitive, *extra_transitive)
                ],
                "parameters": {},
            }
        ),
        encoding="utf-8",
    )
    return npz_path, manifest_path


def tiny_bundle_contract() -> dict[str, object]:
    targets = ("T1", "T2")
    scores = ("G1", "G2")
    columns = np.asarray([0, 2], dtype=np.int64)
    basis = np.asarray([0, 1], dtype=np.int64)
    return {
        "gallery_construction": {
            "targets_per_gallery": 2,
            "galleries": [
                {
                    "gallery": 1,
                    "target_axis_sha256_newline_terminated": STRICT.sha256_strings(
                        targets
                    ),
                    "score_gene_count": 2,
                    "score_axis_sha256_newline_terminated": STRICT.sha256_strings(
                        scores
                    ),
                    "whitener_column_count": 2,
                    "whitener_columns_int64_sha256": STRICT.sha256_int64(columns),
                    "whitener_basis_rows_by_donor": {
                        "k562_gwps": {
                            "count": 2,
                            "int64_sha256": STRICT.sha256_int64(basis),
                        }
                    },
                    "transferable_target_count": 1,
                }
            ],
        },
        "incumbent": {"donor_order": ["k562_gwps"]},
        "held_context_gate": {"surfaces": {}},
    }


def synthetic_source(
    directory: Path,
    *,
    valid: bool = True,
    drop_key: str | None = None,
) -> Path:
    directory.mkdir()
    targets, responses, exportable = 4, 3, 2
    folds, lambdas = 4, 5
    stability_pairs = np.asarray(
        [(left, right) for left in range(folds) for right in range(left + 1, folds)],
        dtype=np.int64,
    )
    arrays = {
        "exported_effects": np.asarray(
            [[1, 0, 0], [0, 1, 0], [0, 0, 0], [1, 1, 0]], dtype=np.float32
        ),
        "coefficients_binary_units": np.zeros((targets, responses), dtype=np.float32),
        "projected_coefficients": np.zeros((targets, responses), dtype=np.float32),
        "nuisance_response_slopes": np.zeros((3, responses)),
        "nuisance_target_slopes": np.zeros((3, targets)),
        "response_basis": np.eye(responses),
        "nullspace_basis": np.zeros((0, targets)),
        "evidence_support": np.ones(targets, dtype=bool),
        "zero_variance": np.zeros(targets, dtype=bool),
        "nuisance_aliased": np.zeros(targets, dtype=bool),
        "fitted": np.asarray([True, True, False, True]),
        "projection_diagonal": np.asarray([1.0, 1.0, 0.0, 1.0]),
        "transferable": np.asarray([True, True, False, True]),
        "target_scales": np.ones(targets),
        "singular_values": np.ones(targets),
        "retained_singular_directions": np.ones(targets, dtype=bool),
        "target_symbols": np.asarray(["A", "B", "C", "D"]),
        "response_symbols": np.asarray(["X", "Y", "Z"]),
        "response_vcc_indices": np.asarray([1, 2, -1], dtype=np.int32),
        "response_exportable": np.asarray([True, True, False]),
        "selected_lambda": np.asarray(3.0),
        "lambda_grid": np.asarray([1.0, 3.0, 10.0, 30.0, 100.0]),
        "held_cell_counts": np.ones(folds, dtype=np.int64),
        "fold_mean_gains": np.ones((folds, lambdas)),
        "mean_held_cell_gains": np.ones(lambdas),
        "fold_selected_effects": np.zeros(
            (folds, targets, responses), dtype=np.float32
        ),
        "fold_zero_variance": np.zeros((folds, targets), dtype=bool),
        "fold_nuisance_aliased": np.zeros((folds, targets), dtype=bool),
        "fold_projection_diagonal": np.ones((folds, targets)),
        "fold_transferable": np.ones((folds, targets), dtype=bool),
        "fold_support_training_channels": np.asarray(
            [
                [2, 3, 4, 6, 7, 8, 10, 11, 12, 14, 15, 16],
                [1, 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 16],
                [1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14, 16],
                [1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15],
            ],
            dtype=np.int64,
        ),
        "fold_support_distinct_called_guides": np.full(
            (folds, targets), 3, dtype=np.int64
        ),
        "fold_support_guides_with_three_untreated_cells": np.full(
            (folds, targets), 2, dtype=np.int64
        ),
        "fold_support_untreated_carrier_cells": np.full(
            (folds, targets), 20, dtype=np.int64
        ),
        "fold_support_training_channel_count": np.full(
            (folds, targets), 4, dtype=np.int64
        ),
        "fold_support_supported": np.ones((folds, targets), dtype=bool),
        "final_support_training_channels": np.arange(1, 17, dtype=np.int64),
        "final_support_distinct_called_guides": np.full(targets, 3, dtype=np.int64),
        "final_support_guides_with_three_untreated_cells": np.full(
            targets, 2, dtype=np.int64
        ),
        "final_support_untreated_carrier_cells": np.full(targets, 20, dtype=np.int64),
        "final_support_training_channel_count": np.full(targets, 4, dtype=np.int64),
        "final_support_supported": np.ones(targets, dtype=bool),
        "fold_unmodeled_indptr": np.arange(folds + 1, dtype=np.int64),
        "fold_unmodeled_target_burden": np.zeros(folds),
        "fold_stability_pairs": stability_pairs,
        "fold_effect_cosine": np.ones((len(stability_pairs), targets)),
        "fold_effect_sign_agreement": np.ones((len(stability_pairs), targets)),
        "projection_nullspace_annihilation_error": np.asarray(0.0),
        "projection_idempotence_error": np.asarray(0.0),
        "projection_gauge_shift_invariance_error": np.asarray(0.0),
    }
    if drop_key is not None:
        arrays.pop(drop_key)
    model = directory / "source-model.npz"
    np.savez_compressed(model, **arrays)
    manifest = {
        "status": "frozen_before_held_outcome_access",
        "experiment": "107-jurkat-high-moi",
        "primary_condition": "untreated",
        "held_truth_opened": False,
        "nadig_expression_or_scores_opened": False,
        "challenge_outcomes_opened": False,
        "model_file": "source-model.npz",
        "model_sha256": sha256(model) if valid else "0" * 64,
        "targets": targets,
        "response_genes": responses,
        "exportable_response_genes": exportable,
        "selected_lambda": 3.0,
        "lambda_grid": [1.0, 3.0, 10.0, 30.0, 100.0],
        "held_channel_folds": [
            [1, 5, 9, 13],
            [2, 6, 10, 14],
            [3, 7, 11, 15],
            [4, 8, 12, 16],
        ],
        "final_support": {"transferable_targets": 3},
        "final_solver": {
            "all_components_converged": True,
            "maximum_relative_residual": 1e-10,
            "relative_tolerance": 1e-8,
        },
        "projection_postconditions": {
            "tolerance": 1e-8,
            "nullspace_annihilation_error": 0.0,
            "projection_idempotence_error": 0.0,
            "gauge_shift_invariance_error": 0.0,
            "validated_after_serialization": True,
        },
    }
    (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return directory


def fallback_evidence() -> dict[str, object]:
    incumbents = [np.arange(12, dtype=np.float64).reshape(3, 4) for _ in range(6)]
    candidates = [value.copy() for value in incumbents]
    transferable = [np.asarray([True, False, True]) for _ in range(6)]
    coordinate_supported = [np.asarray([True, False, True, False]) for _ in range(6)]
    for candidate in candidates:
        candidate[np.ix_([0, 2], [0, 2])] += 0.25
    return {
        "incumbents": incumbents,
        "candidates": candidates,
        "transferable": transferable,
        "coordinate_supported": coordinate_supported,
    }


def held_context_metrics() -> dict[str, dict[str, float]]:
    return {
        "h1": {"pds_delta": 0.0, "directional_credit_delta": 0.0},
        "orion_hek293t": {
            "pds_delta": 0.0,
            "directional_credit_delta": 0.0,
        },
    }


def passing_nested() -> object:
    baseline = [np.full(3, 0.5) for _ in range(6)]
    pds = {
        0.0: baseline,
        0.05: [np.full(3, 0.504) for _ in range(6)],
        0.1: [np.full(3, 0.49) for _ in range(6)],
        0.2: [np.full(3, 0.48) for _ in range(6)],
    }
    return STRICT.nested_blend_selection(pds)


def passing_bootstrap(nested: object) -> dict[str, object]:
    draws = np.full(STRICT.BOOTSTRAP_REPLICATES, 1e-12)
    return {
        "repetitions": STRICT.BOOTSTRAP_REPLICATES,
        "point_estimate": nested.macro_gain,
        "lower_95": float(np.quantile(draws, 0.025, method="linear")),
        "draws": draws,
    }


class ContractAndSourceTests(unittest.TestCase):
    def test_frozen_contract_is_valid_json_and_hash(self) -> None:
        contract = STRICT.validate_contract()
        self.assertEqual(contract["gallery_construction"]["gallery_count"], 6)

    def test_contract_tampering_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "strict_validation_contract.json"
            path.write_bytes(STRICT.CONTRACT_PATH.read_bytes() + b" ")
            with self.assertRaisesRegex(STRICT.StrictValidationError, "bytes changed"):
                STRICT.validate_contract(path)

    def test_loader_schema_exactly_matches_current_fit_source_writer(self) -> None:
        tree = ast.parse((HERE / "fit_source.py").read_text(encoding="utf-8"))
        calls = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "atomic_savez"
        ]
        self.assertEqual(len(calls), 1)
        writer_keys = {
            keyword.arg for keyword in calls[0].keywords if keyword.arg is not None
        }
        self.assertEqual(writer_keys, STRICT.SOURCE_MODEL_KEYS)
        self.assertEqual(len(writer_keys), 50)

    def test_source_loader_validates_manifest_before_npz(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = synthetic_source(Path(temporary) / "source", valid=False)
            with mock.patch.object(STRICT.np, "load") as loader:
                with self.assertRaisesRegex(
                    STRICT.StrictValidationError, "hash changed"
                ):
                    STRICT._load_frozen_source_model(
                        source,
                        expected_targets=4,
                        expected_responses=3,
                        expected_exportable=2,
                        verify_current_implementations=False,
                    )
            loader.assert_not_called()

    def test_production_shaped_current_source_loads(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = synthetic_source(Path(temporary) / "source")
            with np.load(source / "source-model.npz", allow_pickle=False) as archive:
                self.assertEqual(set(archive.files), STRICT.SOURCE_MODEL_KEYS)
                self.assertEqual(len(archive.files), 50)
            loaded = STRICT._load_frozen_source_model(
                source,
                expected_targets=4,
                expected_responses=3,
                expected_exportable=2,
                verify_current_implementations=False,
            )
            np.testing.assert_array_equal(
                loaded.transferable, [True, True, False, True]
            )
            np.testing.assert_array_equal(loaded.exported_effects[2], 0.0)

    def test_source_schema_drift_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = synthetic_source(
                Path(temporary) / "source",
                drop_key="fold_effect_sign_agreement",
            )
            with self.assertRaisesRegex(
                STRICT.StrictValidationError,
                "source model keys changed",
            ):
                STRICT._load_frozen_source_model(
                    source,
                    expected_targets=4,
                    expected_responses=3,
                    expected_exportable=2,
                    verify_current_implementations=False,
                )

    def test_outcome_path_is_never_touched_before_valid_source_manifest(self) -> None:
        contract = STRICT.validate_contract()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = synthetic_source(root / "source", valid=False)
            plan = root / "plan.json"
            plan.write_text("{}", encoding="utf-8")
            outcome = root / "outcome"
            outcome.mkdir()
            (outcome / "delta.npy").write_bytes(b"must not be opened")
            receipt = root / "receipt.json"
            real_is_file = Path.is_file

            def guarded_is_file(path: Path) -> bool:
                if Path(path).name == "delta.npy":
                    raise AssertionError("outcome path touched before authorization")
                return real_is_file(path)

            with (
                mock.patch.object(Path, "is_file", guarded_is_file),
                mock.patch.object(STRICT, "SOURCE_DIRECTORY", source),
                mock.patch.object(STRICT, "NADIG_DIRECTORY", outcome),
                mock.patch.object(STRICT, "PRETRUTH_PLAN_PATH", plan),
                mock.patch.object(
                    STRICT,
                    "validate_contract",
                    return_value=contract,
                ),
                mock.patch.object(
                    STRICT,
                    "OUTCOME_ACCESS_RECEIPT_PATH",
                    receipt,
                ),
            ):
                with self.assertRaisesRegex(STRICT.StrictValidationError, "source"):
                    STRICT.load_nadig_outcome_after_authorization(
                        source,
                        plan,
                        outcome,
                        receipt_path=receipt,
                    )


class AuthorizationBindingTests(unittest.TestCase):
    def test_role_owned_bundles_allow_shared_transitive_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            contract = tiny_bundle_contract()
            donor_arrays = {
                "genes": np.asarray(["G1", "G2"]),
                "perts": np.asarray(["T1"]),
                "delta_shape": np.asarray([1, 2], dtype=np.int64),
                "b_ctrl_shape": np.asarray([2], dtype=np.int64),
            }
            donor_paths = write_role_bundle(
                root,
                "donor.k562_gwps",
                donor_arrays,
                parent_roles=[],
            )
            operator_arrays = {
                "evals": np.asarray([1.0, 2.0], dtype=np.float64),
                "evecs": np.eye(2, dtype=np.float32),
                "columns": np.asarray([0, 2], dtype=np.int64),
                "basis_k562_gwps": np.asarray([0, 1], dtype=np.int64),
            }
            operator_paths = write_role_bundle(
                root,
                "operator.nadig.gallery_1",
                operator_arrays,
                parent_roles=["donor.k562_gwps"],
                extra_transitive=donor_paths,
            )
            incumbent = np.arange(6, dtype=np.float64).reshape(2, 3)
            incumbent_arrays = {
                "targets": np.asarray(["T1", "T2"]),
                "genes": np.asarray(["G1", "GX", "G2"]),
                "score_genes": np.asarray(["G1", "G2"]),
                "score_indices": np.asarray([0, 2], dtype=np.int64),
                "effects": incumbent,
            }
            incumbent_paths = write_role_bundle(
                root,
                "incumbent.nadig.gallery_1",
                incumbent_arrays,
                parent_roles=[
                    "donor.k562_gwps",
                    "operator.nadig.gallery_1",
                ],
                extra_transitive=(*donor_paths, *operator_paths),
            )
            integration_arrays: dict[str, np.ndarray] = {
                "targets": incumbent_arrays["targets"],
                "genes": incumbent_arrays["genes"],
                "score_genes": incumbent_arrays["score_genes"],
                "score_indices": incumbent_arrays["score_indices"],
                "incumbent_effects": incumbent.copy(),
                "source_effects": np.asarray(
                    [[1.0, 0.0, 1.0], [0.0, 0.0, 0.0]], dtype=np.float64
                ),
                "transferable": np.asarray([True, False]),
                "coordinate_supported": np.asarray([True, False, True]),
            }
            for fraction in STRICT.FRACTIONS:
                suffix = f"f{int(round(fraction * 1000)):03d}"
                candidate = incumbent.copy()
                effective = np.zeros(2, dtype=bool)
                if fraction:
                    candidate[0, [0, 2]] += fraction
                    effective[0] = True
                integration_arrays[f"candidate_{suffix}"] = candidate
                integration_arrays[f"effective_{suffix}"] = effective
            integration_paths = write_role_bundle(
                root,
                "integration.nadig.gallery_1",
                integration_arrays,
                parent_roles=[
                    "incumbent.nadig.gallery_1",
                    "operator.nadig.gallery_1",
                ],
                extra_transitive=(*incumbent_paths, *operator_paths),
            )
            manifest_path = root / "inputs.json"
            artifacts = {
                "donor.k562_gwps": [path.name for path in donor_paths],
                "operator.nadig.gallery_1": [path.name for path in operator_paths],
                "incumbent.nadig.gallery_1": [path.name for path in incumbent_paths],
                "integration.nadig.gallery_1": [
                    path.name for path in integration_paths
                ],
            }
            manifest_path.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "experiment": "107-jurkat-high-moi",
                        "status": STRICT.ROLE_BUNDLE_STATUS,
                        "artifacts": artifacts,
                    }
                ),
                encoding="utf-8",
            )
            with mock.patch.object(STRICT, "OFFICIAL_GENE_COUNT", 3):
                bound = STRICT.bind_artifact_manifest(
                    manifest_path,
                    required_roles=tuple(artifacts),
                    root=root,
                    contract=contract,
                )
            self.assertEqual(set(bound["roles"]), set(artifacts))

    def test_artifact_reorder_and_hash_drift_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            transitive = root / "source.txt"
            transitive.write_text("source\n", encoding="utf-8")
            array_path = root / "donor.npz"
            arrays = {
                "genes": np.asarray(["A", "B"]),
                "perts": np.asarray(["P"]),
                "delta_shape": np.asarray([1, 2], dtype=np.int64),
                "b_ctrl_shape": np.asarray([2], dtype=np.int64),
            }
            np.savez(array_path, **arrays)
            bundle_path = root / "donor.json"
            bundle_path.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "experiment": "107-jurkat-high-moi",
                        "status": "frozen_before_experiment_107_nadig_outcome_access",
                        "role": "donor.k562_gwps",
                        "kind": "donor",
                        "npz": STRICT.file_binding(array_path, root=root),
                        "arrays": {
                            name: STRICT.array_binding(value)
                            for name, value in arrays.items()
                        },
                        "parent_roles": [],
                        "transitive_inputs": [
                            STRICT.file_binding(transitive, root=root)
                        ],
                        "parameters": {},
                    }
                ),
                encoding="utf-8",
            )
            manifest_path = root / "inputs.json"
            manifest_path.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "experiment": "107-jurkat-high-moi",
                        "status": "frozen_before_experiment_107_nadig_outcome_access",
                        "artifacts": {"donor.k562_gwps": ["donor.npz", "donor.json"]},
                    }
                ),
                encoding="utf-8",
            )
            required = ("donor.k562_gwps",)
            frozen = STRICT.bind_artifact_manifest(
                manifest_path,
                required_roles=required,
                root=root,
                contract={},
            )
            np.savez(
                array_path,
                genes=np.asarray(["B", "A"]),
                perts=arrays["perts"],
                delta_shape=arrays["delta_shape"],
                b_ctrl_shape=arrays["b_ctrl_shape"],
            )
            with self.assertRaisesRegex(
                STRICT.StrictValidationError,
                "changed",
            ):
                STRICT.verify_artifact_bindings(
                    frozen,
                    manifest_path,
                    required_roles=required,
                    root=root,
                    contract={},
                )

    def test_role_bundle_rejects_array_dtype_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            transitive = root / "source.txt"
            transitive.write_text("source\n", encoding="utf-8")
            array_path = root / "donor.npz"
            arrays = {
                "genes": np.asarray(["A", "B"]),
                "perts": np.asarray(["P"]),
                "delta_shape": np.asarray([1, 2], dtype=np.int32),
                "b_ctrl_shape": np.asarray([2], dtype=np.int64),
            }
            np.savez(array_path, **arrays)
            bundle_path = root / "donor.json"
            bundle_path.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "experiment": "107-jurkat-high-moi",
                        "status": STRICT.ROLE_BUNDLE_STATUS,
                        "role": "donor.k562_gwps",
                        "kind": "donor",
                        "npz": STRICT.file_binding(array_path, root=root),
                        "arrays": {
                            name: STRICT.array_binding(value)
                            for name, value in arrays.items()
                        },
                        "parent_roles": [],
                        "transitive_inputs": [
                            STRICT.file_binding(transitive, root=root)
                        ],
                        "parameters": {},
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(STRICT.StrictValidationError, "dtype changed"):
                STRICT.validate_role_bundle(
                    "donor.k562_gwps",
                    ["donor.npz", "donor.json"],
                    {},
                    root=root,
                )

    def test_operator_bundle_schema_is_strict(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            arrays = {
                "evals": np.asarray([1.0, 2.0], dtype=np.float64),
                "evecs": np.eye(2, dtype=np.float32),
                "columns": np.asarray([0, 2], dtype=np.int64),
                "basis_k562_gwps": np.asarray([0, 1], dtype=np.int64),
            }
            npz_path, manifest_path = write_role_bundle(
                root,
                "operator.nadig.gallery_1",
                arrays,
                parent_roles=["donor.k562_gwps"],
            )
            observed = STRICT.validate_role_bundle(
                "operator.nadig.gallery_1",
                [npz_path.name, manifest_path.name],
                tiny_bundle_contract(),
                root=root,
            )
            self.assertEqual(observed["kind"], "operator")
            arrays["evecs"] = np.eye(3, dtype=np.float32)
            np.savez(npz_path, **arrays)
            manifest = json.loads(manifest_path.read_text())
            manifest["npz"] = STRICT.file_binding(npz_path, root=root)
            manifest["arrays"]["evecs"] = STRICT.array_binding(arrays["evecs"])
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(STRICT.StrictValidationError, "shape changed"):
                STRICT.validate_role_bundle(
                    "operator.nadig.gallery_1",
                    [npz_path.name, manifest_path.name],
                    tiny_bundle_contract(),
                    root=root,
                )

    def test_incumbent_bundle_rejects_fortran_order(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            arrays = {
                "targets": np.asarray(["T1", "T2"]),
                "genes": np.asarray(["G1", "GX", "G2"]),
                "score_genes": np.asarray(["G1", "G2"]),
                "score_indices": np.asarray([0, 2], dtype=np.int64),
                "effects": np.asfortranarray(
                    np.arange(6, dtype=np.float64).reshape(2, 3)
                ),
            }
            npz_path, manifest_path = write_role_bundle(
                root,
                "incumbent.nadig.gallery_1",
                arrays,
                parent_roles=[
                    "donor.k562_gwps",
                    "operator.nadig.gallery_1",
                ],
            )
            with (
                mock.patch.object(STRICT, "OFFICIAL_GENE_COUNT", 3),
                self.assertRaisesRegex(STRICT.StrictValidationError, "not C-order"),
            ):
                STRICT.validate_role_bundle(
                    "incumbent.nadig.gallery_1",
                    [npz_path.name, manifest_path.name],
                    tiny_bundle_contract(),
                    root=root,
                )

    def test_integration_bundle_rejects_unsupported_coordinate_change(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            incumbent = np.arange(6, dtype=np.float64).reshape(2, 3)
            candidate = incumbent.copy()
            candidate[0, 1] += 1.0
            arrays: dict[str, np.ndarray] = {
                "targets": np.asarray(["T1", "T2"]),
                "genes": np.asarray(["G1", "GX", "G2"]),
                "score_genes": np.asarray(["G1", "G2"]),
                "score_indices": np.asarray([0, 2], dtype=np.int64),
                "incumbent_effects": incumbent,
                "source_effects": np.zeros((2, 3), dtype=np.float64),
                "transferable": np.asarray([True, False]),
                "coordinate_supported": np.asarray([True, False, True]),
            }
            for fraction in STRICT.FRACTIONS:
                suffix = f"f{int(round(fraction * 1000)):03d}"
                arrays[f"candidate_{suffix}"] = (
                    incumbent.copy() if fraction == 0.0 else candidate.copy()
                )
                arrays[f"effective_{suffix}"] = (
                    np.zeros(2, dtype=bool)
                    if fraction == 0.0
                    else np.asarray([True, False])
                )
            npz_path, manifest_path = write_role_bundle(
                root,
                "integration.nadig.gallery_1",
                arrays,
                parent_roles=[
                    "incumbent.nadig.gallery_1",
                    "operator.nadig.gallery_1",
                ],
            )
            with (
                mock.patch.object(STRICT, "OFFICIAL_GENE_COUNT", 3),
                self.assertRaisesRegex(
                    STRICT.StrictValidationError,
                    "unsupported target or coordinate changed",
                ),
            ):
                STRICT.validate_role_bundle(
                    "integration.nadig.gallery_1",
                    [npz_path.name, manifest_path.name],
                    tiny_bundle_contract(),
                    root=root,
                )

    def test_operator_missing_scored_coordinate_requires_exact_zero_fallback(self) -> None:
        genes = np.asarray(["G1", "MISSING", "G2"])
        score_indices = np.asarray([0, 1, 2], dtype=np.int64)
        operator_columns = np.asarray([0, 2], dtype=np.int64)
        incumbent = np.zeros((2, 3), dtype=np.float64)
        source = np.zeros_like(incumbent)
        candidates = {
            fraction: incumbent.copy() for fraction in STRICT.FRACTIONS
        }
        parameters = {
            "operator_missing_score_genes": ["MISSING"],
            "operator_missing_score_policy": (
                "retain_in_score_axis_require_zero_incumbent_and_source"
            ),
        }
        self.assertEqual(
            STRICT.assert_operator_missing_score_fallback(
                genes=genes,
                score_indices=score_indices,
                operator_columns=operator_columns,
                incumbent=incumbent,
                source=source,
                candidates=candidates,
                parameters=parameters,
                label="synthetic",
            ),
            ("MISSING",),
        )
        for field in ("incumbent", "source", "candidate"):
            changed_incumbent = incumbent.copy()
            changed_source = source.copy()
            changed_candidates = {
                fraction: value.copy() for fraction, value in candidates.items()
            }
            if field == "incumbent":
                changed_incumbent[0, 1] = 1.0
            elif field == "source":
                changed_source[0, 1] = 1.0
            else:
                changed_candidates[0.1][0, 1] = 1.0
            with self.assertRaises(STRICT.StrictValidationError):
                STRICT.assert_operator_missing_score_fallback(
                    genes=genes,
                    score_indices=score_indices,
                    operator_columns=operator_columns,
                    incumbent=changed_incumbent,
                    source=changed_source,
                    candidates=changed_candidates,
                    parameters=parameters,
                    label=f"synthetic-{field}",
                )

    def test_nadig_axis_reorder_and_delta_hash_drift_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            nadig = root / "nadig"
            nadig.mkdir()
            (nadig / "perts.txt").write_text("A\nB\n", encoding="utf-8")
            (nadig / "genes.txt").write_text("X\nY\n", encoding="utf-8")
            np.save(nadig / "delta.npy", np.eye(2, dtype=np.float32))
            with (
                mock.patch.object(STRICT, "ROOT", root),
                mock.patch.object(STRICT, "NADIG_DIRECTORY", nadig),
            ):
                frozen = STRICT.bind_nadig_artifacts(nadig)
                (nadig / "perts.txt").write_text("B\nA\n", encoding="utf-8")
                with self.assertRaisesRegex(
                    STRICT.StrictValidationError,
                    "target, gene, or delta artifact changed",
                ):
                    STRICT.verify_nadig_artifacts(frozen, nadig)
                (nadig / "perts.txt").write_text("A\nB\n", encoding="utf-8")
                np.save(
                    nadig / "delta.npy",
                    np.asarray([[0.0, 1.0], [1.0, 0.0]], dtype=np.float32),
                )
                with self.assertRaisesRegex(
                    STRICT.StrictValidationError,
                    "target, gene, or delta artifact changed",
                ):
                    STRICT.verify_nadig_artifacts(frozen, nadig)

    def test_immutable_receipt_rejects_replay(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            receipt = Path(temporary) / "receipt.json"
            STRICT.write_immutable_json(receipt, {"status": "first"})
            self.assertEqual(json.loads(receipt.read_text())["status"], "first")
            with self.assertRaisesRegex(
                STRICT.StrictValidationError,
                "replay refused",
            ):
                STRICT.write_immutable_json(receipt, {"status": "second"})

    def test_receipt_exists_before_synthetic_outcome_array_is_opened(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source_directory = root / "source"
            source_directory.mkdir()
            outcome_directory = root / "nadig"
            outcome_directory.mkdir()
            (outcome_directory / "perts.txt").write_text("A\nB\n", encoding="utf-8")
            (outcome_directory / "genes.txt").write_text("X\nY\nZ\n", encoding="utf-8")
            outcome_path = outcome_directory / "delta.npy"
            np.save(outcome_path, np.arange(6, dtype=np.float32).reshape(2, 3))
            plan_path = root / "plan.json"
            receipt_path = root / "receipt.json"
            fake_source = STRICT.FrozenSourceModel(
                directory=source_directory,
                manifest={},
                manifest_sha256="a" * 64,
                model_sha256="b" * 64,
                exported_effects=np.zeros((1, 1)),
                target_symbols=np.asarray(["A"]),
                response_symbols=np.asarray(["X"]),
                response_vcc_indices=np.asarray([0]),
                response_exportable=np.asarray([True]),
                transferable=np.asarray([True]),
            )
            with (
                mock.patch.object(STRICT, "ROOT", root),
                mock.patch.object(STRICT, "SOURCE_DIRECTORY", source_directory),
                mock.patch.object(STRICT, "NADIG_DIRECTORY", outcome_directory),
                mock.patch.object(
                    STRICT,
                    "OUTCOME_ACCESS_RECEIPT_PATH",
                    receipt_path,
                ),
            ):
                nadig_binding = STRICT.bind_nadig_artifacts(outcome_directory)
                plan = {
                    "nadig_target_count": 2,
                    "nadig_gene_count": 3,
                    "nadig_artifacts": nadig_binding,
                    "validation_implementation": {
                        "path": "strict.py",
                        "bytes": 1,
                        "sha256": "c" * 64,
                    },
                    "artifact_bindings": {"binding_sha256": "d" * 64},
                }
                plan_path.write_text(json.dumps(plan), encoding="utf-8")
                real_load = np.load

                def guarded_load(path: Path, *args: object, **kwargs: object):
                    self.assertTrue(receipt_path.is_file())
                    return real_load(path, *args, **kwargs)

                with (
                    mock.patch.object(
                        STRICT,
                        "authorize_outcome_phase",
                        return_value=(fake_source, plan),
                    ),
                    mock.patch.object(STRICT.np, "load", side_effect=guarded_load),
                ):
                    _, _, observed = STRICT.load_nadig_outcome_after_authorization(
                        source_directory,
                        plan_path,
                        outcome_directory,
                        receipt_path=receipt_path,
                    )
                np.testing.assert_array_equal(observed, np.arange(6).reshape(2, 3))
                with (
                    mock.patch.object(
                        STRICT,
                        "authorize_outcome_phase",
                        return_value=(fake_source, plan),
                    ) as authorization,
                    mock.patch.object(STRICT.np, "load") as loader,
                ):
                    with self.assertRaisesRegex(
                        STRICT.StrictValidationError,
                        "replay refused",
                    ):
                        STRICT.load_nadig_outcome_after_authorization(
                            source_directory,
                            plan_path,
                            outcome_directory,
                            receipt_path=receipt_path,
                        )
                authorization.assert_not_called()
                loader.assert_not_called()


class AxisTests(unittest.TestCase):
    def test_donor_axis_can_preserve_duplicates_for_last_occurrence_mapping(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "genes.txt"
            path.write_text("A\nB\nA\n", encoding="utf-8")
            self.assertEqual(
                STRICT._read_axis(  # noqa: SLF001
                    path, "donor gene axis", require_unique=False
                ),
                ("A", "B", "A"),
            )
            with self.assertRaisesRegex(
                STRICT.StrictValidationError, "duplicate symbols"
            ):
                STRICT._read_axis(path, "held gene axis")  # noqa: SLF001

    def test_sha_gallery_partition_is_deterministic_and_excludes_panel(self) -> None:
        source = [f"T{i:03d}" for i in range(20)]
        nadig = list(reversed(source))
        panel = {"T003", "T007"}
        first = STRICT.construct_gallery_targets(
            source, nadig, panel, take_first=12, gallery_size=3
        )
        second = STRICT.construct_gallery_targets(
            source, nadig, panel, take_first=12, gallery_size=3
        )
        self.assertEqual(first, second)
        self.assertEqual(len(first), 4)
        self.assertTrue(all(not (set(gallery) & panel) for gallery in first))
        self.assertEqual(len(set().union(*map(set, first))), 12)

    def test_score_axis_is_lexicographic_and_masked(self) -> None:
        observed = STRICT.score_axis(
            ["Z", "A", "PANEL", "HELD", "B"],
            [True, True, True, True, False],
            ["A", "B", "HELD", "Z", "OTHER", "PANEL"],
            ["PANEL"],
            ["HELD"],
        )
        self.assertEqual(observed, ("A", "Z"))

    def test_basis_selection_uses_one_rng_stream_and_excludes_gallery(self) -> None:
        donor_perts = {
            "k562_gwps": ["a", "b", "c", "d"],
            "orion_hct116": [f"h{i}" for i in range(10)],
            "orion_hek293t": [f"k{i}" for i in range(10)],
        }
        first = STRICT.select_whitener_basis_rows(
            donor_perts, {"b", "h1"}, {"c", "k2"}, row_cap=5
        )
        second = STRICT.select_whitener_basis_rows(
            donor_perts, {"b", "h1"}, {"c", "k2"}, row_cap=5
        )
        for name in STRICT.DONOR_ORDER:
            np.testing.assert_array_equal(first[name], second[name])
            selected_names = {donor_perts[name][index] for index in first[name]}
            self.assertFalse(selected_names & {"b", "h1", "c", "k2"})


class StatisticTests(unittest.TestCase):
    def test_pds_identity_and_swapped_rows(self) -> None:
        truth = np.eye(4)
        np.testing.assert_array_equal(
            STRICT.pds_per_target(truth, truth),
            np.ones(4),
        )
        swapped = truth[[1, 0, 2, 3]]
        self.assertLess(STRICT.pds_per_target(swapped, truth).mean(), 1.0)

    def test_unsupported_and_zero_rows_are_exactly_half(self) -> None:
        truth = np.eye(4)
        prediction = truth.copy()
        prediction[1] = 0
        observed = STRICT.pds_per_target(
            prediction,
            truth,
            np.asarray([True, True, False, True]),
        )
        self.assertEqual(observed[1], 0.5)
        self.assertEqual(observed[2], 0.5)

    def test_nested_selection_uses_five_galleries_and_smallest_tie(self) -> None:
        baseline = [np.full(3, 0.5) for _ in range(6)]
        pds = {
            0.0: baseline,
            0.05: [np.full(3, 0.6) for _ in range(6)],
            0.1: [np.full(3, 0.6) for _ in range(6)],
            0.2: [np.full(3, 0.4) for _ in range(6)],
        }
        result = STRICT.nested_blend_selection(pds)
        self.assertEqual(result.fractions_by_held_gallery, (0.05,) * 6)
        self.assertAlmostEqual(result.macro_gain, 0.1)
        self.assertEqual(STRICT.final_fraction(pds), 0.05)

    def test_nested_selection_does_not_select_on_held_gallery(self) -> None:
        baseline = [np.full(2, 0.5) for _ in range(6)]
        low = [np.full(2, 0.6) for _ in range(6)]
        high = [np.full(2, 0.7) for _ in range(6)]
        low[0] = np.full(2, 1.0)
        high[0] = np.full(2, 0.0)
        pds = {0.0: baseline, 0.05: low, 0.1: high, 0.2: baseline}
        result = STRICT.nested_blend_selection(pds)
        self.assertEqual(result.fractions_by_held_gallery[0], 0.1)

    def test_permutations_are_seeded_deterministic_and_keep_unsupported(self) -> None:
        source = np.arange(18).reshape(6, 3)
        supported = np.asarray([True, False, True, True, False, True])
        first = STRICT.permute_transferable_rows(
            source, supported, permutation_zero_based=7, gallery_one_based=3
        )
        second = STRICT.permute_transferable_rows(
            source, supported, permutation_zero_based=7, gallery_one_based=3
        )
        np.testing.assert_array_equal(first, second)
        np.testing.assert_array_equal(first[~supported], 0)
        self.assertCountEqual(
            map(tuple, first[supported].tolist()),
            map(tuple, source[supported].tolist()),
        )
        different = STRICT.permute_transferable_rows(
            source, supported, permutation_zero_based=8, gallery_one_based=3
        )
        self.assertFalse(np.array_equal(first, different))
        all_permutations = STRICT.deterministic_label_permutations(
            source, supported, gallery_one_based=3
        )
        self.assertEqual(len(all_permutations), 99)

    def test_blend_preserves_all_fallback_bytes_and_norm_matches(self) -> None:
        incumbent = np.asarray(
            [[1.0, -0.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]],
            dtype=np.float64,
        )
        source = np.asarray([[2.0, 9.0, 0.0], [0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])
        candidate, effective = STRICT.blend_with_source(
            incumbent,
            source,
            0.1,
            np.asarray([True, True, False]),
            np.asarray([True, False, True]),
        )
        np.testing.assert_array_equal(effective, [True, False, False])
        STRICT.assert_exact_fallback(
            incumbent, candidate, effective, np.asarray([True, False, True])
        )
        self.assertEqual(
            incumbent[0, 1].view(np.uint64),
            candidate[0, 1].view(np.uint64),
        )
        zero, zero_effective = STRICT.blend_with_source(
            incumbent, source, 0.0, np.ones(3, dtype=bool)
        )
        self.assertTrue(np.array_equal(zero.view(np.uint64), incumbent.view(np.uint64)))
        self.assertFalse(zero_effective.any())

    def test_principal_whitener_uses_supported_block_only(self) -> None:
        values = np.asarray([1.0, 4.0, 9.0])
        vectors = np.eye(3)
        source = np.asarray([[2.0, 3.0]])
        observed = STRICT.apply_supported_principal_whitener(
            source, values, vectors, np.asarray([0, 2]), gamma=1.0
        )
        np.testing.assert_allclose(observed, [[2.0, 1.0]])

    def test_cluster_labels_are_gallery_prefixed_and_deterministic(self) -> None:
        values = np.asarray([[1.0, 0.0], [0.99, 0.01], [-1.0, 0.0], [-0.99, -0.01]])
        first = STRICT.reciprocal_knn_clusters(
            values, gallery_one_based=4, k=1, minimum_cosine=0.8
        )
        second = STRICT.reciprocal_knn_clusters(
            values, gallery_one_based=4, k=1, minimum_cosine=0.8
        )
        np.testing.assert_array_equal(first, second)
        self.assertTrue(all(label.startswith("g4:") for label in first))
        self.assertEqual(len(set(first)), 2)

    def test_bootstrap_is_deterministic_and_repeats_nested_selection(self) -> None:
        baseline = [np.full(4, 0.5) for _ in range(6)]
        pds = {
            0.0: baseline,
            0.05: [np.asarray([0.55, 0.56, 0.54, 0.55]) for _ in range(6)],
            0.1: [np.asarray([0.60, 0.40, 0.60, 0.40]) for _ in range(6)],
            0.2: [np.full(4, 0.45) for _ in range(6)],
        }
        clusters = [
            np.asarray(
                [f"g{gallery}:a", f"g{gallery}:a", f"g{gallery}:b", f"g{gallery}:c"]
            )
            for gallery in range(1, 7)
        ]
        first = STRICT.gallery_stratified_cluster_bootstrap(
            pds, clusters, repetitions=40
        )
        second = STRICT.gallery_stratified_cluster_bootstrap(
            pds, clusters, repetitions=40
        )
        np.testing.assert_array_equal(first["draws"], second["draws"])
        self.assertGreater(first["lower_95"], 0)

    def test_concentration_uses_cluster_sums_over_exact_1800_denominator(self) -> None:
        deltas = [np.full(300, 0.01) for _ in range(6)]
        labels = [
            np.asarray([f"g{gallery}:c{index // 10}" for index in range(300)])
            for gallery in range(1, 7)
        ]
        result = STRICT.gain_concentration(deltas, labels, 0.01)
        self.assertAlmostEqual(result["fraction"], 0.1 / 18.0)
        self.assertTrue(result["passes"])
        failed = STRICT.gain_concentration(deltas, labels, 0.0)
        self.assertFalse(failed["passes"])
        self.assertTrue(math.isinf(failed["fraction"]))

    def test_promotion_gates_keep_strict_and_nonstrict_boundaries_distinct(
        self,
    ) -> None:
        nested = passing_nested()
        bootstrap = passing_bootstrap(nested)
        report = STRICT.promotion_gate_report(
            [0.55] * 6,
            nested,
            bootstrap,
            np.full(99, nested.macro_gain - 1e-12),
            {"fraction": 0.5},
            held_context_metrics(),
            fallback_evidence=fallback_evidence(),
        )
        self.assertTrue(report["passes"])
        failed_bootstrap = {
            "repetitions": 20_000,
            "point_estimate": nested.macro_gain,
            "lower_95": 0.0,
            "draws": np.zeros(20_000),
        }
        failed = STRICT.promotion_gate_report(
            [0.55] * 6,
            nested,
            failed_bootstrap,
            np.full(99, nested.macro_gain),
            {"fraction": 0.5},
            held_context_metrics(),
            fallback_evidence=fallback_evidence(),
        )
        self.assertFalse(failed["passes"])
        self.assertFalse(failed["checks"]["bootstrap_lower_95_strictly_above_0"])
        self.assertFalse(
            failed["checks"]["observed_gain_strictly_exceeds_permutation_maximum"]
        )

    def test_empty_nested_arrays_are_rejected(self) -> None:
        empty = STRICT.NestedSelection(
            fractions_by_held_gallery=(),
            baseline_gallery_means=np.asarray([]),
            selected_gallery_means=np.asarray([]),
            held_gallery_gains=np.asarray([]),
            macro_baseline=0.0,
            macro_selected=0.0,
            macro_gain=0.0,
        )
        with self.assertRaisesRegex(
            STRICT.StrictValidationError,
            "six frozen fractions",
        ):
            STRICT.promotion_gate_report(
                [0.55] * 6,
                empty,
                passing_bootstrap(passing_nested()),
                np.zeros(99),
                {"fraction": 0.0},
                held_context_metrics(),
                fallback_evidence=fallback_evidence(),
            )

    def test_nested_shape_mismatch_is_rejected(self) -> None:
        nested = passing_nested()
        malformed = STRICT.NestedSelection(
            fractions_by_held_gallery=nested.fractions_by_held_gallery,
            baseline_gallery_means=nested.baseline_gallery_means,
            selected_gallery_means=nested.selected_gallery_means[:-1],
            held_gallery_gains=nested.held_gallery_gains,
            macro_baseline=nested.macro_baseline,
            macro_selected=nested.macro_selected,
            macro_gain=nested.macro_gain,
        )
        with self.assertRaisesRegex(
            STRICT.StrictValidationError,
            "six gallery values",
        ):
            STRICT.promotion_gate_report(
                [0.55] * 6,
                malformed,
                passing_bootstrap(nested),
                np.zeros(99),
                {"fraction": 0.0},
                held_context_metrics(),
                fallback_evidence=fallback_evidence(),
            )

    def test_missing_h1_or_hek293t_gate_is_rejected(self) -> None:
        nested = passing_nested()
        with self.assertRaisesRegex(
            STRICT.StrictValidationError,
            "H1 and HEK293T",
        ):
            STRICT.promotion_gate_report(
                [0.55] * 6,
                nested,
                passing_bootstrap(nested),
                np.full(99, nested.macro_gain - 1e-12),
                {"fraction": 0.5},
                {"h1": held_context_metrics()["h1"]},
                fallback_evidence=fallback_evidence(),
            )

    def test_hct116_metrics_are_not_a_strict_promotion_input(self) -> None:
        nested = passing_nested()
        metrics = held_context_metrics()
        metrics["orion_hct116"] = {
            "pds_delta": -999.0,
            "directional_credit_delta": -999.0,
        }
        with self.assertRaisesRegex(
            STRICT.StrictValidationError,
            "H1 and HEK293T",
        ):
            STRICT.promotion_gate_report(
                [0.55] * 6,
                nested,
                passing_bootstrap(nested),
                np.full(99, nested.macro_gain - 1e-12),
                {"fraction": 0.5},
                metrics,
                fallback_evidence=fallback_evidence(),
            )

    def test_h1_and_hek293t_failures_block_promotion(self) -> None:
        nested = passing_nested()
        metrics = held_context_metrics()
        metrics["h1"]["pds_delta"] = -1e-12
        metrics["orion_hek293t"]["directional_credit_delta"] = -1e-12
        report = STRICT.promotion_gate_report(
            [0.55] * 6,
            nested,
            passing_bootstrap(nested),
            np.full(99, nested.macro_gain - 1e-12),
            {"fraction": 0.5},
            metrics,
            fallback_evidence=fallback_evidence(),
        )
        self.assertFalse(report["passes"])
        self.assertFalse(report["checks"]["h1_pds_delta_at_least_0.0"])
        self.assertFalse(
            report["checks"]["orion_hek293t_directional_credit_delta_at_least_0.0"]
        )

    def test_fallback_tamper_is_rejected(self) -> None:
        evidence = fallback_evidence()
        evidence["candidates"][0][1, 1] += 1.0
        with self.assertRaisesRegex(
            STRICT.StrictValidationError,
            "unsupported target or coordinate changed",
        ):
            STRICT.verify_exact_fallback_evidence(evidence)

    def test_hashed_fallback_artifact_tamper_is_rejected(self) -> None:
        evidence = fallback_evidence()
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "fallback.npz"
            arrays = {}
            for gallery in range(1, 7):
                index = gallery - 1
                arrays[f"gallery_{gallery}_incumbent"] = evidence["incumbents"][index]
                arrays[f"gallery_{gallery}_candidate"] = evidence["candidates"][index]
                arrays[f"gallery_{gallery}_transferable"] = evidence["transferable"][
                    index
                ]
                arrays[f"gallery_{gallery}_coordinate_supported"] = evidence[
                    "coordinate_supported"
                ][index]
            np.savez_compressed(path, **arrays)
            expected = sha256(path)
            verified = STRICT.verify_exact_fallback_evidence(
                {"artifact_path": path, "artifact_sha256": expected}
            )
            self.assertEqual(verified["artifact_sha256"], expected)
            with path.open("ab") as handle:
                handle.write(b"tamper")
            with self.assertRaisesRegex(
                STRICT.StrictValidationError,
                "artifact hash changed",
            ):
                STRICT.verify_exact_fallback_evidence(
                    {"artifact_path": path, "artifact_sha256": expected}
                )


if __name__ == "__main__":
    unittest.main()
