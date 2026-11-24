from __future__ import annotations

import copy
import json
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

import run_strict_evaluation as RUNNER
import strict_nadig_validation as STRICT


def binding(root: Path, path: Path) -> dict[str, object]:
    return {
        "path": path.relative_to(root).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": STRICT.sha256_file(path),
    }


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


class AmendmentTests(unittest.TestCase):
    def make_fixture(self, root: Path) -> dict[str, Path]:
        experiment = root / "experiments/107-jurkat-high-moi"
        experiment.mkdir(parents=True)
        runner = experiment / "runner.py"
        validator = experiment / "validator.py"
        receipt = experiment / "receipt.json"
        h1 = experiment / "h1.h5ad"
        h1_ref = experiment / "h1-ref.npz"
        hek = experiment / "hek.h5ad"
        hek_ref = experiment / "hek-ref.npz"
        artifact_binding = "a" * 64
        for path, contents in (
            (runner, b"runner\n"),
            (validator, b"validator\n"),
            (h1, b"opaque h1\n"),
            (h1_ref, b"opaque h1 ref\n"),
            (hek, b"opaque hek\n"),
            (hek_ref, b"opaque hek ref\n"),
        ):
            path.write_bytes(contents)
        write_json(
            receipt,
            {
                "status": "btbd11_validation_and_provenance_repair_complete",
                "artifact_binding_sha256": artifact_binding,
            },
        )
        plan = experiment / "plan.json"
        plan_payload = {
            "status": "pretruth_complete",
            "nadig_outcome_opened": False,
            "artifact_bindings": {"binding_sha256": artifact_binding},
            "validation_implementation": binding(root, validator),
        }
        write_json(plan, plan_payload)
        amendment = experiment / "amendment.json"
        amendment_payload = {
            "schema_version": 1,
            "experiment": RUNNER.EXPERIMENT,
            "status": "frozen_before_experiment_107_outcome_execution",
            "bindings": {
                "pretruth_plan": binding(root, plan),
                "strict_validator": binding(root, validator),
                "strict_artifact_binding_sha256": artifact_binding,
                "runner": binding(root, runner),
                "final_provenance_receipt": binding(root, receipt),
                "held_inputs": {
                    "h1": {
                        "harness": binding(root, h1),
                        "reference_de": binding(root, h1_ref),
                    },
                    "orion_hek293t": {
                        "harness": binding(root, hek),
                        "reference_de": binding(root, hek_ref),
                    },
                },
            },
            "execution_contract": {
                "fractions": list(STRICT.FRACTIONS),
                "gallery_count": 6,
                "targets_per_gallery": 300,
                "permutation_replicates": 99,
                "cluster_bootstrap_replicates": 20_000,
                "held_surfaces": ["h1", "orion_hek293t"],
                "outcome_loader_must_follow_authentication": True,
                "immutable_result_directory": (
                    "experiments/107-jurkat-high-moi/work/"
                    "strict-evaluation-results"
                ),
            },
        }
        write_json(amendment, amendment_payload)
        return {
            "root": root,
            "runner": runner,
            "validator": validator,
            "plan": plan,
            "amendment": amendment,
            "artifact_binding": artifact_binding,
        }

    def authenticate(self, fixture: dict[str, Path]) -> dict[str, object]:
        return RUNNER.authenticate_execution_amendment(
            fixture["amendment"],
            root=fixture["root"],
            runner_path=fixture["runner"],
            plan_path=fixture["plan"],
            validator_path=fixture["validator"],
            expected_artifact_binding=str(fixture["artifact_binding"]),
        )

    def test_authenticates_complete_binding(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = self.make_fixture(Path(temporary))
            observed = self.authenticate(fixture)
            self.assertEqual(
                observed["artifact_binding_sha256"], fixture["artifact_binding"]
            )
            self.assertEqual(
                set(observed["held_paths"]), {"h1", "orion_hek293t"}
            )

    def test_every_primary_binding_mismatch_is_rejected(self) -> None:
        mutations = {
            "runner": lambda fixture: fixture["runner"].write_text("changed\n"),
            "validator": lambda fixture: fixture["validator"].write_text("changed\n"),
            "plan": lambda fixture: fixture["plan"].write_text("{}\n"),
            "artifact": self._mutate_artifact_binding,
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temporary:
                fixture = self.make_fixture(Path(temporary))
                mutate(fixture)
                touched = False

                def outcome_loader(*args: object, **kwargs: object) -> object:
                    nonlocal touched
                    touched = True
                    raise AssertionError("outcome loader must not be called")

                with self.assertRaises(RUNNER.EvaluationDriverError):
                    RUNNER.authorize_then_load_nadig(
                        fixture["amendment"],
                        outcome_loader=outcome_loader,
                        authentication_kwargs={
                            "root": fixture["root"],
                            "runner_path": fixture["runner"],
                            "plan_path": fixture["plan"],
                            "validator_path": fixture["validator"],
                            "expected_artifact_binding": str(
                                fixture["artifact_binding"]
                            ),
                        },
                    )
                self.assertFalse(touched)

    @staticmethod
    def _mutate_artifact_binding(fixture: dict[str, Path]) -> None:
        amendment = json.loads(fixture["amendment"].read_text())
        amendment["bindings"]["strict_artifact_binding_sha256"] = "b" * 64
        write_json(fixture["amendment"], amendment)

    def test_authentication_failure_precedes_outcome_loader(self) -> None:
        touched = False

        def outcome_loader(*args: object, **kwargs: object) -> object:
            nonlocal touched
            touched = True
            raise AssertionError("outcome loader must not be called")

        with mock.patch.object(
            RUNNER,
            "authenticate_execution_amendment",
            side_effect=RUNNER.EvaluationDriverError("bad amendment"),
        ):
            with self.assertRaisesRegex(
                RUNNER.EvaluationDriverError, "bad amendment"
            ):
                RUNNER.authorize_then_load_nadig(outcome_loader=outcome_loader)
        self.assertFalse(touched)


class FrozenEvaluationTests(unittest.TestCase):
    @staticmethod
    def galleries(count: int = 6) -> list[RUNNER.GalleryInput]:
        result = []
        for _ in range(count):
            base = np.zeros((300, 2), dtype=np.float64)
            source = np.ones((300, 2), dtype=np.float64)
            truth = np.ones((300, 2), dtype=np.float64)
            result.append(
                RUNNER.GalleryInput(
                    incumbent=base,
                    source=source,
                    truth=truth,
                    transferable=np.ones(300, dtype=bool),
                    cluster_basis=source,
                )
            )
        return result

    @staticmethod
    def fallback() -> dict[str, tuple[np.ndarray, ...]]:
        zeros = tuple(np.zeros((300, 2), dtype=np.float64) for _ in range(6))
        return {
            "incumbents": zeros,
            "candidates": tuple(item.copy() for item in zeros),
            "transferable": tuple(np.zeros(300, dtype=bool) for _ in range(6)),
            "coordinate_supported": tuple(
                np.ones(2, dtype=bool) for _ in range(6)
            ),
        }

    @staticmethod
    def held() -> dict[str, dict[str, float]]:
        return {
            "h1": {"pds_delta": 0.0, "directional_credit_delta": 0.0},
            "orion_hek293t": {
                "pds_delta": 0.0,
                "directional_credit_delta": 0.0,
            },
        }

    @staticmethod
    def fake_strict() -> SimpleNamespace:
        calls: dict[str, object] = {}
        grid = {
            0.0: tuple(np.full(300, 0.50) for _ in range(6)),
            0.05: tuple(np.full(300, 0.51) for _ in range(6)),
            0.1: tuple(np.full(300, 0.52) for _ in range(6)),
            0.2: tuple(np.full(300, 0.51) for _ in range(6)),
        }

        def pds_fraction_grid(*args: object, **kwargs: object) -> object:
            calls["grid_gallery_count"] = len(args[0])
            return copy.deepcopy(grid)

        def permutation(*args: object, **kwargs: object) -> np.ndarray:
            calls["permutations"] = kwargs["repetitions"]
            return np.zeros(kwargs["repetitions"], dtype=np.float64)

        def clusters(
            values: np.ndarray, *, gallery_one_based: int
        ) -> np.ndarray:
            return np.asarray(
                [f"g{gallery_one_based}:c{index}" for index in range(len(values))]
            )

        def bootstrap(*args: object, **kwargs: object) -> dict[str, object]:
            calls["bootstrap"] = kwargs["repetitions"]
            point = STRICT.nested_blend_selection(args[0]).macro_gain
            draws = np.full(kwargs["repetitions"], point, dtype=np.float64)
            return {
                "repetitions": kwargs["repetitions"],
                "point_estimate": point,
                "lower_95": point,
                "upper_95": point,
                "draws": draws,
            }

        def pds_per_target(
            source: np.ndarray, truth: np.ndarray, support: np.ndarray
        ) -> np.ndarray:
            return np.full(len(source), 0.60, dtype=np.float64)

        namespace = SimpleNamespace(
            FRACTIONS=STRICT.FRACTIONS,
            PERMUTATION_REPLICATES=STRICT.PERMUTATION_REPLICATES,
            BOOTSTRAP_REPLICATES=STRICT.BOOTSTRAP_REPLICATES,
            pds_fraction_grid=pds_fraction_grid,
            nested_blend_selection=STRICT.nested_blend_selection,
            final_fraction=STRICT.final_fraction,
            permutation_null_distribution=permutation,
            reciprocal_knn_clusters=clusters,
            gallery_stratified_cluster_bootstrap=bootstrap,
            gain_concentration=STRICT.gain_concentration,
            pds_per_target=pds_per_target,
            promotion_gate_report=STRICT.promotion_gate_report,
            calls=calls,
        )
        return namespace

    def test_full_frozen_counts_and_outputs_are_requested(self) -> None:
        fake = self.fake_strict()
        result = RUNNER.evaluate_nadig_statistics(
            self.galleries(),
            fallback_evidence=self.fallback(),
            held_context_metrics=self.held(),
            strict_module=fake,
        )
        self.assertEqual(fake.calls["grid_gallery_count"], 6)
        self.assertEqual(fake.calls["permutations"], 99)
        self.assertEqual(fake.calls["bootstrap"], 20_000)
        self.assertEqual(result["selected_fraction"], 0.1)
        self.assertEqual(
            result["nested"].fractions_by_held_gallery, (0.1,) * 6
        )
        self.assertEqual(result["permutation_gains"].shape, (99,))
        self.assertEqual(result["bootstrap"]["draws"].shape, (20_000,))
        self.assertEqual(result["pure_source_gallery_cosine"].shape, (6,))
        self.assertIn("fraction", result["concentration"])
        self.assertEqual(
            result["gate_report"]["fallback_evidence"]["gallery_count"], 6
        )
        self.assertEqual(
            set(result["gate_report"]["held_context_metrics"]),
            {"h1", "orion_hek293t"},
        )

    def test_evaluation_is_deterministic(self) -> None:
        first = RUNNER.evaluate_nadig_statistics(
            self.galleries(),
            fallback_evidence=self.fallback(),
            held_context_metrics=self.held(),
            strict_module=self.fake_strict(),
        )
        second = RUNNER.evaluate_nadig_statistics(
            self.galleries(),
            fallback_evidence=self.fallback(),
            held_context_metrics=self.held(),
            strict_module=self.fake_strict(),
        )
        self.assertEqual(
            RUNNER._nested_payload(first["nested"]),  # noqa: SLF001
            RUNNER._nested_payload(second["nested"]),  # noqa: SLF001
        )
        np.testing.assert_array_equal(
            first["permutation_gains"], second["permutation_gains"]
        )
        np.testing.assert_array_equal(
            first["bootstrap"]["draws"], second["bootstrap"]["draws"]
        )

    def test_exactly_six_galleries_are_required(self) -> None:
        with self.assertRaisesRegex(
            RUNNER.EvaluationDriverError, "exactly six"
        ):
            RUNNER.evaluate_nadig_statistics(
                self.galleries(5),
                fallback_evidence=self.fallback(),
                held_context_metrics=self.held(),
                strict_module=self.fake_strict(),
            )


class DirectionalCreditTests(unittest.TestCase):
    def test_directional_credit_rewards_matching_sign(self) -> None:
        reference = RUNNER.ReferenceDE(
            genes=np.asarray(["A", "B", "C"]),
            gate=np.asarray([True, True, True]),
            gate_indices=np.asarray([0, 1, 2]),
            control_mean_cpm=np.asarray([100.0, 100.0, 100.0]),
            targets=np.asarray(["A"]),
            lfc=np.asarray([[0.0, 1.0, -1.0]]),
            significant=(np.asarray([1, 2]),),
        )
        control = np.log1p(
            RUNNER.TS * np.asarray([0.2, 0.4, 0.4], dtype=np.float64)
        )
        matching = np.asarray([[0.0, 0.2, -0.2]])
        reversed_effect = -matching
        self.assertGreater(
            RUNNER.directional_credit_per_target(
                matching, control, ["A"], reference
            )[0],
            RUNNER.directional_credit_per_target(
                reversed_effect, control, ["A"], reference
            )[0],
        )


class HeldLoaderTests(unittest.TestCase):
    def test_backed_sparse_pseudobulk_loader_on_synthetic_h5ad(self) -> None:
        import anndata as ad
        import pandas as pd
        import scipy.sparse as sp

        work = RUNNER.HERE / "work"
        work.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=work) as temporary:
            path = Path(temporary) / "synthetic.h5ad"
            matrix = sp.csr_matrix(
                np.asarray(
                    [
                        [10, 0, 5],
                        [20, 0, 5],
                        [0, 20, 10],
                        [0, 30, 10],
                        [5, 5, 20],
                        [5, 5, 30],
                    ],
                    dtype=np.float32,
                )
            )
            obs = pd.DataFrame(
                {
                    "target_gene": [
                        "non-targeting",
                        "non-targeting",
                        "T1",
                        "T1",
                        "T2",
                        "T2",
                    ]
                },
                index=[f"cell-{index}" for index in range(6)],
            )
            var = pd.DataFrame(index=["A", "B", "C"])
            ad.AnnData(X=matrix, obs=obs, var=var).write_h5ad(path)
            genes, control, effects = RUNNER._load_held_effects(  # noqa: SLF001
                path, ("T1", "T2")
            )
            np.testing.assert_array_equal(genes, np.asarray(["A", "B", "C"]))
            self.assertEqual(control.shape, (3,))
            self.assertEqual(effects.shape, (2, 3))
            self.assertTrue(np.isfinite(effects).all())


class ImmutableOutputTests(unittest.TestCase):
    def test_outputs_are_read_only_and_replay_is_refused(self) -> None:
        work = RUNNER.HERE / "work"
        work.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=work) as temporary:
            destination = Path(temporary) / "result"
            RUNNER.write_immutable_results(
                destination,
                summary={"finite": 1.0},
                arrays={"values": np.arange(3)},
                input_bindings=[],
            )
            for name in (
                RUNNER.RESULT_ARRAYS_NAME,
                RUNNER.RESULT_SUMMARY_NAME,
                RUNNER.RESULT_MANIFEST_NAME,
            ):
                self.assertEqual(
                    os.stat(destination / name).st_mode & 0o222,
                    0,
                )
            with self.assertRaisesRegex(
                RUNNER.EvaluationDriverError, "replay refused"
            ):
                RUNNER.write_immutable_results(
                    destination,
                    summary={"finite": 1.0},
                    arrays={"values": np.arange(3)},
                    input_bindings=[],
                )


if __name__ == "__main__":
    unittest.main()
