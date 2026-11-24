from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np

PATH = Path(__file__).with_name("run.py")
SPEC = importlib.util.spec_from_file_location("experiment139_run", PATH)
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


class Experiment139Tests(unittest.TestCase):
    def test_metadata_contract_does_not_open_effect_values(self) -> None:
        result = RUN.metadata_report()
        self.assertFalse(result["effect_values_opened"])
        self.assertEqual(result["development_transferable"], 609)
        self.assertEqual(result["test_transferable"], 409)
        self.assertEqual(result["shared_exportable_landmarks"], 138)
        self.assertEqual(result["response_genes"], 4992)

    def test_contiguous_fold_sizes(self) -> None:
        folds = RUN.contiguous_folds(613)
        self.assertEqual(
            [int(np.sum(folds == x)) for x in range(5)],
            [123, 123, 123, 122, 122],
        )

    def test_routing_weights_are_simplex(self) -> None:
        similarity = np.asarray([[0.2, 0.1, -0.1, 0.0, 0.3]])
        for arm in RUN.ARMS:
            weights = RUN.routing_weights(similarity, arm)
            self.assertTrue(np.all(weights >= 0))
            np.testing.assert_allclose(weights.sum(axis=1), 1.0, atol=1e-14)

    def test_identity_is_exact_equal_mean(self) -> None:
        sources = np.arange(60, dtype=np.float64).reshape(5, 3, 4)
        evidence = RUN.LandmarkEvidence(np.ones((3, 2)), np.ones(3, dtype=bool))
        candidate, _ = RUN.routed_candidate(
            sources, np.asarray([0, 1]), evidence, RUN.ARM_BY_NAME["identity"]
        )
        np.testing.assert_array_equal(candidate, np.mean(sources, axis=0))

    def test_candidate_is_convex_and_fallback_exact(self) -> None:
        rng = np.random.default_rng(4)
        sources = rng.normal(size=(5, 3, 6))
        evidence = RUN.LandmarkEvidence(
            rng.normal(size=(3, 3)), np.asarray([True, False, True])
        )
        old = RUN.MIN_NONZERO_LANDMARKS
        RUN.MIN_NONZERO_LANDMARKS = 1
        try:
            candidate, report = RUN.routed_candidate(
                sources,
                np.asarray([0, 2, 4]),
                evidence,
                RUN.ARM_BY_NAME["rho050_temp010"],
            )
        finally:
            RUN.MIN_NONZERO_LANDMARKS = old
        identity = np.mean(sources, axis=0)
        np.testing.assert_array_equal(candidate[1], identity[1])
        self.assertTrue(np.all(candidate >= sources.min(axis=0) - 1e-12))
        self.assertTrue(np.all(candidate <= sources.max(axis=0) + 1e-12))
        self.assertTrue(report["simplex"])

    def test_positive_rescaling_of_landmark_does_not_change_candidate(self) -> None:
        rng = np.random.default_rng(8)
        sources = rng.normal(size=(5, 2, 5))
        evidence = RUN.LandmarkEvidence(rng.normal(size=(2, 3)), np.ones(2, dtype=bool))
        old = RUN.MIN_NONZERO_LANDMARKS
        RUN.MIN_NONZERO_LANDMARKS = 1
        try:
            first, _ = RUN.routed_candidate(
                sources, np.asarray([0, 2, 4]), evidence, RUN.ARMS[1]
            )
            second, _ = RUN.routed_candidate(
                sources,
                np.asarray([0, 2, 4]),
                RUN.LandmarkEvidence(evidence.values * 100.0, evidence.supported),
                RUN.ARMS[1],
            )
        finally:
            RUN.MIN_NONZERO_LANDMARKS = old
        np.testing.assert_allclose(first, second, rtol=1e-12, atol=1e-12)

    def test_target_null_preserves_unsupported_rows(self) -> None:
        evidence = RUN.LandmarkEvidence(
            np.arange(12, dtype=np.float64).reshape(3, 4),
            np.asarray([True, False, True]),
        )
        changed = RUN.target_null(evidence, np.asarray([1, 0]))
        np.testing.assert_array_equal(changed.values[0], evidence.values[2])
        np.testing.assert_array_equal(changed.values[1], evidence.values[1])
        np.testing.assert_array_equal(changed.supported, evidence.supported)

    def test_landmark_null_is_coordinate_only(self) -> None:
        evidence = RUN.LandmarkEvidence(
            np.arange(6, dtype=np.float64).reshape(2, 3), np.ones(2, dtype=bool)
        )
        changed = RUN.landmark_null(evidence, np.asarray([2, 0, 1]))
        np.testing.assert_array_equal(changed.values, evidence.values[:, [2, 0, 1]])
        np.testing.assert_array_equal(changed.supported, evidence.supported)

    def test_inner_gallery_excludes_held_decoys(self) -> None:
        candidate = np.asarray(
            [
                [0.8, 0.1, 0.9],
                [0.1, 0.8, 0.9],
                [0.0, 0.0, 1.0],
            ]
        )
        identity = np.asarray(
            [
                [0.8, 0.1, 0.0],
                [0.1, 0.8, 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        inner = np.asarray([0, 1])
        full = np.arange(3)
        np.testing.assert_array_equal(
            RUN.pds_from_similarity(candidate, inner, inner),
            RUN.pds_from_similarity(identity, inner, inner),
        )
        self.assertTrue(
            np.all(
                RUN.pds_from_similarity(candidate, inner, full)
                < RUN.pds_from_similarity(identity, inner, full)
            )
        )

    def test_candidate_freeze_semantics_cannot_override_result(self) -> None:
        freeze = {
            "final_arm": "rho025_temp025",
            "development_passes": True,
            "test_authorized": True,
        }
        result = {"final_arm": "identity", "passes": False}
        with self.assertRaises(RUN.ContractError):
            RUN.validate_candidate_semantics(freeze, result)

    def test_identity_wins_exact_tie(self) -> None:
        scores = {
            arm.name: {
                dataset: {
                    "pds": np.zeros(8),
                    "signed_credit": np.zeros(8),
                }
                for dataset in RUN.PARENT.SCORED_DATASETS
            }
            for arm in RUN.ARMS
        }
        selected, _ = RUN.select_arm(scores, np.arange(8))
        self.assertEqual(selected.name, "identity")

    def test_zero_positive_gain_fails_concentration(self) -> None:
        self.assertEqual(RUN.positive_share([0.0, -1.0, 0.0]), 1.0)


if __name__ == "__main__":
    unittest.main()
