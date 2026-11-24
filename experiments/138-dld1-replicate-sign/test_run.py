from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np


MODULE_PATH = Path(__file__).with_name("run.py")
SPEC = importlib.util.spec_from_file_location("experiment138_run", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


class Experiment138Tests(unittest.TestCase):
    def proposal(
        self,
        low1: np.ndarray,
        low2: np.ndarray,
        measured1: np.ndarray | None = None,
        measured2: np.ndarray | None = None,
        supported: np.ndarray | None = None,
    ) -> RUN.Proposal:
        shape = low1.shape
        return RUN.Proposal(
            low1=np.asarray(low1, dtype=np.float64),
            low2=np.asarray(low2, dtype=np.float64),
            measured1=np.ones(shape, dtype=bool)
            if measured1 is None
            else np.asarray(measured1, dtype=bool),
            measured2=np.ones(shape, dtype=bool)
            if measured2 is None
            else np.asarray(measured2, dtype=bool),
            supported=np.ones(shape[0], dtype=bool)
            if supported is None
            else np.asarray(supported, dtype=bool),
        )

    def test_contiguous_folds_preserve_order_and_sizes(self) -> None:
        labels = RUN.contiguous_folds(613)
        self.assertEqual(
            [int(np.sum(labels == fold)) for fold in range(5)],
            [123, 123, 123, 122, 122],
        )
        for fold in range(5):
            positions = np.flatnonzero(labels == fold)
            np.testing.assert_array_equal(
                positions, np.arange(positions[0], positions[-1] + 1)
            )

    def test_actionable_requires_both_top_k_and_same_nonzero_sign(self) -> None:
        low1 = np.asarray([[9.0, 8.0, 7.0, 6.0, 5.0]])
        low2 = np.asarray([[9.0, -8.0, 1.0, 6.0, 0.0]])
        action, sign = RUN.actionable_coordinates(self.proposal(low1, low2), 3)
        np.testing.assert_array_equal(
            action, np.asarray([[True, False, False, False, False]])
        )
        np.testing.assert_array_equal(sign, np.asarray([[1, 0, 0, 0, 0]]))

    def test_top_k_ties_use_response_axis_order(self) -> None:
        magnitude = np.asarray([[3.0, 3.0, 3.0, 2.0]])
        eligible = np.ones_like(magnitude, dtype=bool)
        result = RUN.stable_top_k_mask(magnitude, eligible, 2)
        np.testing.assert_array_equal(result, [[True, True, False, False]])

    def test_transform_preserves_sign_norm_and_exact_fallback(self) -> None:
        incumbent = np.asarray(
            [
                [2.0, -3.0, 4.0, -5.0],
                [1.0, -2.0, 3.0, -4.0],
            ]
        )
        proposal = self.proposal(
            np.asarray([[9.0, -8.0, -7.0, 6.0], [5.0, 4.0, 3.0, 2.0]]),
            np.asarray([[8.0, -7.0, -6.0, 5.0], [4.0, 3.0, 2.0, 1.0]]),
            supported=np.asarray([True, False]),
        )
        transformed, report = RUN.transform_incumbent(
            incumbent, proposal, RUN.Arm("test", 3, 0.05)
        )
        np.testing.assert_array_equal(transformed[1], incumbent[1])
        np.testing.assert_array_equal(np.signbit(transformed), np.signbit(incumbent))
        np.testing.assert_allclose(
            np.linalg.norm(transformed, axis=1),
            np.linalg.norm(incumbent, axis=1),
            rtol=1e-12,
            atol=1e-12,
        )
        self.assertTrue(report["fallback_exact"])
        self.assertTrue(report["norm_preserved"])

    def test_identity_is_bit_exact(self) -> None:
        incumbent = np.asarray([[1.0, -0.0, 2.0]])
        proposal = self.proposal(incumbent, incumbent)
        transformed, _ = RUN.transform_incumbent(
            incumbent, proposal, RUN.ARM_BY_NAME["identity"]
        )
        self.assertTrue(
            np.array_equal(transformed.view(np.uint64), incumbent.view(np.uint64))
        )

    def test_target_null_moves_rows_only_within_active_set(self) -> None:
        low1 = np.arange(12, dtype=np.float64).reshape(3, 4)
        proposal = self.proposal(
            low1,
            low1 + 1,
            supported=np.asarray([True, False, True]),
        )
        permuted = RUN.permute_target_labels(
            proposal,
            np.asarray([0, 2]),
            np.asarray([1, 0]),
        )
        np.testing.assert_array_equal(permuted.low1[0], proposal.low1[2])
        np.testing.assert_array_equal(permuted.low1[2], proposal.low1[0])
        np.testing.assert_array_equal(permuted.low1[1], proposal.low1[1])
        np.testing.assert_array_equal(permuted.supported, proposal.supported)

    def test_response_null_moves_values_and_masks_together(self) -> None:
        low1 = np.asarray([[1.0, 2.0, 3.0]])
        measured = np.asarray([[True, False, True]])
        proposal = self.proposal(
            low1,
            low1 + 10,
            measured1=measured,
            measured2=~measured,
        )
        permutation = np.asarray([2, 0, 1])
        permuted = RUN.permute_response_coordinates(proposal, permutation)
        np.testing.assert_array_equal(permuted.low1, low1[:, permutation])
        np.testing.assert_array_equal(permuted.measured1, measured[:, permutation])
        np.testing.assert_array_equal(permuted.measured2, (~measured)[:, permutation])

    def test_complete_gallery_pds_for_query_subset(self) -> None:
        truth = np.eye(4, dtype=np.float64)
        query = np.asarray([1, 3])
        np.testing.assert_allclose(RUN.pds_queries(truth, truth, query), 1.0)
        np.testing.assert_allclose(
            RUN.pds_queries(np.zeros_like(truth), truth, query), 0.5
        )

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
            "selected_final_arm": "k050_beta0p025",
            "development_passes": True,
            "untouched_test_access_authorized": True,
        }
        result = {"final_arm": "identity", "passes": False}
        with self.assertRaises(RUN.ContractError):
            RUN.validate_candidate_freeze_semantics(freeze, result)

    def test_identity_wins_exact_selection_tie(self) -> None:
        score = {
            arm.name: {
                dataset: {
                    "pds": np.zeros(10),
                    "signed_credit": np.zeros(10),
                }
                for dataset in RUN.PARENT.SCORED_DATASETS
            }
            for arm in RUN.ARMS
        }
        selected, _ = RUN.select_arm(score, np.arange(10))
        self.assertEqual(selected.name, "identity")

    def test_negative_dataset_makes_arm_ineligible(self) -> None:
        summary = {
            "datasets": {
                dataset: {"pds": 0.01, "signed_credit": 0.01}
                for dataset in RUN.PARENT.SCORED_DATASETS
            },
            "macro": {"pds": 0.01, "signed_credit": 0.01},
        }
        summary["datasets"]["rpe1"]["pds"] = -1e-6
        self.assertFalse(RUN.eligible_for_selection(summary))

    def test_concentration(self) -> None:
        self.assertAlmostEqual(RUN.positive_share([1.0, 1.0, 1.0]), 1.0 / 3.0)
        self.assertEqual(RUN.positive_share([0.0, -1.0]), 1.0)
        self.assertEqual(RUN.positive_share([2.0, 1.0, -3.0]), 2.0 / 3.0)

    def test_metadata_contract_without_effect_arrays(self) -> None:
        report = RUN.metadata_report()
        self.assertFalse(report["effect_arrays_opened"])
        self.assertEqual(report["development_targets"], 613)
        self.assertEqual(report["development_dld1_supported"], 140)
        self.assertEqual(report["development_exact_fallback"], 473)
        self.assertEqual(report["test_targets"], 410)
        self.assertEqual(report["test_dld1_supported"], 94)
        self.assertEqual(report["test_exact_fallback"], 316)
        self.assertEqual(report["response_genes"], 4992)


if __name__ == "__main__":
    unittest.main()
