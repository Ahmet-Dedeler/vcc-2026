from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np


PATH = Path(__file__).with_name("run.py")
SPEC = importlib.util.spec_from_file_location("experiment140_run", PATH)
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


class Experiment140Tests(unittest.TestCase):
    def test_orion_negr1_is_legacy_only_and_never_candidate_control(self) -> None:
        self.assertTrue(RUN.is_orion_control("Non-Targeting"))
        self.assertTrue(RUN.is_legacy_orion_control("Non-Targeting"))
        self.assertTrue(RUN.is_legacy_orion_control("NEGR1"))
        self.assertFalse(RUN.is_orion_control("NEGR1"))
        self.assertFalse(RUN.is_orion_control("non-targeting"))
        self.assertFalse(RUN.is_orion_control("Non-Targeting "))
        self.assertEqual(
            RUN.ORION_CONTROL_SCAN_LABELS,
            ("Non-Targeting", "NEGR1"),
        )
        source = PATH.read_text()
        self.assertIn(
            '("gene_target", "=", label)',
            source,
        )
        self.assertIn("release_arrow_memory()", source)

    def test_union_pruned_selection_matches_live_a15(self) -> None:
        report = RUN.basis_selection_report(RUN.PRODUCTION_DONORS)
        self.assertEqual(
            report["rows_kept"],
            {
                "k562_gwps": 7784,
                "orion_hct116": 10807,
                "orion_hek293t": 10803,
            },
        )
        self.assertEqual(report["rows_kept_total"], 29394)
        self.assertEqual(
            report["row_index_sha256"],
            {
                "k562_gwps": (
                    "adf7cf446e9b3633fc05d07b66461a168a33325a8282a57403b092061e42f55a"
                ),
                "orion_hct116": (
                    "8d8001a999dc8270f39051998990e89ac233163d3465943d896690e23d65da36"
                ),
                "orion_hek293t": (
                    "40b54949f1626fa5f46d2b96b44e7ec08a46bf238bb42e39549a340cd536c462"
                ),
            },
        )
        k562_held = RUN.basis_selection_report(RUN.OUTER_DONORS["k562"])
        self.assertEqual(k562_held["rows_kept_total"], 21610)
        self.assertEqual(
            k562_held["rows_kept"],
            {"orion_hct116": 10807, "orion_hek293t": 10803},
        )

    def test_fit_operator_uses_frozen_union_pruned_selection(self) -> None:
        source = PATH.read_text()
        self.assertIn("selection = basis_selection(donors)", source)
        self.assertIn("selected_rows = selection[name]", source)
        self.assertNotIn(
            "total_rows = sum(len(geometry[name].perts) for name in donors)",
            source,
        )
        self.assertNotIn("delta.mean(axis=0, dtype=np.float64)", source)
        self.assertIn("delta.mean(axis=0), dtype=np.float32", source)
        self.assertIn("centered.astype(np.float64)", source)

    def test_umi_weight_derivation_matches_sum_then_normalize(self) -> None:
        batch_sums = np.asarray(
            [
                [90.0, 10.0],
                [10.0, 30.0],
            ]
        )
        batch_umi = batch_sums.sum(axis=1)
        normalized = batch_sums / batch_umi[:, None]
        umi_weighted = (batch_umi / batch_umi.sum()) @ normalized
        pooled = batch_sums.sum(axis=0) / batch_sums.sum()
        count_weighted = np.mean(normalized, axis=0)
        np.testing.assert_allclose(umi_weighted, pooled, rtol=0.0, atol=1e-15)
        self.assertFalse(np.allclose(count_weighted, pooled, rtol=0.0, atol=1e-8))

    def test_eta_zero_is_bit_exact_existing_row(self) -> None:
        base = np.arange(24, dtype=np.float32).reshape(6, 4)
        correction = np.linspace(-1, 1, 24, dtype=np.float32).reshape(6, 4)
        supported = np.asarray([True, True, False, True, False, True])
        output = RUN.blend_rows(base, correction, supported, 0.0)
        np.testing.assert_array_equal(output, base)

    def test_unsupported_fallback_is_bit_exact_at_nonzero_eta(self) -> None:
        base = np.arange(24, dtype=np.float64).reshape(6, 4)
        correction = np.ones_like(base)
        supported = np.asarray([True, False, True, False, True, True])
        output = RUN.blend_rows(base, correction, supported, 1.0)
        np.testing.assert_array_equal(output[~supported], base[~supported])
        np.testing.assert_array_equal(output[supported], base[supported] + 1.0)

    def test_supported_centering_is_zero_mean_and_fallback_compatible(self) -> None:
        correction = np.asarray(
            [
                [1.0, 3.0],
                [5.0, 7.0],
                [99.0, 99.0],
            ]
        )
        support = np.asarray([True, True, False])
        centered = np.zeros_like(correction)
        centered[support] = correction[support] - correction[support].mean(axis=0)
        np.testing.assert_allclose(centered[support].sum(axis=0), 0.0)
        np.testing.assert_array_equal(centered[~support], 0.0)

    def test_fold_sizes_are_frozen(self) -> None:
        folds = RUN.contiguous_folds(613)
        self.assertEqual(
            [int(np.sum(folds == fold)) for fold in range(5)],
            [123, 123, 123, 122, 122],
        )

    def test_target_and_batch_nulls_are_derangements(self) -> None:
        target_rng = np.random.default_rng(RUN.TARGET_NULL_SEED)
        batch_rng = np.random.default_rng(RUN.BATCH_NULL_SEED)
        target = RUN.stable_derangements(target_rng, 3, 31)
        batch = np.stack(
            [RUN.stable_derangements(batch_rng, 3, 31) for _ in range(3)],
            axis=1,
        )
        base = np.arange(31)
        self.assertTrue(np.all(target != base))
        self.assertTrue(np.all(batch != base[None, None, :]))
        self.assertFalse(np.array_equal(batch[:, 0], batch[:, 1]))

    def test_complete_gallery_subset_has_no_held_decoy_leak_in_selection(self) -> None:
        truth = np.asarray(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 1.0]]
        )
        identity = truth.copy()
        candidate = truth.copy()
        candidate[0] = [1.0, 0.0, 2.0]
        inner = np.asarray([0, 1])
        full = np.arange(3)
        np.testing.assert_array_equal(
            RUN.pds_subset(candidate, truth, np.asarray([0]), inner),
            RUN.pds_subset(identity, truth, np.asarray([0]), inner),
        )
        self.assertTrue(
            RUN.pds_subset(candidate, truth, np.asarray([0]), full)[0]
            < RUN.pds_subset(identity, truth, np.asarray([0]), full)[0]
        )

    def test_identity_wins_exact_selection_tie(self) -> None:
        scores = {
            eta: {
                dataset: {
                    "pds": np.zeros(8),
                    "signed_credit": np.zeros(8),
                }
                for dataset in RUN.SCORED_DATASETS
            }
            for eta in RUN.ETAS
        }
        selected, _ = RUN.select_eta(scores, np.arange(8))
        self.assertEqual(selected, 0.0)

    def test_selection_eligibility_is_by_biological_family(self) -> None:
        scores = {
            eta: {
                dataset: {
                    "pds": np.zeros(8),
                    "signed_credit": np.zeros(8),
                }
                for dataset in RUN.SCORED_DATASETS
            }
            for eta in RUN.ETAS
        }
        for metric in ("pds", "signed_credit"):
            scores[0.25]["k562_gwps"][metric][:] = -0.02
            scores[0.25]["k562_essential"][metric][:] = 0.04
            scores[0.25]["rpe1"][metric][:] = 0.01
            scores[0.25]["nadig_jurkat"][metric][:] = 0.01
        selected, summaries = RUN.select_eta(scores, np.arange(8))
        self.assertEqual(selected, 0.25)
        self.assertGreaterEqual(
            summaries[0.25]["families"]["pds"]["k562"],
            0.0,
        )
        self.assertLess(
            summaries[0.25]["datasets"]["k562_gwps"]["pds"],
            0.0,
        )

    def test_component_modes_are_frozen_and_disjoint(self) -> None:
        corrections = {
            "non_k562": {
                "k562_gwps": np.ones((2, 2)),
                "orion_hct116": 2 * np.ones((2, 2)),
                "orion_hek293t": 3 * np.ones((2, 2)),
            },
            "k562": {
                "orion_hct116": 4 * np.ones((2, 2)),
                "orion_hek293t": 5 * np.ones((2, 2)),
            },
        }
        self.assertEqual(
            set(RUN.component_corrections(corrections, "orion_only")["non_k562"]),
            {"orion_hct116", "orion_hek293t"},
        )
        self.assertEqual(
            set(RUN.component_corrections(corrections, "k562_only")["non_k562"]),
            {"k562_gwps"},
        )
        self.assertEqual(
            RUN.component_corrections(corrections, "k562_only")["k562"],
            {},
        )

    def test_result_checks_enforce_both_nested_and_fixed_inference(self) -> None:
        family = {"k562": 0.004, "rpe1": 0.004, "jurkat": 0.004}
        summary = {
            "families": {
                "pds": family,
                "signed_credit": family,
            },
            "macro": {
                "pds": 0.004,
                "signed_credit": 0.004,
            },
        }
        nested = {"final_eta": 0.25, "oof_summary": summary}
        boot = {
            scope: {
                metric: {"lower_95": 0.001}
                for metric in ("pds", "signed_credit")
            }
            for scope in ("nested", "fixed")
        }
        nulls = {
            kind: {
                scope: {
                    metric: {"p95": 0.001}
                    for metric in ("pds", "signed_credit")
                }
                for scope in ("nested", "fixed")
            }
            for kind in ("target_label", "batch_mixture")
        }
        concentration = {
            scope: {
                metric: {"family": 0.4, "target_block": 0.4}
                for metric in ("pds", "signed_credit")
            }
            for scope in ("nested", "fixed")
        }
        checks = RUN.result_checks(
            nested,
            summary,
            boot,
            nulls,
            concentration,
            {"identity_reproduction_passes": True},
        )
        self.assertTrue(all(checks.values()))
        self.assertIn("fixed_bootstrap_pds_lower_positive", checks)
        self.assertIn("fixed_target_label_null_signed_credit", checks)
        self.assertIn("fixed_target_block_concentration_pds", checks)

    def test_zero_positive_gain_fails_concentration(self) -> None:
        self.assertEqual(RUN.positive_share([0.0, -1.0, 0.0]), 1.0)

    def test_no_untouched_test_execution_path(self) -> None:
        choices = {
            "plan",
            "freeze",
            "prepare",
            "identity",
            "authorize-development",
            "develop",
            "audit",
        }
        source = PATH.read_text()
        self.assertNotIn("evaluate-test", choices)
        self.assertNotIn("load_test_truth", source)
        self.assertNotIn("untouched-test-result", source)


if __name__ == "__main__":
    unittest.main()
