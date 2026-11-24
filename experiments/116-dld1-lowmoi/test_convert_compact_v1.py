from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import convert_compact_v1 as conversion  # noqa: E402


class ProtocolTests(unittest.TestCase):
    def test_protocol_is_source_only_and_forbids_outcome_actions(self) -> None:
        protocol = conversion.load_protocol()
        self.assertTrue(all(protocol["forbidden"].values()))
        self.assertTrue(
            protocol["author_value_preservation"][
                "store_values_without_centering_scaling_clipping_or_imputation"
            ]
        )
        self.assertFalse(
            protocol["duplicate_promoter_policy"]["selection_by_observed_effect"]
        )

    def test_real_metadata_mapping_matches_frozen_counts(self) -> None:
        protocol = conversion.load_protocol()
        genes = conversion.read_axis(protocol["inputs"]["canonical_gene_axis"])
        panel = conversion.read_axis(protocol["inputs"]["vcc_panel_axis"])
        labels = conversion.read_target_labels(protocol)
        _, groups, counts = conversion.map_promoter_labels(labels, genes)
        conversion.validate_metadata_mapping(
            protocol,
            labels,
            genes,
            groups,
            counts,
            panel,
        )
        self.assertEqual(len(groups), 4715)


class MappingTests(unittest.TestCase):
    def test_exact_target_mapping_and_no_alias_repair(self) -> None:
        labels = [
            "NTC_1",
            "ABC_P1",
            "ABC_ENST1(W)",
            "ABCD_P1",
            "abc_P1",
            "BARE",
        ]
        mapped, groups, counts = conversion.map_promoter_labels(
            labels,
            ["ABC", "AB", "BARE"],
        )
        self.assertEqual(mapped, [None, "ABC", "ABC", None, None, "BARE"])
        self.assertEqual(groups, {"ABC": [1, 2], "BARE": [5]})
        self.assertEqual(counts, {"ntc": 1, "mapped": 3, "unsupported": 2})

    def test_value_parser_preserves_zero_and_missingness(self) -> None:
        zero, measured = conversion.parse_value("-0.0")
        self.assertEqual(zero, -0.0)
        self.assertTrue(measured)
        missing, measured = conversion.parse_value("NA")
        self.assertTrue(math.isnan(missing))
        self.assertFalse(measured)
        with self.assertRaises(conversion.ConversionError):
            conversion.parse_value("inf")
        with self.assertRaises(conversion.ConversionError):
            conversion.parse_value("not-a-number")

    def test_duplicate_promoters_use_coordinatewise_finite_median(self) -> None:
        values = np.asarray(
            [
                [1.0, 3.0, 10.0],
                [np.nan, 5.0, 20.0],
                [np.nan, np.nan, 30.0],
            ]
        )
        result = conversion.aggregate_block(
            values,
            np.asarray([1]),
            np.asarray([2]),
            [(0, np.asarray([0, 1]))],
            2,
        )
        np.testing.assert_allclose(
            result,
            np.asarray([[2.0, 10.0], [5.0, 20.0], [np.nan, 30.0]]),
            equal_nan=True,
        )

    def test_panel_coordinates_are_zero_and_unmeasured(self) -> None:
        indices = np.asarray([0, 1, 2])
        aggregated = np.asarray(
            [[1.0, np.nan], [2.0, 3.0], [4.0, 5.0]],
            dtype=np.float64,
        )
        values, measured = conversion.apply_transfer_exclusion(
            indices,
            aggregated,
            np.asarray([False, True, False]),
        )
        np.testing.assert_array_equal(values[1], np.asarray([0.0, 0.0]))
        np.testing.assert_array_equal(measured[1], np.asarray([False, False]))
        self.assertFalse(measured[0, 1])
        self.assertEqual(values[0, 1], 0.0)


class DiagnosticTests(unittest.TestCase):
    def test_vector_metrics(self) -> None:
        a = np.asarray([1.0, -2.0, 3.0])
        b = np.asarray([2.0, -4.0, 6.0])
        mask = np.asarray([True, True, True])
        result = conversion.vector_metrics(a, b, mask, mask, top_k=3)
        self.assertAlmostEqual(result["cosine"], 1.0)
        self.assertAlmostEqual(result["pearson"], 1.0)
        self.assertAlmostEqual(result["top100_sign_agreement"], 1.0)

    def test_pds_is_one_for_identity_best_gallery(self) -> None:
        similarity = np.asarray(
            [
                [1.0, 0.1, 0.2],
                [0.3, 1.0, 0.4],
                [0.2, 0.1, 1.0],
            ]
        )
        np.testing.assert_allclose(conversion.pds_scores(similarity), 1.0)

    def test_pairwise_masked_cosine_respects_masks(self) -> None:
        query = np.asarray([[1.0, 10.0], [0.0, 1.0]])
        gallery = np.asarray([[2.0, -10.0], [1.0, 0.0]])
        query_mask = np.asarray([[True, False], [True, True]])
        gallery_mask = np.asarray([[True, False], [True, True]])
        result = conversion.pairwise_masked_cosine(
            query,
            gallery,
            query_mask,
            gallery_mask,
        )
        self.assertAlmostEqual(result[0, 0], 1.0)
        self.assertAlmostEqual(result[1, 1], 0.0)


if __name__ == "__main__":
    unittest.main()
