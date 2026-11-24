from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.stats import binom


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import guide_calling  # noqa: E402


def dense_bh_canonical(p_values: np.ndarray) -> np.ndarray:
    """Independent dense BH reference with stable canonical-index ties."""

    values = np.asarray(p_values, dtype=np.float64)
    order = np.argsort(values, kind="stable")
    ranked = values[order] * values.size / np.arange(1, values.size + 1)
    adjusted_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(values)
    adjusted[order] = np.minimum(adjusted_sorted, 1.0)
    return adjusted


def make_batch(
    cells: list[tuple[int, list[int], list[int]]],
) -> guide_calling.SparseGuideBatch:
    cell_indices = np.asarray([cell for cell, _, _ in cells], dtype=np.int64)
    lengths = np.asarray([len(guides) for _, guides, _ in cells], dtype=np.int64)
    indptr = np.r_[0, np.cumsum(lengths)]
    guide_indices = np.asarray(
        [guide for _, guides, _ in cells for guide in guides],
        dtype=np.int64,
    )
    observed_umi = np.asarray(
        [count for _, _, counts in cells for count in counts],
        dtype=np.int64,
    )
    return guide_calling.SparseGuideBatch(
        cell_indices=cell_indices,
        indptr=indptr,
        guide_indices=guide_indices,
        observed_umi=observed_umi,
    )


class SparseGuideCallingTests(unittest.TestCase):
    def test_randomized_sparse_q_values_equal_dense_bh(self) -> None:
        rng = np.random.default_rng(107_2026)
        for hypothesis_count in (1, 2, 7, 31, 127):
            for _ in range(100):
                raw = rng.gamma(shape=0.7, scale=1.0, size=hypothesis_count)
                zero = rng.random(hypothesis_count) < 0.15
                raw[zero] = 0.0
                if not np.any(raw > 0):
                    raw[0] = 1.0
                probability = raw / raw.sum()

                observed_pool = np.flatnonzero(probability > 0)
                observed_count = int(rng.integers(0, min(observed_pool.size, 20) + 1))
                if observed_count:
                    guides = rng.choice(
                        observed_pool,
                        size=observed_count,
                        replace=False,
                    )
                    counts = rng.integers(1, 20, size=observed_count)
                    permutation = rng.permutation(observed_count)
                    guides = guides[permutation]
                    counts = counts[permutation]
                else:
                    guides = np.empty(0, dtype=np.int64)
                    counts = np.empty(0, dtype=np.int64)
                total = int(counts.sum())

                sparse_q = guide_calling.sparse_observed_q_values(
                    guides,
                    counts,
                    total_guide_umi=total,
                    channel_probability=probability,
                    hypothesis_count=hypothesis_count,
                )
                dense_p = np.ones(hypothesis_count, dtype=np.float64)
                if guides.size:
                    dense_p[guides] = binom.sf(
                        counts - 1,
                        total,
                        probability[guides],
                    )
                dense_q = dense_bh_canonical(dense_p)
                np.testing.assert_allclose(
                    sparse_q,
                    dense_q[guides],
                    rtol=0,
                    atol=0,
                )

    def test_stable_canonical_ties_match_dense_reference(self) -> None:
        hypothesis_count = 10
        probability = np.full(hypothesis_count, 0.1)
        guides = np.array([8, 2, 5, 0], dtype=np.int64)
        counts = np.array([2, 2, 2, 2], dtype=np.int64)
        total = int(counts.sum())

        observed = guide_calling.sparse_observed_q_values(
            guides,
            counts,
            total_guide_umi=total,
            channel_probability=probability,
            hypothesis_count=hypothesis_count,
        )
        dense_p = np.ones(hypothesis_count)
        dense_p[guides] = binom.sf(counts - 1, total, probability[guides])
        expected = dense_bh_canonical(dense_p)
        np.testing.assert_array_equal(observed, expected[guides])

    def test_exact_83401_hypothesis_denominator_matches_dense_bh(self) -> None:
        hypothesis_count = guide_calling.GUIDE_HYPOTHESES
        probability = np.full(hypothesis_count, 0.4 / (hypothesis_count - 3))
        probability[[2, 40, 80_000]] = [0.1, 0.2, 0.3]
        guides = np.array([80_000, 17, 40, 2], dtype=np.int64)
        counts = np.array([16, 8, 14, 12], dtype=np.int64)
        total = int(counts.sum())

        sparse_q = guide_calling.sparse_observed_q_values(
            guides,
            counts,
            total_guide_umi=total,
            channel_probability=probability,
            hypothesis_count=hypothesis_count,
        )
        dense_p = np.ones(hypothesis_count)
        dense_p[guides] = binom.sf(counts - 1, total, probability[guides])
        dense_q = dense_bh_canonical(dense_p)
        np.testing.assert_array_equal(sparse_q, dense_q[guides])

    def test_channel_null_uses_all_retained_hash_singlets(self) -> None:
        batches = [
            make_batch([(2, [0, 3], [2, 1]), (5, [], [])]),
            make_batch([(11, [1, 3], [4, 2])]),
        ]
        fitted = guide_calling.fit_channel_guide_null(
            iter(batches),
            hypothesis_count=5,
        )
        self.assertEqual(fitted.retained_cell_count, 3)
        self.assertEqual(fitted.total_guide_umi, 9)
        np.testing.assert_array_equal(fitted.guide_umi, [2, 4, 0, 3, 0])
        np.testing.assert_array_equal(
            fitted.probability,
            np.array([2, 4, 0, 3, 0], dtype=np.float64) / 9,
        )

    def test_streaming_calls_are_chunk_invariant(self) -> None:
        cells = [
            (1, [0, 1], [8, 1]),
            (4, [0, 2], [1, 9]),
            (8, [], []),
            (10, [3, 4], [7, 2]),
        ]
        one_batch = [make_batch(cells)]
        split_batches = [make_batch(cells[:2]), make_batch(cells[2:])]
        first_null = guide_calling.fit_channel_guide_null(
            one_batch,
            hypothesis_count=6,
        )
        second_null = guide_calling.fit_channel_guide_null(
            split_batches,
            hypothesis_count=6,
        )
        np.testing.assert_array_equal(first_null.guide_umi, second_null.guide_umi)
        np.testing.assert_array_equal(
            first_null.probability,
            second_null.probability,
        )

        target_map = np.array([0, -1, 1, 2, 2, -1], dtype=np.int64)
        ntc = np.array([False, True, False, False, False, True])
        first = list(
            guide_calling.iter_channel_guide_calls(
                one_batch,
                channel_null=first_null,
                guide_to_target=target_map,
                non_targeting_guides=ntc,
                fdr=1.0,
            )
        )
        second = list(
            guide_calling.iter_channel_guide_calls(
                split_batches,
                channel_null=second_null,
                guide_to_target=target_map,
                non_targeting_guides=ntc,
                fdr=1.0,
            )
        )
        self.assertEqual([call.cell_index for call in first], [1, 4, 8, 10])
        for left, right in zip(first, second, strict=True):
            self.assertEqual(left.cell_index, right.cell_index)
            self.assertEqual(left.total_guide_umi, right.total_guide_umi)
            self.assertEqual(
                left.total_called_guide_burden,
                right.total_called_guide_burden,
            )
            self.assertEqual(
                left.non_targeting_guide_burden,
                right.non_targeting_guide_burden,
            )
            np.testing.assert_array_equal(
                left.called_guide_indices,
                right.called_guide_indices,
            )
            np.testing.assert_array_equal(left.called_q_values, right.called_q_values)
            np.testing.assert_array_equal(left.target_indices, right.target_indices)

    def test_all_called_guides_burden_but_only_unique_targets_collapse(self) -> None:
        hypothesis_count = 7
        # Guide 0 is NTC, guide 1 is composite/unmapped, guides 2 and 3 map
        # to the same target, and guide 4 maps uniquely to another target.
        target_map = np.array([-1, -1, 4, 4, 9, -1, -1], dtype=np.int64)
        ntc = np.array([True, False, False, False, False, False, False])
        probability = np.full(hypothesis_count, 1 / hypothesis_count)
        call = guide_calling.call_sparse_cell(
            19,
            np.array([4, 1, 3, 0, 2]),
            np.array([5, 5, 5, 5, 5]),
            channel_probability=probability,
            guide_to_target=target_map,
            non_targeting_guides=ntc,
            hypothesis_count=hypothesis_count,
            fdr=1.0,
        )
        np.testing.assert_array_equal(call.called_guide_indices, [0, 1, 2, 3, 4])
        self.assertEqual(call.total_called_guide_burden, 5)
        self.assertEqual(call.non_targeting_guide_burden, 1)
        np.testing.assert_array_equal(call.target_indices, [4, 9])

    def test_q_threshold_is_inclusive_and_zero_observations_never_call(self) -> None:
        hypothesis_count = 5
        target_map = np.full(hypothesis_count, -1, dtype=np.int64)
        ntc = np.zeros(hypothesis_count, dtype=bool)
        probability = np.full(hypothesis_count, 0.2)
        q_values = guide_calling.sparse_observed_q_values(
            np.array([1]),
            np.array([1]),
            total_guide_umi=1,
            channel_probability=probability,
            hypothesis_count=hypothesis_count,
        )
        call = guide_calling.call_sparse_cell(
            0,
            np.array([1]),
            np.array([1]),
            channel_probability=probability,
            guide_to_target=target_map,
            non_targeting_guides=ntc,
            hypothesis_count=hypothesis_count,
            fdr=float(q_values[0]),
        )
        np.testing.assert_array_equal(call.called_guide_indices, [1])
        self.assertNotIn(0, call.called_guide_indices)
        self.assertNotIn(2, call.called_guide_indices)

    def test_sparse_batch_rejects_invalid_records(self) -> None:
        with self.assertRaisesRegex(
            guide_calling.GuideCallingError,
            "unique within a cell",
        ):
            make_batch([(0, [1, 1], [2, 3])])
        with self.assertRaisesRegex(
            guide_calling.GuideCallingError,
            "strictly increasing",
        ):
            make_batch([(3, [], []), (2, [], [])])
        with self.assertRaisesRegex(
            guide_calling.GuideCallingError,
            "positive",
        ):
            make_batch([(0, [1], [0])])

    def test_channel_null_hard_failures(self) -> None:
        with self.assertRaisesRegex(
            guide_calling.GuideCallingError,
            "no retained",
        ):
            guide_calling.fit_channel_guide_null([], hypothesis_count=3)
        with self.assertRaisesRegex(
            guide_calling.GuideCallingError,
            "no guide observations",
        ):
            guide_calling.fit_channel_guide_null(
                [make_batch([(0, [], [])])],
                hypothesis_count=3,
            )
        with self.assertRaisesRegex(
            guide_calling.GuideCallingError,
            "canonical guide axis",
        ):
            guide_calling.fit_channel_guide_null(
                [make_batch([(0, [3], [1])])],
                hypothesis_count=3,
            )
        with self.assertRaisesRegex(
            guide_calling.GuideCallingError,
            "globally strictly increasing",
        ):
            guide_calling.fit_channel_guide_null(
                [
                    make_batch([(2, [0], [1])]),
                    make_batch([(2, [0], [1])]),
                ],
                hypothesis_count=3,
            )

    def test_positive_observation_with_zero_channel_probability_is_error(self) -> None:
        with self.assertRaisesRegex(
            guide_calling.GuideCallingError,
            "zero channel probability",
        ):
            guide_calling.sparse_observed_q_values(
                np.array([1]),
                np.array([2]),
                total_guide_umi=2,
                channel_probability=np.array([1.0, 0.0, 0.0]),
                hypothesis_count=3,
            )

    def test_call_rejects_bad_mapping_and_probability_contracts(self) -> None:
        guides = np.array([0])
        counts = np.array([1])
        probability = np.array([0.5, 0.5])
        with self.assertRaisesRegex(
            guide_calling.GuideCallingError,
            "guide_to_target axis",
        ):
            guide_calling.call_sparse_cell(
                0,
                guides,
                counts,
                channel_probability=probability,
                guide_to_target=np.array([-1]),
                non_targeting_guides=np.array([False, False]),
                hypothesis_count=2,
            )
        with self.assertRaisesRegex(
            guide_calling.GuideCallingError,
            "boolean canonical",
        ):
            guide_calling.call_sparse_cell(
                0,
                guides,
                counts,
                channel_probability=probability,
                guide_to_target=np.array([-1, -1]),
                non_targeting_guides=np.array([0, 0]),
                hypothesis_count=2,
            )
        with self.assertRaisesRegex(
            guide_calling.GuideCallingError,
            "sum to one",
        ):
            guide_calling.sparse_observed_q_values(
                guides,
                counts,
                total_guide_umi=1,
                channel_probability=np.array([0.4, 0.4]),
                hypothesis_count=2,
            )

    def test_canonical_axis_requires_all_250_ntc_records(self) -> None:
        hypothesis_count = guide_calling.GUIDE_HYPOTHESES
        probability = np.full(hypothesis_count, 1 / hypothesis_count)
        with self.assertRaisesRegex(
            guide_calling.GuideCallingError,
            "exactly 250 non-targeting",
        ):
            guide_calling.call_sparse_cell(
                0,
                np.array([0]),
                np.array([1]),
                channel_probability=probability,
                guide_to_target=np.full(hypothesis_count, -1),
                non_targeting_guides=np.zeros(hypothesis_count, dtype=bool),
            )

    def test_second_pass_must_match_null_fit_pass(self) -> None:
        fitted = guide_calling.fit_channel_guide_null(
            [make_batch([(0, [0], [2]), (1, [1], [2])])],
            hypothesis_count=3,
        )
        calls = guide_calling.iter_channel_guide_calls(
            [make_batch([(0, [0], [2]), (1, [1], [3])])],
            channel_null=fitted,
            guide_to_target=np.array([-1, -1, -1]),
            non_targeting_guides=np.zeros(3, dtype=bool),
            fdr=1.0,
        )
        with self.assertRaisesRegex(
            guide_calling.GuideCallingError,
            "differ from channel null-fit pass",
        ):
            list(calls)


if __name__ == "__main__":
    unittest.main()
