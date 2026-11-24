from __future__ import annotations

import math
import io
import itertools
import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.stats import binom


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import preprocess  # noqa: E402


class PreprocessTests(unittest.TestCase):
    @staticmethod
    def matrix_stream(entries: list[tuple[int, int, int]]) -> io.BytesIO:
        lines = [
            "%%MatrixMarket matrix coordinate integer general",
            "% synthetic",
            f"3 4 {len(entries)}",
            *(f"{row} {column} {value}" for row, column, value in entries),
        ]
        return io.BytesIO(("\n".join(lines) + "\n").encode())

    def test_streaming_transcriptome_passes(self) -> None:
        entries = [(1, 1, 2), (3, 1, 1), (2, 2, 5), (1, 4, 3)]
        summary = preprocess.scan_transcriptome_summary(
            self.matrix_stream(entries),
            mitochondrial_features=np.array([2]),
        )
        self.assertEqual(summary.coordinate_order, "column_major")
        np.testing.assert_array_equal(summary.totals, [3, 5, 0, 3])
        np.testing.assert_array_equal(summary.detected_features, [2, 1, 0, 1])
        np.testing.assert_array_equal(summary.mitochondrial_counts, [1, 0, 0, 0])
        np.testing.assert_array_equal(summary.feature_totals, [5, 5, 1])

        ambient = preprocess.ambient_feature_counts(
            self.matrix_stream(entries),
            totals=summary.totals,
            lower=3,
        )
        np.testing.assert_array_equal(ambient, [5, 0, 1])

        probability = np.array([0.5, 0.25, 0.25])
        tested = np.array([False, True, False, True])
        observed = preprocess.observed_data_log_probabilities(
            self.matrix_stream(entries),
            totals=summary.totals,
            ambient_probability=probability,
            tested=tested,
        )
        expected = np.array(
            [
                5 * np.log(0.25) - math.lgamma(6),
                3 * np.log(0.5) - math.lgamma(4),
            ]
        )
        np.testing.assert_allclose(observed, expected, rtol=0, atol=1e-14)

    def test_streaming_transcriptome_rejects_duplicate_or_unsorted(self) -> None:
        entries = [(1, 1, 2), (1, 1, 1), (2, 2, 5)]
        with self.assertRaisesRegex(preprocess.PreprocessError, "canonical order"):
            preprocess.scan_transcriptome_summary(self.matrix_stream(entries))

    def test_good_turing_is_positive_and_normalized(self) -> None:
        counts = np.array([0, 0, 1, 2, 2, 4, 10, 10, 10], dtype=np.int64)
        proportions = preprocess.safe_good_turing_proportions(counts)
        self.assertTrue(np.all(proportions > 0))
        np.testing.assert_allclose(proportions.sum(), 1.0, rtol=0, atol=1e-15)
        self.assertEqual(proportions[0], proportions[1])
        self.assertNotEqual(proportions[2], proportions[3])

    def test_observed_log_probability_allows_globally_zero_features(self) -> None:
        entries = [(1, 1, 2), (2, 2, 3), (1, 3, 1), (2, 3, 1)]
        totals = np.array([2, 3, 2, 0], dtype=np.int64)
        probability = np.array([0.4, 0.6, 0.0], dtype=np.float64)
        observed = preprocess.observed_data_log_probabilities(
            self.matrix_stream(entries),
            totals=totals,
            ambient_probability=probability,
            tested=np.array([True, True, True, False]),
        )
        self.assertEqual(observed.shape, (3,))

    def test_multinomial_log_probability_matches_scipy(self) -> None:
        probability = np.array([0.1, 0.2, 0.3, 0.4], dtype=np.float64)
        indices = np.array([0, 2, 3], dtype=np.int64)
        counts = np.array([2, 1, 3], dtype=np.int64)
        observed = preprocess.multinomial_log_probability(
            indices,
            counts,
            6,
            probability,
        )
        expected = np.log(
            math.factorial(6)
            / (
                math.factorial(2)
                * math.factorial(0)
                * math.factorial(1)
                * math.factorial(3)
            )
            * 0.1**2
            * 0.3
            * 0.4**3
        )
        np.testing.assert_allclose(observed, expected, rtol=0, atol=1e-14)

    def test_sparse_guide_bh_includes_all_hypotheses(self) -> None:
        hypothesis_count = 10
        indices = np.array([1, 4, 8], dtype=np.int64)
        counts = np.array([7, 2, 1], dtype=np.int64)
        probability = np.full(hypothesis_count, 0.1, dtype=np.float64)
        sparse = preprocess.sparse_guide_qvalues(
            indices,
            counts,
            total_guide_umi=10,
            channel_probability=probability,
            hypothesis_count=hypothesis_count,
        )
        dense_p = np.ones(hypothesis_count, dtype=np.float64)
        dense_p[indices] = binom.sf(counts - 1, 10, probability[indices])
        dense = preprocess.benjamini_hochberg(dense_p)
        np.testing.assert_allclose(sparse, dense[indices], rtol=0, atol=0)

    def test_emptydrops_monte_carlo_is_chunk_invariant(self) -> None:
        ambient = np.array([0.5, 0.5], dtype=np.float64)
        totals = np.array([2, 2, 3], dtype=np.int64)
        observed = np.array(
            [
                preprocess.multinomial_data_log_probability(
                    np.array([0]),
                    np.array([2]),
                    ambient,
                ),
                preprocess.multinomial_data_log_probability(
                    np.array([0, 1]),
                    np.array([1, 1]),
                    ambient,
                ),
                preprocess.multinomial_data_log_probability(
                    np.array([0, 1]),
                    np.array([1, 2]),
                    ambient,
                ),
            ]
        )
        first = preprocess.emptydrops_monte_carlo(
            totals,
            observed,
            ambient,
            channel_number=3,
            iterations=2_000,
            iteration_block_size=1,
        )
        second = preprocess.emptydrops_monte_carlo(
            totals,
            observed,
            ambient,
            channel_number=3,
            iterations=2_000,
            iteration_block_size=37,
        )
        np.testing.assert_array_equal(first.exceedances, second.exceedances)
        np.testing.assert_array_equal(first.p_value, second.p_value)
        np.testing.assert_array_equal(first.limited, second.limited)
        self.assertLess(abs(first.p_value[0] - 0.5), 0.04)
        self.assertEqual(first.p_value[1], 1.0)

    def test_emptydrops_monte_carlo_matches_independent_exact_enumeration(
        self,
    ) -> None:
        ambient = np.array([0.25, 0.75], dtype=np.float64)
        checkpoints = np.array([2, 3], dtype=np.int64)
        paths = np.array(
            list(itertools.product(range(ambient.size), repeat=checkpoints[-1])),
            dtype=np.int32,
        )

        # Independent count-form oracle for every possible categorical path.
        exact_scores = np.empty((paths.shape[0], checkpoints.size), dtype=np.float64)
        path_probability = np.empty(paths.shape[0], dtype=np.float64)
        for path_index, path in enumerate(paths):
            path_probability[path_index] = float(np.prod(ambient[path]))
            for checkpoint_index, total in enumerate(checkpoints):
                counts = np.bincount(path[:total], minlength=ambient.size)
                exact_scores[path_index, checkpoint_index] = float(
                    np.dot(counts, np.log(ambient))
                    - np.sum([math.lgamma(int(count) + 1) for count in counts])
                )

        implementation_scores = preprocess.simulated_prefix_log_probabilities_batch(
            paths,
            np.log(ambient),
            checkpoints,
        )
        np.testing.assert_allclose(
            implementation_scores,
            exact_scores,
            rtol=0,
            atol=1e-15,
        )
        np.testing.assert_allclose(path_probability.sum(), 1.0, rtol=0, atol=1e-15)

        observed = np.array(
            [
                2 * np.log(ambient[0]) - math.lgamma(3),
                3 * np.log(ambient[0]) - math.lgamma(4),
            ],
            dtype=np.float64,
        )
        exact_tail_probability = np.array(
            [
                path_probability[exact_scores[:, index] <= observed[index]].sum()
                for index in range(checkpoints.size)
            ]
        )
        np.testing.assert_allclose(
            exact_tail_probability,
            [ambient[0] ** 2, ambient[0] ** 3],
            rtol=0,
            atol=1e-15,
        )

        iterations = 257
        channel = 6
        result = preprocess.emptydrops_monte_carlo(
            checkpoints,
            observed,
            ambient,
            channel_number=channel,
            iterations=iterations,
            iteration_block_size=17,
        )

        # Reconstruct the frozen independent PCG64 streams, but evaluate each
        # draw with the count-form oracle above rather than the production
        # prefix-likelihood implementation.
        expected_exceedances = np.zeros(checkpoints.size, dtype=np.int64)
        for iteration in range(iterations):
            seed = np.random.SeedSequence(
                [preprocess.EMPTYDROPS_BASE_SEED, channel, iteration]
            )
            generator = np.random.Generator(np.random.PCG64(seed))
            path = generator.choice(
                ambient.size,
                size=int(checkpoints[-1]),
                p=ambient,
            )
            for index, total in enumerate(checkpoints):
                counts = np.bincount(path[:total], minlength=ambient.size)
                score = float(
                    np.dot(counts, np.log(ambient))
                    - np.sum([math.lgamma(int(count) + 1) for count in counts])
                )
                expected_exceedances[index] += score <= observed[index]

        np.testing.assert_array_equal(result.exceedances, expected_exceedances)
        np.testing.assert_allclose(
            result.p_value,
            (expected_exceedances + 1) / (iterations + 1),
            rtol=0,
            atol=0,
        )

    def test_batched_prefix_probabilities_match_scalar_reference(self) -> None:
        choices = np.array(
            [
                [0, 1, 1, 2, 0, 2],
                [2, 2, 1, 0, 1, 2],
                [1, 0, 1, 0, 2, 1],
            ],
            dtype=np.int32,
        )
        log_probability = np.log(np.array([0.2, 0.3, 0.5]))
        checkpoints = np.array([2, 4, 6], dtype=np.int64)
        expected = np.stack(
            [
                preprocess._simulated_prefix_log_probabilities(
                    row,
                    log_probability,
                    checkpoints,
                )
                for row in choices
            ]
        )
        observed = preprocess.simulated_prefix_log_probabilities_batch(
            choices,
            log_probability,
            checkpoints,
        )
        np.testing.assert_allclose(observed, expected, rtol=0, atol=1e-15)

    def test_emptydrops_forced_retain_precedes_bh(self) -> None:
        totals = np.array([50, 101, 150, 1_000], dtype=np.int64)
        p_values = np.array([np.nan, 0.03, 0.02, 0.9], dtype=np.float64)
        fdr = preprocess.emptydrops_fdr(
            totals,
            p_values,
            lower=100,
            retain=1_000,
        )
        self.assertTrue(np.isnan(fdr[0]))
        expected = preprocess.benjamini_hochberg(np.array([0.03, 0.02, 0.0]))
        np.testing.assert_allclose(fdr[1:], expected, rtol=0, atol=0)

    def test_demultiplex_rounds_peel_only_negatives(self) -> None:
        seen: list[tuple[tuple[str, ...], int]] = []

        def classify(active: tuple[str, ...], round_number: int) -> dict[str, str]:
            seen.append((active, round_number))
            if round_number == 1:
                return {"a": "Negative", "b": "Doublet", "c": "hash1"}
            if round_number == 2:
                return {"b": "Negative", "c": "hash1"}
            return {"c": "hash1"}

        calls = preprocess.run_demultiplex_rounds(("a", "b", "c"), classify)
        self.assertEqual(
            seen,
            [(("a", "b", "c"), 1), (("b", "c"), 2), (("c",), 3)],
        )
        self.assertEqual(
            calls,
            {"a": "Negative", "b": "Negative", "c": "hash1"},
        )

    def test_demultiplex_rounds_stops_when_no_new_negative(self) -> None:
        rounds: list[int] = []

        def classify(active: tuple[str, ...], round_number: int) -> dict[str, str]:
            rounds.append(round_number)
            return {cell: "hash1" for cell in active}

        calls = preprocess.run_demultiplex_rounds(("a", "b"), classify)
        self.assertEqual(rounds, [1])
        self.assertEqual(calls, {"a": "hash1", "b": "hash1"})

    def test_within_group_demean_has_zero_group_means(self) -> None:
        values = np.array(
            [[1.0, 2.0], [3.0, 6.0], [10.0, 4.0], [14.0, 8.0]],
            dtype=np.float64,
        )
        groups = np.array([1, 1, 2, 2])
        centered = preprocess.within_group_demean(values, groups)
        np.testing.assert_allclose(
            centered[groups == 1].mean(axis=0),
            0.0,
            atol=1e-15,
        )
        np.testing.assert_allclose(
            centered[groups == 2].mean(axis=0),
            0.0,
            atol=1e-15,
        )

    def test_nuisance_residualizes_response_and_target_design(self) -> None:
        nuisance_train = np.array([[1.0], [2.0], [3.0], [4.0]])
        response_train = np.c_[
            2 * nuisance_train[:, 0],
            [1.0, -1.0, 1.0, -1.0],
        ]
        target_train = np.c_[3 * nuisance_train[:, 0], [0.0, 1.0, 0.0, 1.0]]
        nuisance_held = np.array([[5.0], [6.0]])
        response_held = np.c_[2 * nuisance_held[:, 0], [2.0, -2.0]]
        target_held = np.c_[3 * nuisance_held[:, 0], [1.0, 0.0]]

        response_resid, target_resid, held_response, held_target = (
            preprocess.residualize_training_and_held(
                nuisance_train,
                response_train,
                target_train,
                nuisance_held,
                response_held,
                target_held,
            )
        )
        np.testing.assert_allclose(response_resid[:, 0], 0.0, atol=1e-12)
        np.testing.assert_allclose(target_resid[:, 0], 0.0, atol=1e-12)
        np.testing.assert_allclose(held_response[:, 0], 0.0, atol=1e-12)
        np.testing.assert_allclose(held_target[:, 0], 0.0, atol=1e-12)

    def test_identifiable_projection_removes_nullspace_shift(self) -> None:
        design = np.array(
            [
                [1.0, 0.0, 1.0],
                [0.0, 1.0, 1.0],
                [1.0, 1.0, 2.0],
                [-1.0, -1.0, -2.0],
            ]
        )
        coefficients = np.array([[1.0], [2.0], [3.0]])
        null = np.array([[1.0], [1.0], [-1.0]])
        projected, _ = preprocess.project_identifiable_effect(design, coefficients)
        shifted, _ = preprocess.project_identifiable_effect(
            design,
            coefficients + 17.0 * null,
        )
        np.testing.assert_allclose(projected, shifted, rtol=0, atol=1e-12)

    def test_masked_blend_preserves_unsupported_bits(self) -> None:
        incumbent = np.array(
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            dtype=np.float64,
        )
        source = np.array(
            [[9.0, 8.0, 7.0], [6.0, 5.0, 4.0]],
            dtype=np.float64,
        )
        supported = np.array([[True, False, False], [False, True, False]])
        candidate = preprocess.masked_blend_exact_fallback(
            incumbent,
            source,
            supported,
            0.25,
        )
        self.assertEqual(
            candidate[0, 1:].tobytes(),
            incumbent[0, 1:].tobytes(),
        )
        self.assertEqual(
            candidate[1, [0, 2]].tobytes(),
            incumbent[1, [0, 2]].tobytes(),
        )
        np.testing.assert_allclose(candidate[supported], [3.0, 5.0])


if __name__ == "__main__":
    unittest.main()
