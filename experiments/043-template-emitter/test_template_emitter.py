"""Tiny synthetic tests; deliberately no cell-eval2 invocation or real-data load."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))

from template_emitter import (  # noqa: E402
    arithmetic_mean_cpm,
    emit_template_fixed_margins,
    integer_column_totals,
    scorer_pseudobulk,
    template_scaled_cpm_target,
)


class TemplateEmitterTest(unittest.TestCase):
    def setUp(self) -> None:
        # Variable-depth, sparse cells with gene-specific depth correlation.
        self.template = sp.csr_matrix(
            np.array(
                [
                    [8, 0, 2, 0, 0, 0],
                    [0, 8, 0, 4, 0, 0],
                    [8, 0, 0, 0, 8, 0],
                    [0, 5, 5, 0, 0, 10],
                    [0, 0, 8, 8, 8, 0],
                    [12, 0, 0, 6, 0, 12],
                ],
                dtype=np.int64,
            )
        )

    def test_control_profile_is_an_exact_fixed_point(self) -> None:
        profile = np.asarray(self.template.sum(axis=0)).ravel()
        emitted, diag = emit_template_fixed_margins(
            self.template, profile=profile, seed=7
        )
        np.testing.assert_array_equal(emitted.toarray(), self.template.toarray())
        self.assertEqual(diag.moved_counts, 0)
        self.assertEqual(diag.introduced_nnz, 0)
        self.assertEqual(diag.removed_template_nnz, 0)
        self.assertEqual(diag.template_cpm_relative_l1, 0.0)

    def test_perturbed_profile_has_exact_margins_and_pseudobulk(self) -> None:
        base = np.asarray(self.template.sum(axis=0)).ravel().astype(float)
        profile = base * np.array([1.35, 0.75, 1.20, 0.60, 1.10, 0.85])
        rows = np.asarray(self.template.sum(axis=1)).ravel().astype(np.int64)
        cols = integer_column_totals(profile, int(rows.sum()), seed=12)

        emitted, diag = emit_template_fixed_margins(
            self.template, column_totals=cols, row_margins=rows, seed=11
        )

        np.testing.assert_array_equal(
            np.asarray(emitted.sum(axis=1)).ravel(), rows
        )
        np.testing.assert_array_equal(
            np.asarray(emitted.sum(axis=0)).ravel(), cols
        )
        expected_pb = np.log1p(50_000.0 * cols / cols.sum())
        np.testing.assert_array_equal(scorer_pseudobulk(emitted), expected_pb)
        self.assertTrue(np.issubdtype(emitted.dtype, np.integer))
        self.assertGreaterEqual(emitted.data.min(initial=0), 0)
        self.assertEqual(diag.moved_counts, diag.row_l1_before_repair // 2)
        self.assertLessEqual(diag.introduced_nnz, diag.moved_counts)

        # The output has valid arithmetic CPM means and exactly the template's
        # library-size distribution, unlike the old constant-depth even emitter.
        cpm = arithmetic_mean_cpm(emitted)
        self.assertTrue(np.all(np.isfinite(cpm)))
        np.testing.assert_array_equal(
            np.sort(np.asarray(emitted.sum(axis=1)).ravel()), np.sort(rows)
        )

        # On this deliberately depth-confounded fixture, preserving the control
        # template is materially closer to its depth-aware CPM target than evenly
        # spreading each gene over constant-depth synthetic cells.
        target_cpm = template_scaled_cpm_target(self.template, cols, rows)
        even = np.zeros(self.template.shape, dtype=np.int64)
        for gene, total in enumerate(cols):
            base_count, remainder = divmod(int(total), self.template.shape[0])
            even[:, gene] = base_count
            even[:remainder, gene] += 1
        even_cpm = arithmetic_mean_cpm(sp.csr_matrix(even))
        even_rel_l1 = np.abs(even_cpm - target_cpm).sum() / np.abs(target_cpm).sum()
        self.assertLess(diag.template_cpm_relative_l1, 0.02)
        self.assertLess(diag.template_cpm_relative_l1, even_rel_l1)

    def test_random_sparse_margin_invariants(self) -> None:
        rng = np.random.default_rng(20260822)
        for case in range(50):
            n_cells = int(rng.integers(2, 10))
            n_genes = int(rng.integers(2, 25))
            template = rng.poisson(2.0, size=(n_cells, n_genes)).astype(np.int64)
            template[rng.random(template.shape) < 0.72] = 0
            for row in np.flatnonzero(template.sum(axis=1) == 0):
                template[row, rng.integers(n_genes)] = int(rng.integers(1, 8))
            profile = rng.lognormal(0.0, 1.0, n_genes)
            profile[rng.random(n_genes) < 0.25] = 0
            if profile.sum() == 0:
                profile[0] = 1

            emitted, _ = emit_template_fixed_margins(
                sp.csr_matrix(template), profile=profile, seed=case
            )
            np.testing.assert_array_equal(
                np.asarray(emitted.sum(axis=1)).ravel(), template.sum(axis=1)
            )
            self.assertEqual(int(emitted.sum()), int(template.sum()))
            self.assertGreaterEqual(emitted.data.min(initial=0), 0)

    def test_aligned_template_can_use_full_axis_library_sizes(self) -> None:
        aligned = self.template[:, :5]
        full_libs = np.asarray(self.template.sum(axis=1)).ravel().astype(np.int64)
        cols = np.array([21, 9, 14, 7, 13], dtype=np.int64)
        rows = integer_column_totals(full_libs, int(cols.sum()), seed=32)
        emitted, _ = emit_template_fixed_margins(
            aligned,
            column_totals=cols,
            row_margins=rows,
            template_library_sizes=full_libs,
            seed=31,
        )
        np.testing.assert_array_equal(
            np.asarray(emitted.sum(axis=1)).ravel(), rows
        )
        np.testing.assert_array_equal(
            np.asarray(emitted.sum(axis=0)).ravel(), cols
        )


if __name__ == "__main__":
    unittest.main()
