#!/usr/bin/env python3
"""Synthetic tests for the Experiment 107 source-model freeze boundary."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from contextlib import ExitStack
from pathlib import Path
from unittest import mock

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import fit_source  # noqa: E402
from test_model_input import SyntheticModelInputFixture, write_json  # noqa: E402


class SourceFitTests(unittest.TestCase):
    def test_binary_csr_preserves_sparse_shape_and_empty_rows(self) -> None:
        matrix = fit_source._binary_csr(
            np.array([0, 1, 1, 3]),
            np.array([2, 0, 3]),
            rows=3,
            columns=4,
            name="test",
        )
        self.assertEqual(matrix.shape, (3, 4))
        self.assertEqual(matrix.nnz, 3)
        np.testing.assert_array_equal(
            matrix.toarray(),
            [
                [0, 0, 1, 0],
                [0, 0, 0, 0],
                [1, 0, 0, 1],
            ],
        )

    def test_solver_report_requires_convergence_record(self) -> None:
        with self.assertRaisesRegex(fit_source.SourceFitError, "diagnostics"):
            fit_source._solver_report(None)
        diagnostics = fit_source.modeling.RidgeSolveDiagnostics(
            ridge_lambda=3.0,
            iterations=7,
            operator_applications=9,
            warm_started=True,
            converged=np.array([True, True]),
            rhs_norms=np.array([2.0, 1.0]),
            residual_norms=np.array([1e-10, 2e-10]),
            relative_residual_norms=np.array([5e-11, 2e-10]),
            absolute_tolerance=1e-10,
            relative_tolerance=1e-8,
        )
        observed = fit_source._solver_report(diagnostics)
        self.assertTrue(observed["all_components_converged"])
        self.assertEqual(observed["iterations"], 7)
        self.assertEqual(observed["maximum_relative_residual"], 2e-10)

    def test_memory_guard_rejects_invalid_limit(self) -> None:
        for limit in (0.0, -1.0, fit_source.MAX_RSS_GIB + 1):
            with self.assertRaisesRegex(fit_source.SourceFitError, "memory limit"):
                fit_source.guard_memory("test", limit)

    def test_fold_diagnostics_pack_and_stability_are_recomputable(self) -> None:
        indptr, values = fit_source._pack_fold_vectors(
            (
                np.array([0.0, 2.0]),
                np.array([3.0]),
                np.array([], dtype=np.float64),
            )
        )
        np.testing.assert_array_equal(indptr, [0, 2, 3, 3])
        np.testing.assert_array_equal(values, [0.0, 2.0, 3.0])

        effects = np.zeros((4, 3, 4), dtype=np.float64)
        effects[:, 0] = np.array([1.0, -1.0, 2.0, -2.0])
        effects[:, 1] = np.array([1.0, 1.0, -1.0, -1.0])
        effects[1, 0] *= -1.0
        transferable = np.ones((4, 3), dtype=bool)
        transferable[:, 2] = False
        pairs, cosine, sign = fit_source._fold_stability(
            effects,
            transferable,
            np.array([True, True, False, False]),
        )
        self.assertEqual(pairs.shape, (6, 2))
        self.assertEqual(cosine.shape, (6, 3))
        self.assertEqual(sign.shape, (6, 3))
        self.assertTrue(np.all(np.isnan(cosine[:, 2])))
        self.assertTrue(np.all(np.isnan(sign[:, 2])))
        self.assertAlmostEqual(cosine[0, 0], -1.0)
        self.assertAlmostEqual(sign[0, 0], 0.0)

    def _fixture_patches(
        self,
        fixture: SyntheticModelInputFixture,
    ) -> tuple[mock._patch, ...]:
        return (
            mock.patch.object(fit_source, "WORK", fixture.work),
            mock.patch.object(fit_source.model_input, "WORK", fixture.work),
            mock.patch.object(
                fit_source,
                "EXPECTED_TARGET_COUNT",
                fixture.target_count,
            ),
            mock.patch.object(
                fit_source,
                "EXPECTED_GUIDE_COUNT",
                fixture.guide_count,
            ),
            mock.patch.object(
                fit_source,
                "EXPECTED_EXPORTABLE_RESPONSE_COUNT",
                fixture.exportable_count,
            ),
            mock.patch.object(
                fit_source,
                "EXPECTED_ORIGINAL_RESPONSE_COUNT",
                fixture.response_count,
            ),
            mock.patch.object(
                fit_source.model_input,
                "EXPECTED_VCC_MAPPED_RESPONSE_COUNT",
                fixture.vcc_mapped_response_count,
            ),
            mock.patch.object(
                fit_source.source_artifacts,
                "verify_source_artifacts",
                return_value=fixture.source_verification,
            ),
        )

    def test_production_shaped_build_load_and_freeze_path(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as directory:
            fixture = SyntheticModelInputFixture(Path(directory))
            fixture.build()
            patches = self._fixture_patches(fixture)
            with ExitStack() as stack:
                for patch in patches:
                    stack.enter_context(patch)
                loaded = fit_source.load_frozen_input(
                    fixture.output,
                    source_verification=fixture.source_verification,
                )
                self.assertEqual(loaded[1].shape, (128, fixture.target_count))
                self.assertEqual(
                    loaded[0].shape,
                    (128, fixture.exportable_count),
                )
                output = fixture.root / "frozen-source-model"
                manifest = fit_source.fit_and_freeze(
                    fixture.output,
                    output,
                )

            self.assertEqual(
                manifest["status"],
                "frozen_before_held_outcome_access",
            )
            self.assertEqual(
                manifest["active_memory_monitor"]["status"],
                "complete",
            )
            self.assertTrue(
                manifest["final_support"][
                    "full_per_target_fold_and_final_arrays_in_model_file"
                ]
            )
            self.assertTrue(
                manifest["projection_postconditions"]["validated_after_serialization"]
            )
            with np.load(output / "source-model.npz", allow_pickle=False) as values:
                for name in (
                    "fold_support_distinct_called_guides",
                    "fold_support_guides_with_three_untreated_cells",
                    "fold_support_untreated_carrier_cells",
                    "fold_support_training_channel_count",
                    "fold_support_supported",
                    "final_support_distinct_called_guides",
                    "final_support_supported",
                ):
                    self.assertIn(name, values.files)
                self.assertEqual(
                    values["fold_support_supported"].shape,
                    (
                        len(fit_source.modeling.FROZEN_HELD_CHANNEL_FOLDS),
                        fixture.target_count,
                    ),
                )
                np.testing.assert_array_equal(
                    values["final_support_supported"],
                    np.ones(fixture.target_count, dtype=bool),
                )
                self.assertEqual(
                    values["exported_effects"].shape,
                    (fixture.target_count, fixture.response_count),
                )
                for name in (
                    "exported_effects",
                    "coefficients_binary_units",
                    "projected_coefficients",
                    "nuisance_response_slopes",
                    "response_basis",
                    "fold_selected_effects",
                ):
                    np.testing.assert_array_equal(values[name][..., -1], 0.0)
                np.testing.assert_array_equal(
                    values["response_exportable"],
                    [True, True, True, False],
                )
            self.assertEqual(
                manifest["effect_semantics"],
                (
                    "identifiable_high_moi_conditional_target_contrast_"
                    "minimum_norm_row_space_component_not_singleton_marginal_effect"
                ),
            )
            self.assertTrue(manifest["design_audit_validated_before_response_fit"])
            self.assertEqual(
                manifest["modeled_response_genes"],
                fixture.exportable_count,
            )
            self.assertEqual(manifest["response_genes"], fixture.response_count)

    def test_forbidden_response_changes_cannot_affect_selection_or_effects(
        self,
    ) -> None:
        with (
            tempfile.TemporaryDirectory(dir=HERE) as first_directory,
            tempfile.TemporaryDirectory(dir=HERE) as second_directory,
        ):
            first_fixture = SyntheticModelInputFixture(Path(first_directory))
            second_fixture = SyntheticModelInputFixture(
                Path(second_directory),
                forbidden_response_shift=100,
            )
            first_fixture.build()
            second_fixture.build()

            results = []
            for fixture in (first_fixture, second_fixture):
                with ExitStack() as stack:
                    for patch in self._fixture_patches(fixture):
                        stack.enter_context(patch)
                    loaded = fit_source.load_frozen_input(
                        fixture.output,
                        source_verification=fixture.source_verification,
                    )
                    results.append(
                        fit_source.modeling.select_production_estimator(
                            loaded[0],
                            loaded[1],
                            loaded[2],
                            loaded[3],
                            loaded[4],
                            design_audit=loaded[7],
                        )
                    )

            self.assertEqual(results[0].selected_lambda, results[1].selected_lambda)
            np.testing.assert_array_equal(
                results[0].estimator.exported_effects,
                results[1].estimator.exported_effects,
            )
            np.testing.assert_array_equal(
                results[0].fold_mean_gains,
                results[1].fold_mean_gains,
            )

    def test_load_frozen_input_rejects_stale_assignment_summary(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as directory:
            fixture = SyntheticModelInputFixture(Path(directory))
            fixture.build()
            summary_path = fixture.work / "assignment-summary.json"
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            summary["retained_cells"] = 999
            write_json(summary_path, summary)
            patches = self._fixture_patches(fixture)
            with ExitStack() as stack:
                for patch in patches:
                    stack.enter_context(patch)
                with self.assertRaisesRegex(
                    fit_source.SourceFitError,
                    "assignment summary is stale",
                ):
                    fit_source.load_frozen_input(
                        fixture.output,
                        source_verification=fixture.source_verification,
                    )

    def test_load_frozen_input_rejects_stale_response_report(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as directory:
            fixture = SyntheticModelInputFixture(Path(directory))
            fixture.build()
            response_report_path = fixture.work / "channel05" / "response-report.json"
            response_report = json.loads(
                response_report_path.read_text(encoding="utf-8")
            )
            response_report["changed"] = True
            write_json(response_report_path, response_report)
            patches = self._fixture_patches(fixture)
            with ExitStack() as stack:
                for patch in patches:
                    stack.enter_context(patch)
                with self.assertRaisesRegex(
                    fit_source.SourceFitError,
                    "response report changed after model input",
                ):
                    fit_source.load_frozen_input(
                        fixture.output,
                        source_verification=fixture.source_verification,
                    )

    def test_load_frozen_input_rejects_coherently_rehashed_eligibility_drift(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as directory:
            fixture = SyntheticModelInputFixture(Path(directory))
            fixture.build()
            with np.load(fixture.output, allow_pickle=False) as values:
                arrays = {name: values[name] for name in values.files}
            arrays["precomputed_eligible"] = arrays["precomputed_eligible"].copy()
            arrays["precomputed_eligible"][0] = False
            np.savez_compressed(fixture.output, **arrays)
            report_path = fixture.output.with_suffix(".json")
            report = json.loads(report_path.read_text(encoding="utf-8"))
            report["output_sha256"] = fit_source.sha256_file(fixture.output)
            write_json(report_path, report)

            patches = self._fixture_patches(fixture)
            with ExitStack() as stack:
                for patch in patches:
                    stack.enter_context(patch)
                with self.assertRaisesRegex(
                    fit_source.SourceFitError,
                    "frozen untreated target support changed",
                ):
                    fit_source.load_frozen_input(
                        fixture.output,
                        source_verification=fixture.source_verification,
                    )


if __name__ == "__main__":
    unittest.main()
