from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np


MODULE_PATH = Path(__file__).with_name("develop.py")
SPEC = importlib.util.spec_from_file_location("experiment142_develop", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
DEVELOP = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = DEVELOP
SPEC.loader.exec_module(DEVELOP)


class Experiment142DevelopmentTests(unittest.TestCase):
    def test_block_grouped_folds_are_balanced_and_keep_blocks_intact(self) -> None:
        sizes = [4, 4, 3, 3, 2, 2, 1, 1]
        blocks = np.concatenate(
            [np.full(size, block, dtype=np.int16) for block, size in enumerate(sizes)]
        )
        folds = DEVELOP.block_grouped_folds(
            blocks, fold_count=5, expected_sizes=(4, 4, 4, 4, 4)
        )
        self.assertEqual(
            [int(np.sum(folds == fold)) for fold in range(5)],
            [4, 4, 4, 4, 4],
        )
        for block in np.unique(blocks):
            self.assertEqual(len(np.unique(folds[blocks == block])), 1)
        np.testing.assert_array_equal(
            folds,
            DEVELOP.block_grouped_folds(
                blocks, fold_count=5, expected_sizes=(4, 4, 4, 4, 4)
            ),
        )

    def test_eta_zero_returns_exact_frozen_baseline(self) -> None:
        slots = np.arange(24, dtype=np.float32).reshape(3, 2, 4)
        baseline = np.mean(slots.astype(np.float64), axis=1)
        weights = np.full((3, 2), 0.5)
        actual = DEVELOP.mix_slots(slots, weights, 0.0, baseline)
        np.testing.assert_array_equal(actual, baseline)

    def test_pm_replacement_preserves_missing_targets_and_genes(self) -> None:
        centered = np.asarray(
            [
                [2.0, 4.0, 6.0, 8.0],
                [1.0, 3.0, 5.0, 7.0],
            ]
        )
        pm = {
            "target_columns": np.asarray([0, 2], dtype=np.int64),
            "source_columns": np.asarray([1, 0], dtype=np.int64),
            "mean": np.asarray([10.0, 20.0]),
            "scale": 2.0,
            "delta": np.asarray([[14.0, 26.0], [99.0, 99.0]]),
            "pert_index": {"both": 0, "pm_only": 1},
        }
        actual, report = DEVELOP.composite_k562_standardized_rows(
            centered,
            ("both", "k_only"),
            ("A", "B", "C", "D"),
            ("A", "B", "C", "D"),
            2.0,
            pm,
        )
        base = centered / 2.0
        # The helper's mean vector is aligned with target/source column pairs.
        pm_standardized = np.asarray([(26.0 - 10.0) / 2.0, (14.0 - 20.0) / 2.0])
        expected = base.copy()
        expected[0, [0, 2]] = 0.5 * base[0, [0, 2]] + 0.5 * pm_standardized
        np.testing.assert_allclose(actual, expected)
        np.testing.assert_array_equal(actual[0, [1, 3]], base[0, [1, 3]])
        np.testing.assert_array_equal(actual[1], base[1])
        self.assertEqual(report["targets_modified"], 1)
        self.assertEqual(report["targets_preserved_exactly"], 1)

    def test_target_pairing_null_weights_use_same_draw_across_donors(self) -> None:
        null_z = np.asarray(
            [
                [[-1.0, 1.0], [0.0, 0.0]],
                [[1.0, -1.0], [0.5, -0.5]],
            ]
        )
        with patch.object(DEVELOP, "ETAS", (0.0, 0.5)):
            weights = DEVELOP.target_pairing_null_weight_grid(null_z)
        np.testing.assert_array_equal(weights[:, 0], 0.5)
        self.assertGreater(weights[0, 1, 0, 1], weights[0, 1, 0, 0])
        self.assertGreater(weights[1, 1, 0, 0], weights[1, 1, 0, 1])

    def test_false_pair_reliability_is_recomputed_before_weighting(self) -> None:
        half = np.asarray(
            [
                [
                    [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
                    [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
                ]
            ],
            dtype=np.float32,
        )
        valid = np.ones((1, 3), dtype=bool)
        permutations = np.asarray([[1, 2, 0], [2, 0, 1]], dtype=np.int16)
        with (
            patch.object(DEVELOP, "NULL_DRAWS", 2),
            patch.object(DEVELOP, "MIN_DIRECTIONAL_VALUES", 2),
        ):
            z = DEVELOP.recompute_target_null_z(
                half,
                valid,
                np.ones(2),
                np.eye(2, dtype=np.float32),
                permutations,
                np.full(3, -1.0),
            )
        self.assertEqual(z.shape, (2, 3))
        self.assertTrue(np.all(np.isfinite(z)))
        self.assertTrue(np.all((z >= -1.0) & (z <= 1.0)))

    def test_donor_assignment_nulls_permute_observed_z_before_normalizing(self) -> None:
        observed = np.asarray([[-1.0, 0.0, 1.0], [0.5, -0.5, 0.0]])
        permutations = np.asarray(
            [
                [[2, 1, 0], [1, 2, 0]],
                [[1, 0, 2], [2, 0, 1]],
            ],
            dtype=np.int8,
        )
        with patch.object(DEVELOP, "ETAS", (0.0, 0.5)):
            weights = DEVELOP.donor_assignment_null_weight_grid(observed, permutations)
        expected_z = np.take_along_axis(observed, permutations[0], axis=1)
        np.testing.assert_allclose(
            weights[0, 1], DEVELOP.normalized_weights(expected_z, 0.5)
        )

    def test_each_ancestry_control_selects_independently(self) -> None:
        self.assertEqual(
            DEVELOP.ANCESTRY_CONTROLS,
            (
                "cell_count",
                "total_umi",
                "occupied_batch_count",
                "pooled_effect_norm",
            ),
        )
        source = MODULE_PATH.read_text()
        self.assertIn("for control in ANCESTRY_CONTROLS:", source)
        self.assertIn("score_weight_family(", source)
        self.assertNotIn("uses_primary_selected_etas_without_retuning", source)

    def test_eta_selection_uses_family_eligibility_and_lexicographic_ties(self) -> None:
        indices = np.arange(4)
        truth = {name: np.eye(4, dtype=np.float64) for name in DEVELOP.SCORED_DATASETS}
        geometries = {
            group: {"baseline": np.eye(4, dtype=np.float64)}
            for group in DEVELOP.GROUP_SLOTS
        }
        predictions = {
            eta: {group: np.eye(4, dtype=np.float64) for group in DEVELOP.GROUP_SLOTS}
            for eta in DEVELOP.ETAS
        }
        scores = {}
        for eta in DEVELOP.ETAS:
            scores[eta] = {
                dataset: {
                    "pds": np.zeros(4),
                    "signed_credit": np.zeros(4),
                }
                for dataset in DEVELOP.SCORED_DATASETS
            }
        # Exact ties must select the smaller eta, hence identity.
        selected, _ = DEVELOP.select_eta(
            indices,
            predictions,
            scores,
            geometries,
            truth,
            inner_gallery=False,
        )
        self.assertEqual(selected, 0.0)
        # K562 is eligible after averaging its two screens, despite one negative screen.
        for metric in ("pds", "signed_credit"):
            scores[0.25]["k562_gwps"][metric][:] = -0.02
            scores[0.25]["k562_essential"][metric][:] = 0.04
            scores[0.25]["rpe1"][metric][:] = 0.01
            scores[0.25]["nadig_jurkat"][metric][:] = 0.01
        selected, _ = DEVELOP.select_eta(
            indices,
            predictions,
            scores,
            geometries,
            truth,
            inner_gallery=False,
        )
        self.assertEqual(selected, 0.25)

        control_scores = {
            eta: {
                dataset: {
                    "pds": np.zeros(4),
                    "signed_credit": np.zeros(4),
                }
                for dataset in DEVELOP.SCORED_DATASETS
            }
            for eta in DEVELOP.ETAS
        }
        for metric in ("pds", "signed_credit"):
            for dataset in DEVELOP.SCORED_DATASETS:
                control_scores[0.5][dataset][metric][:] = 0.02
        control_selected, _ = DEVELOP.select_eta(
            indices,
            predictions,
            control_scores,
            geometries,
            truth,
            inner_gallery=False,
        )
        self.assertEqual(control_selected, 0.5)

    def test_fold_family_gate_checks_all_fifteen_cells_for_both_metrics(self) -> None:
        folds = np.repeat(np.arange(5), 2)
        score = {
            dataset: {metric: np.full(10, 0.01) for metric in ("pds", "signed_credit")}
            for dataset in DEVELOP.SCORED_DATASETS
        }
        cells = DEVELOP.fold_family_cells(score, folds)
        self.assertEqual(
            sum(
                len(cells[metric][family])
                for metric in ("pds", "signed_credit")
                for family in DEVELOP.FAMILIES
            ),
            30,
        )
        self.assertTrue(
            all(
                value >= 0
                for metric in cells.values()
                for family in metric.values()
                for value in family
            )
        )

    def test_block_bootstrap_and_concentration(self) -> None:
        values = np.asarray([1.0, 1.0, -1.0, -1.0])
        blocks = np.asarray([0, 0, 1, 1])
        draws = np.asarray([[0, 1], [0, 0], [1, 1]])
        np.testing.assert_allclose(
            DEVELOP.block_bootstrap(values, blocks, draws), [0.0, 1.0, -1.0]
        )
        self.assertEqual(DEVELOP.positive_share([0.0, -1.0]), 1.0)
        self.assertAlmostEqual(DEVELOP.positive_share([1.0, 1.0]), 0.5)
        self.assertGreater(DEVELOP.positive_share([2.0, 1.0]), 0.5)
        self.assertAlmostEqual(DEVELOP.positive_share([5e-13, 5e-13]), 0.5)

    def test_failed_lock_acquisition_preserves_existing_owner_lock(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            lock = Path(directory) / "owned.lock"
            lock.write_text("other-process\n")
            with self.assertRaises(FileExistsError):
                with DEVELOP.exclusive_lock(lock):
                    self.fail("existing lock must not be acquired")
            self.assertEqual(lock.read_text(), "other-process\n")

    def test_prefreeze_paths_bind_frozen_axes_and_source_authorization(self) -> None:
        paths = set(DEVELOP.required_prefreeze_paths())
        self.assertIn(DEVELOP.DEVELOPMENT_TARGETS.resolve(), paths)
        self.assertIn(DEVELOP.RESPONSE_GENES.resolve(), paths)
        self.assertIn(DEVELOP.EXP132_AUTHORIZATION.resolve(), paths)
        self.assertIn((DEVELOP.EXP132_PLAN / "plan-manifest.json").resolve(), paths)

    def test_source_plan_is_tied_to_current_development_axes(self) -> None:
        DEVELOP.validate_source_plan_axes()

    def test_ledger_binding_checks_axis_hash_and_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "axis.txt"
            path.write_text("A\nB\n")
            record = {
                "path": str(path.relative_to("/")),
                "bytes": path.stat().st_size,
                "sha256": DEVELOP.sha256_file(path),
            }
            bindings = {path.resolve(): record["sha256"]}
            with patch.object(DEVELOP, "ROOT", Path("/")):
                DEVELOP.validate_ledger_binding(
                    record,
                    path,
                    bindings,
                    verify_file=True,
                )
                bad = dict(record, sha256="0" * 64)
                with self.assertRaises(DEVELOP.ContractError):
                    DEVELOP.validate_ledger_binding(
                        bad,
                        path,
                        bindings,
                        verify_file=True,
                    )

    def test_prefreeze_authorization_is_mandatory_and_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            missing = Path(directory) / "missing.json"
            with patch.object(DEVELOP, "DEVELOPMENT_PREFREEZE_AUTHORIZATION", missing):
                with self.assertRaises(DEVELOP.ContractError):
                    DEVELOP.validate_prefreeze_authorization()

            authorization = Path(directory) / "authorization.json"
            authorization.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "development_preoutcome_prepare_authorized": False,
                        "development_613_authorized": False,
                        "untouched_410_test_authorized": False,
                        "count_generation_authorized": False,
                        "cell_eval2_authorized": False,
                        "challenge_build_authorized": False,
                        "submission_authorized": False,
                        "bindings": {},
                    }
                )
            )
            with (
                patch.object(
                    DEVELOP, "DEVELOPMENT_PREFREEZE_AUTHORIZATION", authorization
                ),
                self.assertRaises(DEVELOP.ContractError),
            ):
                DEVELOP.validate_prefreeze_authorization()

            bound = Path(directory) / "bound.txt"
            bound.write_text("frozen")
            authorization.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "development_preoutcome_prepare_authorized": True,
                        "development_613_authorized": False,
                        "untouched_410_test_authorized": False,
                        "count_generation_authorized": False,
                        "cell_eval2_authorized": False,
                        "challenge_build_authorized": False,
                        "submission_authorized": False,
                        "bindings": {str(bound): "not-the-real-hash"},
                    }
                )
            )
            with (
                patch.object(
                    DEVELOP, "DEVELOPMENT_PREFREEZE_AUTHORIZATION", authorization
                ),
                patch.object(
                    DEVELOP, "required_prefreeze_paths", return_value=(bound,)
                ),
                self.assertRaises(DEVELOP.ContractError),
            ):
                DEVELOP.validate_prefreeze_authorization()

    def test_artifact_commands_require_authorized_explicit_provenance(self) -> None:
        authorization = {
            "bindings": {
                "prepare_provenance_utc": DEVELOP.USER_PROVENANCE_ANCHOR,
            }
        }
        with patch.dict("os.environ", {}, clear=True):
            with self.assertRaises(DEVELOP.ContractError):
                DEVELOP.authorized_provenance_utc(authorization, "prepare")
        with patch.dict(
            "os.environ",
            {"VCC_PROVENANCE_UTC": "2026-08-31T21:45:00+03:00"},
            clear=True,
        ):
            self.assertEqual(
                DEVELOP.authorized_provenance_utc(authorization, "prepare"),
                DEVELOP.USER_PROVENANCE_ANCHOR,
            )
        with patch.dict(
            "os.environ",
            {"VCC_PROVENANCE_UTC": "2026-08-31T18:46:00+00:00"},
            clear=True,
        ):
            with self.assertRaises(DEVELOP.ContractError):
                DEVELOP.authorized_provenance_utc(authorization, "prepare")

    def test_prepare_and_audit_have_no_truth_loader_or_sealed_test_path(self) -> None:
        source = MODULE_PATH.read_text()
        prepare_body = source[
            source.index("def prepare(") : source.index("def validate_preoutcome(")
        ]
        audit_body = source[source.index("def audit(") : source.index("def main(")]
        self.assertNotIn("load_development_truth", prepare_body)
        self.assertNotIn("import_experiment132", prepare_body)
        self.assertNotIn("import_experiment061", source)
        self.assertNotIn("load_development_truth", audit_body)
        self.assertNotIn("test_targets.txt", source)
        self.assertNotIn("cell-eval2", source)
        self.assertNotIn("build_submission", source)
        self.assertEqual(
            set(DEVELOP.main.__code__.co_consts)
            & {"evaluate-test", "count", "build", "submit"},
            set(),
        )

    def test_truth_file_verification_happens_after_outcome_access_lock(self) -> None:
        source = MODULE_PATH.read_text()
        develop_body = source[source.index("def develop(") : source.index("def audit(")]
        lock_position = develop_body.index("create_outcome_access_lock(authorization)")
        verification_position = develop_body.index(
            "validate_development_authorization(verify_truth_files=True)"
        )
        truth_position = develop_body.index("truth = load_development_truth(limit)")
        self.assertLess(lock_position, verification_position)
        self.assertLess(verification_position, truth_position)


if __name__ == "__main__":
    unittest.main()
