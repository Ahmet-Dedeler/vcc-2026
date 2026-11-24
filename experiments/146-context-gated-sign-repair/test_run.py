from __future__ import annotations

import importlib.util
import inspect
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp146_run", HERE / "run.py")
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


class Experiment146Tests(unittest.TestCase):
    def test_protocol_and_single_candidate_are_frozen(self) -> None:
        protocol = RUN.protocol()
        self.assertEqual(protocol["source_schema"]["evaluation_targets"], 81)
        self.assertEqual(
            tuple(protocol["source_schema"]["sealed_410_overlap_exclusions"]),
            RUN.SEALED_410_OVERLAP_TARGETS,
        )
        self.assertEqual(protocol["target_folds"]["count"], 5)
        self.assertEqual(
            tuple(protocol["target_folds"]["expected_sizes"]), (17, 16, 16, 16, 16)
        )
        self.assertEqual(len(RUN.CANDIDATES), 1)
        self.assertEqual(RUN.FIXED_CANDIDATE.score_quantile, 0.995)
        self.assertEqual(RUN.TRAINING_ROW_CAP, 500_000)
        self.assertEqual(RUN.HGB_MAX_ITER, 100)
        self.assertEqual(RUN.SOURCE_TARGET_SUM, 50_000.0)
        self.assertEqual(protocol["response_axis"]["source_count"], 4_974)
        self.assertEqual(
            protocol["response_axis"]["source_sha256"],
            RUN.EXPECTED_SOURCE_RESPONSE_AXIS_SHA256,
        )
        self.assertEqual(
            tuple(protocol["response_axis"]["fdr_axis_exclusions"]),
            RUN.FDR_AXIS_EXCLUSIONS,
        )
        self.assertEqual(
            RUN.PREPARE_AUTHORIZATION.name,
            "SOURCE_ONLY_PREPARE_AUTHORIZATION.v4.json",
        )
        self.assertEqual(
            protocol["source_schema"]["fdr_missingness"][
                "nonfinite_entries_on_available_pairs"
            ],
            RUN.EXPECTED_FDR_NONFINITE_COUNT,
        )
        self.assertEqual(
            protocol["source_schema"]["fdr_missingness"]["all_missing_block_sha256"],
            RUN.EXPECTED_FDR_ALL_MISSING_BLOCK_SHA256,
        )

    def test_target_folds_are_deterministic_balanced_and_whole(self) -> None:
        targets = [f"T{index:02d}" for index in range(81)]
        first = RUN.deterministic_target_folds(targets)
        second = RUN.deterministic_target_folds(list(targets))
        self.assertTrue(np.array_equal(first, second))
        self.assertEqual(
            [int(np.sum(first == fold)) for fold in range(5)],
            [17, 16, 16, 16, 16],
        )
        for target in targets:
            index = targets.index(target)
            self.assertEqual(first[index], second[index])

    def test_metric_masks_are_distinct_and_exact(self) -> None:
        genes = ["A", "B", "C", "VCC", "OTHER"]
        source_targets = ["A", "B", "C"]
        pds, signed = RUN.metric_masks(genes, source_targets, ["VCC"])
        self.assertEqual(pds.tolist(), [False, False, False, False, True])
        self.assertEqual(
            signed.tolist(),
            [
                [False, True, True, False, True],
                [True, False, True, False, True],
                [True, True, False, False, True],
            ],
        )
        self.assertFalse(np.array_equal(pds, signed[0]))

    def test_metric_masks_fail_closed_on_duplicates(self) -> None:
        with self.assertRaises(RUN.ContractError):
            RUN.metric_masks(["A", "A"], ["A"], [])

    def test_multilevel_csv_index_reader_is_streamed_and_exact(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "wide.csv"
            path.write_text(
                ",context,context\ntarget,A,B\nresponse_id,,\nG1,0.1,0.2\nG2,0.3,0.4\n"
            )
            self.assertEqual(
                RUN.read_multilevel_csv_index(path),
                ("response_id", "G1", "G2"),
            )
            path.write_text(",context,context\ntarget,A,B\nG1,0.1,0.2\nG1,0.3,0.4\n")
            with self.assertRaisesRegex(RUN.ContractError, "duplicate CSV index"):
                RUN.read_multilevel_csv_index(path)

    def test_fdr_missingness_summary_keeps_missing_values_unlabelled(self) -> None:
        fdr = np.full((2, 3, 4), np.nan, dtype=np.float32)
        available = np.asarray(
            [
                [True, True, False],
                [True, False, True],
            ],
            dtype=bool,
        )
        fdr[0, 0] = [0.01, np.nan, 0.2, 0.3]
        fdr[1, 0] = [np.nan, 0.02, 0.2, 0.3]
        fdr[1, 2] = [0.01, 0.02, 0.2, 0.3]
        audit = RUN.summarize_fdr_missingness(
            fdr,
            available,
            ["C0", "C1"],
            ["T0", "T1", "T2"],
        )
        self.assertEqual(audit["available_context_target_pairs"], 4)
        self.assertEqual(audit["nonfinite_entries_on_available_pairs"], 6)
        self.assertEqual(audit["all_missing_available_blocks"], 1)
        self.assertEqual(
            audit["all_missing_blocks"],
            [{"context": "C0", "target": "T1"}],
        )
        self.assertEqual(audit["policy"], "nonfinite_is_unlabelled_never_imputed")
        fdr[0, 2, 0] = 0.1
        with self.assertRaisesRegex(
            RUN.ContractError,
            "unavailable source pairs contain FDR values",
        ):
            RUN.summarize_fdr_missingness(
                fdr,
                available,
                ["C0", "C1"],
                ["T0", "T1", "T2"],
            )

    def test_exact_nine_target_exclusion_yields_frozen_81(self) -> None:
        manifest = RUN.read_json(RUN.DEPMAP_MANIFEST)
        controls = set(manifest["dataset_contract"]["control_targets"])
        experimental = tuple(
            target
            for target in manifest["dataset_contract"]["targets"]
            if target not in controls
        )
        targets = RUN.permitted_evaluation_targets(experimental)
        self.assertEqual(len(experimental), 90)
        self.assertEqual(len(targets), 81)
        self.assertTrue(set(RUN.SEALED_410_OVERLAP_TARGETS).isdisjoint(targets))
        self.assertEqual(RUN.sha256_lines(targets), RUN.EXPECTED_EVALUATION_81_SHA256)

    def test_six_batch_by_five_target_outer_design(self) -> None:
        batches = (
            "a",
            "a",
            "b",
            "b",
            "c",
            "c",
            "d",
            "d",
            "e",
            "e",
            "f",
            "f",
            "f",
            "f",
            "f",
            "f",
        )
        batch_folds = RUN.experimental_batch_folds(batches)
        self.assertEqual(len(batch_folds), 6)
        self.assertEqual(len(batch_folds) * RUN.TARGET_FOLD_COUNT, 30)
        targets = [f"T{index:02d}" for index in range(81)]
        folds = RUN.deterministic_target_folds(targets)
        held_contexts = set(batch_folds[2].tolist())
        outer_fold = 3
        train_contexts = set(range(16)) - held_contexts
        self.assertTrue(held_contexts.isdisjoint(train_contexts))
        held_targets = set(np.flatnonzero(folds == outer_fold).tolist())
        for inner_fold in [fold for fold in range(5) if fold != outer_fold]:
            fit_targets = set(
                np.flatnonzero((folds != outer_fold) & (folds != inner_fold)).tolist()
            )
            selection_targets = set(np.flatnonzero(folds == inner_fold).tolist())
            self.assertTrue(held_targets.isdisjoint(fit_targets))
            self.assertTrue(held_targets.isdisjoint(selection_targets))
        refit_targets = set(np.flatnonzero(folds != outer_fold).tolist())
        self.assertTrue(held_targets.isdisjoint(refit_targets))
        self.assertIn(len(refit_targets), (64, 65))

    def test_complete_gallery_context_filter_removes_partial_galleries(self) -> None:
        available = np.ones((4, 5), dtype=bool)
        available[1, 3] = False
        available[3, 1] = False
        selected = RUN.complete_gallery_contexts(
            available,
            context_indices=[0, 1, 2, 3],
            target_indices=[1, 2, 3],
        )
        self.assertTrue(np.array_equal(selected, [0, 2]))
        self.assertTrue(np.all(available[np.ix_(selected, [1, 2, 3])]))

    def test_sign_repair_preserves_magnitude_norm_and_energy(self) -> None:
        base = np.asarray([10.0, 2.0, 1.0, 0.5])
        probability = np.asarray([0.99, 0.98, 0.97, 0.96])
        eligible = np.ones(4, dtype=bool)
        candidate = RUN.Candidate(
            score_quantile=0.0,
            energy_cap=0.05,
            pds_axis_energy_cap=0.05,
            selected_count_cap=2,
            signed_l1_mass_cap=0.5,
        )
        repaired, selected, audit = RUN.apply_sign_repair(
            base,
            probability,
            eligible,
            candidate,
            pds_mask=np.ones(4, dtype=bool),
            threshold=0.0,
        )
        self.assertTrue(np.array_equal(np.abs(repaired), np.abs(base)))
        self.assertAlmostEqual(float(repaired @ repaired), float(base @ base))
        self.assertLessEqual(audit.energy_fraction, candidate.energy_cap + 1e-12)
        self.assertLessEqual(
            audit.pds_axis_energy_fraction,
            candidate.pds_axis_energy_cap + 1e-12,
        )
        self.assertLessEqual(audit.selected_count, candidate.selected_count_cap)
        self.assertLessEqual(
            audit.signed_l1_mass_fraction,
            candidate.signed_l1_mass_cap + 1e-12,
        )
        self.assertFalse(selected[0])
        self.assertTrue(np.any(selected[1:]))

    def test_identity_fallback_is_bit_exact(self) -> None:
        base = np.asarray([1.0, -2.0, 3.0])
        repaired, selected, audit = RUN.apply_sign_repair(
            base,
            np.ones(3),
            np.ones(3, dtype=bool),
            None,
        )
        self.assertTrue(np.array_equal(repaired, base))
        self.assertFalse(np.any(selected))
        self.assertEqual(audit.energy_fraction, 0.0)
        self.assertEqual(audit.angle_degrees, 0.0)

    def test_zero_tail_uses_the_exact_flip_selected_mask(self) -> None:
        arrays = {
            "baseline": np.asarray([[4.0, 2.0, 1.0, 0.5]], dtype=np.float32),
            "signed_masks": np.ones((1, 4), dtype=bool),
            "pds_mask": np.ones(4, dtype=bool),
            "donor_available": np.ones((3, 1), dtype=bool),
        }
        probabilities = np.asarray([[0.99, 0.98, 0.97, 0.96]], dtype=np.float32)
        candidate = RUN.Candidate(
            score_quantile=0.0,
            energy_cap=0.25,
            pds_axis_energy_cap=0.25,
            selected_count_cap=2,
            signed_l1_mass_cap=0.5,
        )
        flipped, flip_mask, _ = RUN.repair_probability_matrix(
            arrays=arrays,
            target_indices=[0],
            probabilities=probabilities,
            candidate=candidate,
            action="flip",
        )
        zeroed, zero_mask, _ = RUN.repair_probability_matrix(
            arrays=arrays,
            target_indices=[0],
            probabilities=probabilities,
            candidate=candidate,
            action="zero",
        )
        self.assertTrue(np.array_equal(flip_mask, zero_mask))
        self.assertTrue(np.any(zero_mask))
        self.assertTrue(np.all(zeroed[zero_mask] == 0.0))
        self.assertTrue(
            np.array_equal(zeroed[~zero_mask], arrays["baseline"][~zero_mask])
        )
        self.assertTrue(
            np.array_equal(flipped[flip_mask], -arrays["baseline"][flip_mask])
        )

    def test_signed_credit_uses_row_specific_mask(self) -> None:
        prediction = np.asarray([[2.0, 5.0, -1.0], [2.0, 5.0, -1.0]])
        truth = np.asarray([[1.0, -1.0, -1.0], [1.0, -1.0, -1.0]])
        masks = np.asarray([[False, True, True], [True, False, True]])
        observed = RUN.signed_credit_per_target(prediction, truth, masks)
        self.assertTrue(np.allclose(observed, [-4.0 / 6.0, 3.0 / 3.0]))

    def test_pds_mask_changes_retrieval_geometry(self) -> None:
        prediction = np.asarray([[1.0, 0.0, 1.0], [0.0, 1.0, 1.0]])
        truth = prediction.copy()
        score = RUN.pds_per_target(
            prediction,
            truth,
            np.asarray([True, True, False]),
        )
        self.assertTrue(np.array_equal(score, np.ones(2)))

    def test_fixed_candidate_selection_falls_back_to_identity(self) -> None:
        records = {
            arm: [
                {
                    "batch": batch,
                    "pds_delta": 0.0,
                    "signed_credit_delta": -0.001,
                }
                for batch in ("a", "b")
            ]
            for arm in RUN.ALL_ARMS
        }
        selected, audit = RUN.select_fixed_primary(records)
        self.assertIsNone(selected)
        self.assertEqual(audit["selected"], "identity")

    def test_feature_block_and_fixed_estimators_run_on_synthetic_data(self) -> None:
        rng = np.random.default_rng(146)
        gene_count = 40
        donor_count = 3
        donor_effects = rng.normal(size=(donor_count, 2, gene_count)).astype(np.float32)
        donor_controls = rng.normal(size=(donor_count, gene_count)).astype(np.float32)
        features = RUN.build_feature_block(
            context_control=rng.normal(size=gene_count),
            context_score=rng.normal(size=RUN.CONTROL_PC_RANK),
            donor_controls=donor_controls,
            donor_scores=rng.normal(size=(donor_count, RUN.CONTROL_PC_RANK)),
            donor_effects=donor_effects,
            donor_available=np.ones((donor_count, 2), dtype=bool),
            baseline=np.mean(donor_effects, axis=0),
            target_index=0,
        )
        self.assertEqual(features.shape, (gene_count, len(RUN.FEATURE_NAMES)))
        self.assertTrue(np.all(np.isfinite(features)))

        x = np.vstack([features, features + 0.05])
        y = np.tile(np.asarray([0, 1], dtype=np.uint8), gene_count)
        weight = np.ones(len(y))
        models = RUN.fit_estimators(x, y, weight, 123)
        self.assertEqual(set(models), set(RUN.MODEL_ARMS))
        self.assertEqual(models["full_hgb"]["estimator"].max_iter, 100)
        self.assertEqual(models["donor_hgb"]["estimator"].max_iter, 100)
        self.assertTrue(
            np.array_equal(
                models["donor_hgb"]["columns"],
                RUN.DONOR_ONLY_FEATURE_INDICES,
            )
        )
        for name, model in models.items():
            probability = RUN.predict_wrong_probability(model, name, features)
            self.assertEqual(probability.shape, (gene_count,))
            self.assertTrue(np.all((probability >= 0) & (probability <= 1)))

    def test_availability_stratified_nulls_preserve_blocks(self) -> None:
        targets = np.arange(6)
        available = np.asarray(
            [
                [True, True, True, True, True, True],
                [True, True, False, False, False, True],
                [True, True, True, True, True, False],
                [True, True, True, True, True, True],
            ],
            dtype=bool,
        )
        mappings, audit = RUN.availability_stratified_target_permutations(
            targets,
            available,
            np.arange(4),
            draws=7,
            seed=9,
        )
        self.assertEqual(len(mappings), 7)
        self.assertTrue(audit["availability_pattern_stratified"])
        self.assertTrue(audit["all_training_blocks_preserved"])
        self.assertEqual(audit["singleton_target_indices"], [5])
        for mapping in mappings:
            self.assertEqual(set(mapping), set(targets.tolist()))
            self.assertEqual(set(mapping.values()), set(targets.tolist()))
            self.assertEqual(mapping[5], 5)
            self.assertTrue(all(mapping[target] != target for target in range(5)))
            for source, target in mapping.items():
                self.assertTrue(
                    np.array_equal(available[:, source], available[:, target])
                )

    def test_permutation_validator_rejects_changed_availability(self) -> None:
        available = np.asarray(
            [
                [True, True, True],
                [True, False, True],
            ],
            dtype=bool,
        )
        RUN.validate_target_truth_permutation(
            {0: 2, 1: 1, 2: 0},
            [0, 1, 2],
            available,
            [0, 1],
        )
        with self.assertRaisesRegex(
            RUN.ContractError,
            "availability pattern",
        ):
            RUN.validate_target_truth_permutation(
                {0: 1, 1: 0, 2: 2},
                [0, 1, 2],
                available,
                [0, 1],
            )

    def test_real_orion_duplicate_axes_use_last_occurrence(self) -> None:
        for name in ("orion_hct116", "orion_hek293t"):
            genes = RUN.read_lines(
                RUN.DONOR_ROOTS[name] / "genes.txt",
                unique=False,
            )
            index = RUN.last_occurrence_index(genes)
            self.assertLess(len(index), len(genes))
            self.assertEqual(
                index["MKKS"],
                max(position for position, gene in enumerate(genes) if gene == "MKKS"),
            )

    def test_exact_five_single_donor_targets_are_allowed_for_features(self) -> None:
        manifest = RUN.read_json(RUN.DEPMAP_MANIFEST)
        controls = set(manifest["dataset_contract"]["control_targets"])
        experimental = tuple(
            target
            for target in manifest["dataset_contract"]["targets"]
            if target not in controls
        )
        targets = RUN.permitted_evaluation_targets(experimental)
        donor_target_sets = {
            name: set(RUN.read_lines(root / "perts.txt"))
            for name, root in RUN.DONOR_ROOTS.items()
        }
        donor_available = np.asarray(
            [
                [target in donor_target_sets[name] for target in targets]
                for name in RUN.DONOR_ROOTS
            ],
            dtype=bool,
        )
        audit = RUN.validate_donor_support(donor_available, targets)
        self.assertEqual(
            tuple(audit["single_donor_targets"]),
            RUN.SINGLE_DONOR_TARGETS,
        )
        self.assertEqual(audit["support_histogram"], {"1": 5, "2": 4, "3": 72})

        rng = np.random.default_rng(146)
        donor_effects = rng.normal(size=(3, 1, 12)).astype(np.float32)
        donor_controls = rng.normal(size=(3, 12)).astype(np.float32)
        one_donor = np.asarray([[True], [False], [False]])
        features = RUN.build_feature_block(
            context_control=rng.normal(size=12),
            context_score=rng.normal(size=RUN.CONTROL_PC_RANK),
            donor_controls=donor_controls,
            donor_scores=rng.normal(size=(3, RUN.CONTROL_PC_RANK)),
            donor_effects=donor_effects,
            donor_available=one_donor,
            baseline=donor_effects[0],
            target_index=0,
        )
        self.assertTrue(np.all(np.isfinite(features)))
        self.assertTrue(np.allclose(features[:, 8], 1.0 / 3.0))
        self.assertTrue(np.all(features[:, 19] == 0.0))
        agreement = RUN.donor_agreement_probability(
            {
                "baseline": donor_effects[0],
                "donor_available": one_donor,
                "donor_effects": donor_effects,
            },
            [0],
        )
        self.assertTrue(np.array_equal(agreement, np.zeros_like(agreement)))
        repaired, selected, _ = RUN.repair_probability_matrix(
            arrays={
                "baseline": donor_effects[0],
                "donor_available": one_donor,
                "signed_masks": np.ones((1, 12), dtype=bool),
                "pds_mask": np.ones(12, dtype=bool),
            },
            target_indices=[0],
            probabilities=np.ones((1, 12), dtype=np.float32),
            candidate=RUN.FIXED_CANDIDATE,
        )
        self.assertFalse(np.any(selected))
        self.assertTrue(np.array_equal(repaired, donor_effects[0]))

    def test_source_reconstruction_uses_target_sum_50000(self) -> None:
        profile = RUN.reconstruct_log1p_profile(np.asarray([1.0, 2.0, 7.0]))
        self.assertAlmostEqual(float(np.expm1(profile).sum()), 50_000.0, places=2)

    def test_complete_delta_hashing_is_forbidden(self) -> None:
        with self.assertRaisesRegex(RUN.ContractError, "refuse to hash"):
            RUN.sha256_file(Path("/tmp/delta.npy"))
        bindings = RUN.prepare_binding_paths()
        self.assertFalse(any(path.name == "delta.npy" for path in bindings))
        self.assertNotIn(RUN.DEPMAP_PSEUDOBULK, bindings)
        self.assertNotIn(RUN.DEPMAP_FDR, bindings)

    def test_selected_slice_requires_lock_and_reads_only_requested_cells(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "delta.npy"
            matrix = np.arange(80, dtype=np.float32).reshape(10, 8)
            np.save(path, matrix)
            ledger = {
                "files": {
                    "synthetic/delta.npy": {
                        "bytes": path.stat().st_size,
                        "sha256": "a" * 64,
                    }
                }
            }
            with mock.patch.object(
                RUN,
                "repo_relative",
                return_value="synthetic/delta.npy",
            ):
                with self.assertRaisesRegex(RUN.ContractError, "outside source-access"):
                    RUN.selected_matrix_slice(path, [2, 7], [1, 5], ledger=ledger)
                with RUN.source_array_access():
                    selected = RUN.selected_matrix_slice(
                        path,
                        [2, 7],
                        [1, 5],
                        ledger=ledger,
                    )
            self.assertTrue(np.array_equal(selected, matrix[np.ix_([2, 7], [1, 5])]))
            self.assertFalse(isinstance(selected, np.memmap))

    def test_depmap_outcomes_require_source_access_lock(self) -> None:
        with self.assertRaisesRegex(RUN.ContractError, "outside source-access lock"):
            RUN.load_depmap_source()

    def test_historical_whitener_excludes_evaluation_target_coordinates(self) -> None:
        donor_effects = np.asarray(
            [
                [[4.0, 3.0, 1.0, 0.0], [2.0, 1.0, 0.0, 1.0]],
                [[3.0, 4.0, 1.0, 1.0], [1.0, 2.0, 1.0, 0.0]],
            ],
            dtype=np.float32,
        )
        available = np.ones((2, 2), dtype=bool)
        whitener, audit = RUN.fit_historical_panel_whitener(
            donor_effects,
            available,
            ["T0", "T1", "G2", "G3"],
            ["T0", "T1"],
            np.asarray([False, False, True, True]),
        )
        transformed = whitener.apply(np.asarray([[8.0, 9.0, 1.0, 2.0]]))
        self.assertTrue(np.array_equal(transformed[:, :2], [[8.0, 9.0]]))
        self.assertTrue(np.all(np.isfinite(transformed)))
        self.assertTrue(audit["evaluation_target_coordinates_excluded"])

    def test_score_uses_complete_supplied_gallery(self) -> None:
        arrays = {
            "available": np.ones((1, 4), dtype=bool),
            "baseline": np.asarray(
                [[0.7, 0.7], [0.0, 1.0], [0.8, 0.6], [0.9, 0.4]],
                dtype=np.float32,
            ),
            "effects": np.asarray(
                [
                    [
                        [1.0, 0.0],
                        [0.0, 1.0],
                        [0.8, 0.6],
                        [0.9, 0.4],
                    ]
                ],
                dtype=np.float32,
            ),
            "fdr": np.zeros((1, 4, 2), dtype=np.float32),
            "signed_masks": np.ones((4, 2), dtype=bool),
            "pds_mask": np.ones(2, dtype=bool),
            "donor_available": np.ones((3, 4), dtype=bool),
        }
        scored = RUN.score_candidate_context(
            arrays=arrays,
            context_index=0,
            target_indices=[0, 1],
            gallery_indices=[0, 1, 2, 3],
            probabilities=np.zeros((2, 2), dtype=np.float32),
            candidate=None,
        )
        self.assertEqual(scored["gallery_size"], 4)
        self.assertEqual(len(scored["pds_delta"]), 2)

    def test_inner_selection_uses_complete_gallery_not_mini_galleries(self) -> None:
        source = inspect.getsource(RUN.inner_complete_gallery_evaluation)
        self.assertIn("gallery_indices=inner_targets", source)
        self.assertNotIn("gallery_indices=validation_targets", source)
        self.assertIn("len(inner_targets) in (64, 65)", source)

    def test_bootstraps_and_fdr_sensitivity_are_frozen(self) -> None:
        target_values = np.linspace(-0.1, 0.2, 81)
        target_folds = RUN.deterministic_target_folds(
            [f"T{index:02d}" for index in range(81)]
        )
        first = RUN.target_block_bootstrap(
            target_values,
            target_folds,
            draws=20,
            seed=7,
        )
        second = RUN.target_block_bootstrap(
            target_values,
            target_folds,
            draws=20,
            seed=7,
        )
        self.assertTrue(np.array_equal(first, second))
        self.assertEqual(RUN.FDR_SENSITIVITY, (0.01, 0.05, 0.10))
        self.assertEqual(RUN.PERMUTATION_NULL_DRAWS, 9)

    def test_target_metric_uses_equal_batch_weighting(self) -> None:
        records = [
            {
                "batch": "a",
                "available_targets": np.asarray([0]),
                "metrics": {
                    "full_hgb": {
                        "signed_credit_delta": np.asarray([1.0]),
                    }
                },
            },
            {
                "batch": "a",
                "available_targets": np.asarray([0]),
                "metrics": {
                    "full_hgb": {
                        "signed_credit_delta": np.asarray([1.0]),
                    }
                },
            },
            {
                "batch": "b",
                "available_targets": np.asarray([0]),
                "metrics": {
                    "full_hgb": {
                        "signed_credit_delta": np.asarray([5.0]),
                    }
                },
            },
        ]
        observed = RUN.equal_batch_target_metric(
            records,
            arm="full_hgb",
            metric="signed_credit_delta",
            target_count=1,
        )
        self.assertTrue(np.allclose(observed, [3.0]))

    def test_resource_report_contract(self) -> None:
        report = {
            "status": "complete",
            "max_rss_gib": 18.0,
            "min_free_percent": 25,
            "peak_rss_gib": 4.0,
            "minimum_free_percent": 61,
            "samples": 3,
        }
        RUN.validate_memory_report(report)
        bad = {**report, "minimum_free_percent": 24}
        with self.assertRaises(RUN.ContractError):
            RUN.validate_memory_report(bad)

    def test_authorization_rejects_forbidden_permission(self) -> None:
        bindings = {"x": "abc"}
        document = {
            "schema_version": 1,
            "experiment": RUN.HERE.name,
            "stage": "source-only-prepare",
            "status": "authorized",
            "confirmation": RUN.PREPARE_CONFIRMATION,
            "bindings": bindings,
            "permissions": {
                "source_prepare": True,
                "source_evaluate": False,
                "development_613": True,
                "sealed_410": False,
                "challenge_control_expression": False,
                "cell_eval2": False,
                "count_generation": False,
                "challenge_build": False,
                "submission": False,
                "writes_outside_experiment_146": False,
            },
        }
        with self.assertRaises(RUN.ContractError):
            RUN.validate_authorization_document(
                document,
                stage="source-only-prepare",
                confirmation=RUN.PREPARE_CONFIRMATION,
                bindings=bindings,
            )

    def test_evaluation_authorization_binds_memory_guard(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            authorization = root / "authorization.json"
            prepared_manifest = root / "artifact-manifest.json"
            prepared_manifest.write_text("{}\n")
            captured: list[Path] = []

            def fake_bindings(paths: object) -> dict[str, str]:
                captured.extend(list(paths))
                return {"frozen": "hash"}

            authorization.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "experiment": RUN.HERE.name,
                        "stage": "source-only-evaluate",
                        "status": "authorized",
                        "confirmation": RUN.EVALUATE_CONFIRMATION,
                        "bindings": {"frozen": "hash"},
                        "permissions": {
                            "source_prepare": False,
                            "source_evaluate": True,
                            "development_613": False,
                            "sealed_410": False,
                            "challenge_control_expression": False,
                            "cell_eval2": False,
                            "count_generation": False,
                            "challenge_build": False,
                            "submission": False,
                            "authorization_creation": False,
                            "writes_outside_experiment_146": False,
                        },
                    }
                )
            )
            with (
                mock.patch.object(RUN, "EVALUATE_AUTHORIZATION", authorization),
                mock.patch.object(RUN, "expected_bindings", side_effect=fake_bindings),
            ):
                RUN.validate_stage_authorization(
                    "source-only-evaluate",
                    RUN.EVALUATE_CONFIRMATION,
                    prepared_manifest,
                )
            self.assertIn(RUN.MEMORY_GUARD, captured)

    def test_evaluation_loads_prepared_arrays_inside_memory_guard(self) -> None:
        source = inspect.getsource(RUN.evaluate_source)
        self.assertLess(
            source.index('with process_tree_memory_guard("evaluate")'),
            source.index("plan, arrays = validate_prepared()"),
        )

    def test_prepared_validation_reasserts_response_axis_and_fdr_missingness(
        self,
    ) -> None:
        source = inspect.getsource(RUN.validate_prepared)
        self.assertIn("EXPECTED_SOURCE_RESPONSE_GENE_COUNT", source)
        self.assertIn("EXPECTED_SOURCE_RESPONSE_AXIS_SHA256", source)
        self.assertIn("FDR_AXIS_EXCLUSIONS", source)
        self.assertIn("summarize_fdr_missingness", source)
        self.assertIn("validate_fdr_missingness", source)

    def test_missing_prepare_authorization_fails_before_hashing_sources(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            missing = Path(temporary) / "missing-authorization.json"
            with (
                mock.patch.object(RUN, "PREPARE_AUTHORIZATION", missing),
                mock.patch.object(
                    RUN,
                    "expected_bindings",
                    side_effect=AssertionError("source hashing must not begin"),
                ),
            ):
                with self.assertRaisesRegex(RUN.ContractError, "missing authorization"):
                    RUN.validate_stage_authorization(
                        "source-only-prepare",
                        RUN.PREPARE_CONFIRMATION,
                    )

    def test_manifest_detects_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "value.txt").write_text("before\n")
            manifest = RUN.directory_manifest(directory)
            (directory / "artifact-manifest.json").write_text(json.dumps(manifest))
            RUN.validate_directory_manifest(directory)
            (directory / "value.txt").write_text("after\n")
            with self.assertRaises(RUN.ContractError):
                RUN.validate_directory_manifest(directory)

    def test_cli_has_no_development_or_submission_entrypoint(self) -> None:
        text = (HERE / "run.py").read_text()
        self.assertIn('subparsers.add_parser("audit"', text)
        self.assertIn('"prepare", help=', text)
        self.assertIn('"evaluate", help=', text)
        self.assertNotIn('add_parser("develop"', text)
        self.assertNotIn('add_parser("submit"', text)
        self.assertNotIn("import cell_eval2", text)
        self.assertNotIn("write_submission(", text)
        self.assertNotIn("emit_counts(", text)
        self.assertNotIn("context_A.h5ad", text)
        self.assertNotIn("context_B.h5ad", text)
        self.assertNotIn("context_C.h5ad", text)

    def test_required_controls_nulls_and_resource_policy_are_present(self) -> None:
        protocol = RUN.protocol()
        controls = set(protocol["estimator"]["matched_controls"])
        self.assertEqual(len(controls), 5)
        self.assertTrue(any("zero-tail" in value for value in controls))
        self.assertTrue(any("donor-agreement" in value for value in controls))
        self.assertTrue(any("full feature" in value for value in controls))
        self.assertTrue(any("donor-only features" in value for value in controls))
        self.assertTrue(any("donor-only HistGradient" in value for value in controls))
        self.assertEqual(
            protocol["source_metrics"]["whole_target_refit_permutation"]["draws"],
            RUN.PERMUTATION_NULL_DRAWS,
        )
        self.assertTrue(
            protocol["source_metrics"]["whole_target_refit_permutation"][
                "availability_pattern_stratified"
            ]
        )
        self.assertEqual(
            protocol["source_metrics"]["whole_target_refit_permutation"][
                "singleton_class_policy"
            ],
            "fixed_and_recorded",
        )
        self.assertEqual(
            tuple(protocol["baseline"]["single_donor_targets"]),
            RUN.SINGLE_DONOR_TARGETS,
        )
        self.assertIn("target_block_bootstrap", protocol["source_metrics"])
        self.assertEqual(
            tuple(protocol["source_metrics"]["fdr_sensitivity"]),
            RUN.FDR_SENSITIVITY,
        )
        self.assertEqual(
            protocol["resource_policy"]["max_process_tree_rss_gib"],
            18.0,
        )
        self.assertEqual(
            protocol["resource_policy"]["min_system_free_percent"],
            25,
        )


if __name__ == "__main__":
    unittest.main()
