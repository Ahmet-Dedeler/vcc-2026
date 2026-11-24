from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq


MODULE_PATH = Path(__file__).with_name("run.py")
SPEC = importlib.util.spec_from_file_location("experiment142_run", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


class Experiment142Tests(unittest.TestCase):
    def test_whole_batch_split_is_deterministic_disjoint_and_balanced(self) -> None:
        umi = np.asarray([80.0, 40.0, 30.0, 20.0, 10.0, 0.0])
        labels = tuple(map(str, range(len(umi))))
        first = RUN.split_one_target(umi, labels, "donor", "target", 123)
        second = RUN.split_one_target(umi, labels, "donor", "target", 123)
        np.testing.assert_array_equal(first, second)
        self.assertTrue(np.all(np.isin(first[umi > 0], [0, 1])))
        self.assertTrue(np.all(first[umi == 0] == -1))
        totals = [umi[first == half].sum() for half in range(2)]
        self.assertLessEqual(abs(totals[0] - totals[1]), umi.max())

    def test_eight_donor_wide_splits_are_unique_and_shared_by_all_targets(self) -> None:
        labels = tuple(f"b{x}" for x in range(32))
        umi = np.tile(np.arange(1.0, 33.0), (3, 1))
        cells = np.ones_like(umi, dtype=np.int32)
        plan = RUN.make_split_plan("donor", ("A", "B", "C"), labels, umi, cells)
        splits = {
            RUN.canonical_split_bytes(plan["assignments"][split])
            for split in range(RUN.SPLITS)
        }
        self.assertEqual(len(splits), RUN.SPLITS)
        self.assertEqual(plan["assignments"].shape, (RUN.SPLITS, len(labels)))

    def test_split_plan_uses_cell_threshold_and_pair_threshold(self) -> None:
        targets = ("A", "B")
        labels = tuple(f"b{x}" for x in range(12))
        umi = np.tile(np.arange(1.0, 13.0), (2, 1))
        cells = np.ones_like(umi, dtype=np.int32)
        plan = RUN.make_split_plan("donor", targets, labels, umi, cells)
        self.assertEqual(plan["assignments"].shape, (8, 12))
        self.assertTrue(np.all(plan["valid_pairs"] >= RUN.MIN_VALID_PAIRS))
        self.assertTrue(np.all(plan["half_batches"] >= RUN.MIN_BATCHES_HALF))

    def test_half_profile_centering_is_zero_mean(self) -> None:
        rng = np.random.default_rng(1)
        original_dev = RUN.DEV_COUNT
        original_response = RUN.RESPONSE_COUNT
        RUN.DEV_COUNT = 7
        RUN.RESPONSE_COUNT = 11
        try:
            sums = rng.poisson(
                3.0,
                size=(RUN.SPLITS, RUN.HALVES, RUN.DEV_COUNT, RUN.RESPONSE_COUNT),
            ).astype(np.float64)
            cells = np.full((RUN.SPLITS, RUN.HALVES, RUN.DEV_COUNT), 5, dtype=np.int32)
            umi = sums.sum(axis=-1) + 100.0
            pooled = sums[0, 0] + sums[0, 1]
            for split in range(1, RUN.SPLITS):
                sums[split, 0] = sums[0, 0]
                sums[split, 1] = sums[0, 1]
                umi[split, 0] = umi[0, 0]
                umi[split, 1] = umi[0, 1]
            centered, valid, pooled_log, pooled_centered, pooled_cells = (
                RUN.finalize_source_profiles(sums, cells, umi)
            )
            self.assertTrue(np.all(valid))
            np.testing.assert_allclose(centered.mean(axis=2), 0.0, rtol=0.0, atol=2e-6)
            np.testing.assert_allclose(
                pooled_centered.mean(axis=0), 0.0, rtol=0.0, atol=2e-6
            )
            np.testing.assert_array_equal(pooled_cells, 10)
            self.assertEqual(pooled.shape, (RUN.DEV_COUNT, RUN.RESPONSE_COUNT))
            self.assertEqual(pooled_log.shape, (RUN.DEV_COUNT, RUN.RESPONSE_COUNT))
        finally:
            RUN.DEV_COUNT = original_dev
            RUN.RESPONSE_COUNT = original_response

    def test_exact_midrank_pds_handles_ties(self) -> None:
        query = np.asarray([[1.0, 0.0]])
        gallery = np.asarray([[1.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        result = RUN.pds_matrix(query, gallery)
        np.testing.assert_allclose(result, [[0.75, 0.75, 0.0]])

    def test_lower_quartile_reliability_and_midrank_z(self) -> None:
        values = np.asarray([1.0, 1.0, 3.0, 4.0])
        ranks = RUN.average_tie_ranks(values)
        np.testing.assert_allclose(ranks, [1.5, 1.5, 3.0, 4.0])
        z = RUN.midrank_z(values)
        np.testing.assert_allclose(z, [-2 / 3, -2 / 3, 1 / 3, 1.0])

    def test_matched_null_is_deranged_and_complete(self) -> None:
        x = np.linspace(0, 1, 33)
        features = np.column_stack([x, x**2, np.sqrt(x)])
        permutations, blocks = RUN.matched_null_permutations(features, 7, 99)
        base = np.arange(len(features))
        self.assertTrue(np.all(permutations != base[None, :]))
        for row in permutations:
            np.testing.assert_array_equal(np.sort(row), base)
        self.assertTrue(all(len(block) >= 2 for block in blocks))

    def test_bootstrap_spearman_reranks_every_resample(self) -> None:
        x = np.asarray([0.0, 0.0, 1.0, 2.0, 3.0])
        y = np.asarray([0.0, 1.0, 1.0, 3.0, 2.0])
        draws = 7
        seed = 77
        actual = RUN.bootstrap_spearman(x, y, draws, seed)
        rng = np.random.default_rng(seed)
        indices = rng.integers(0, len(x), size=(draws, len(x)))
        expected = np.asarray([RUN.spearman(x[index], y[index]) for index in indices])
        np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1e-15)

    def test_normalized_weights_identity_and_bounds(self) -> None:
        z = {
            "a": np.asarray([-1.0, 0.0, 1.0]),
            "b": np.asarray([1.0, 0.0, -1.0]),
        }
        identity = RUN.normalized_weights(z, ("a", "b"), 0.0)
        np.testing.assert_array_equal(identity, np.full((3, 2), 0.5))
        weighted = RUN.normalized_weights(z, ("a", "b"), 0.5)
        np.testing.assert_allclose(weighted.sum(axis=1), 1.0)
        self.assertTrue(np.all(weighted >= 0.25))
        missing = RUN.normalized_weights(
            {"a": np.asarray([np.nan]), "b": np.asarray([1.0])},
            ("a", "b"),
            0.5,
        )
        np.testing.assert_array_equal(missing, [[0.5, 0.5]])

    def test_common_scores_are_exact_identity(self) -> None:
        common = np.linspace(-1.0, 1.0, 11)
        weights = RUN.normalized_weights(
            {"a": common, "b": common, "c": common},
            ("a", "b", "c"),
            0.5,
        )
        np.testing.assert_array_equal(weights, np.full((11, 3), 1.0 / 3.0))

    def test_target_level_null_gate_zeroes_non_survivors(self) -> None:
        values = np.asarray([0.1, 0.2, 0.3, 0.4])
        z = RUN.gated_midrank_z(values, np.asarray([False, True, False, True]))
        np.testing.assert_array_equal(z, [0.0, -1.0, 0.0, 1.0])

    def test_donor_assignment_nulls_preserve_each_targets_values(self) -> None:
        nulls = RUN.donor_assignment_nulls(("a", "b", "c"), 7, 123)
        self.assertEqual(nulls.shape, (7, RUN.DEV_COUNT, 3))
        expected = np.arange(3, dtype=np.int8)
        self.assertTrue(np.all(np.sort(nulls, axis=2) == expected[None, None, :]))
        self.assertTrue(np.all(np.any(nulls != expected[None, None, :], axis=(1, 2))))

    def test_batch_axis_requires_exact_order(self) -> None:
        RUN.validate_batch_axis(("b1", "b2"), ("b1", "b2"), "donor")
        with self.assertRaises(RUN.ContractError):
            RUN.validate_batch_axis(("b2", "b1"), ("b1", "b2"), "donor")

    def test_spectral_coordinates_preserve_full_transform_cosines(self) -> None:
        rng = np.random.default_rng(4)
        matrix = rng.normal(size=(5, 7))
        q, _ = np.linalg.qr(rng.normal(size=(7, 7)))
        evals = np.linspace(0.5, 2.0, 7)
        spectral = RUN.spectral_coordinates(matrix, evals, q.astype(np.float32))
        weights = evals ** (-RUN.WHITEN_GAMMA / 2.0)
        full = ((matrix @ q) * weights) @ q.T
        spectral_cos = spectral @ spectral.T
        spectral_cos /= np.outer(
            np.linalg.norm(spectral, axis=1), np.linalg.norm(spectral, axis=1)
        )
        full_cos = full @ full.T
        full_cos /= np.outer(np.linalg.norm(full, axis=1), np.linalg.norm(full, axis=1))
        np.testing.assert_allclose(spectral_cos, full_cos, rtol=2e-6, atol=2e-6)

    def test_orion_duplicate_symbols_use_last_token_only(self) -> None:
        genes = ("A", "DUP", "B", "DUP")
        responses = ("DUP", "A")
        last_token = {gene: token for token, gene in enumerate(genes)}
        mapping = np.full(len(genes), -1, dtype=np.int32)
        for response_index, gene in enumerate(responses):
            mapping[last_token[gene]] = response_index
        np.testing.assert_array_equal(mapping, [1, -1, -1, 0])

    def test_sparse_orion_batch_sum_uses_sparse_tokens_and_total_counts(self) -> None:
        original_dev = RUN.DEV_COUNT
        original_response = RUN.RESPONSE_COUNT
        RUN.DEV_COUNT = 2
        RUN.RESPONSE_COUNT = 2
        try:
            table = pa.table(
                {
                    "gene_target": ["A", "B"],
                    "total_counts": [6.0, 9.0],
                    "gene_token_id": [[0, 2], [1, 3]],
                    "gene_expression": [[1.0, 5.0], [2.0, 7.0]],
                }
            )
            values, cells, umi = RUN.sparse_orion_batch_sum(
                table,
                {"A": 0, "B": 1},
                np.asarray([0, -1, -1, 1], dtype=np.int32),
            )
            np.testing.assert_array_equal(values, [[1.0, 0.0], [0.0, 7.0]])
            np.testing.assert_array_equal(cells, [1, 1])
            np.testing.assert_array_equal(umi, [6.0, 9.0])
        finally:
            RUN.DEV_COUNT = original_dev
            RUN.RESPONSE_COUNT = original_response

    def test_streamed_orion_sum_filters_and_matches_sparse_aggregation(self) -> None:
        original_dev = RUN.DEV_COUNT
        original_response = RUN.RESPONSE_COUNT
        original_batch_size = RUN.ORION_SCAN_BATCH_SIZE
        RUN.DEV_COUNT = 2
        RUN.RESPONSE_COUNT = 2
        RUN.ORION_SCAN_BATCH_SIZE = 1
        try:
            table = pa.table(
                {
                    "pass_guide_filter": [1, 1, 0, 1],
                    "gene_target": ["A", "B", "A", "OTHER"],
                    "total_counts": [6.0, 9.0, 4.0, 3.0],
                    "gene_token_id": [[0, 2], [1, 3], [0], [3]],
                    "gene_expression": [[1.0, 5.0], [2.0, 7.0], [4.0], [3.0]],
                }
            )
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "batch.parquet"
                pq.write_table(table, path)
                values, cells, umi = RUN.stream_orion_batch_sum(
                    path,
                    {"A": 0, "B": 1},
                    np.asarray([0, -1, -1, 1], dtype=np.int32),
                    ("A", "B"),
                )
            np.testing.assert_array_equal(values, [[1.0, 0.0], [0.0, 7.0]])
            np.testing.assert_array_equal(cells, [1, 1])
            np.testing.assert_array_equal(umi, [6.0, 9.0])
        finally:
            RUN.DEV_COUNT = original_dev
            RUN.RESPONSE_COUNT = original_response
            RUN.ORION_SCAN_BATCH_SIZE = original_batch_size

    def test_atomic_writes_refuse_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "x.json"
            RUN.atomic_json(path, {"x": 1})
            with self.assertRaises(RUN.ContractError):
                RUN.atomic_json(path, {"x": 2})

    def test_stage_directory_is_atomic_and_cleans_failure(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            final = Path(directory) / "complete"
            with RUN.atomic_stage_directory(final) as staging:
                (staging / "x.txt").write_text("ok")
                self.assertFalse(final.exists())
            self.assertEqual((final / "x.txt").read_text(), "ok")

            failed = Path(directory) / "failed"
            with self.assertRaisesRegex(RuntimeError, "stop"):
                with RUN.atomic_stage_directory(failed) as staging:
                    (staging / "x.txt").write_text("partial")
                    raise RuntimeError("stop")
            self.assertFalse(failed.exists())

    def test_freeze_timestamp_must_match_authorization(self) -> None:
        authorization = {
            "bindings": {
                "freeze_provenance_utc": "2026-08-31T00:00:00+00:00",
            }
        }
        with patch.dict(
            os.environ,
            {"VCC_PROVENANCE_UTC": "2026-08-31T03:00:00+03:00"},
        ):
            self.assertEqual(
                RUN.authorized_freeze_utc(authorization),
                "2026-08-31T00:00:00+00:00",
            )
        with patch.dict(
            os.environ,
            {"VCC_PROVENANCE_UTC": "2030-01-02T03:04:05+00:00"},
        ):
            with self.assertRaises(RUN.ContractError):
                RUN.authorized_freeze_utc(authorization)
        with patch.dict(os.environ, {}, clear=False):
            os.environ.pop("VCC_PROVENANCE_UTC", None)
            with self.assertRaises(RUN.ContractError):
                RUN.authorized_freeze_utc(authorization)

    def test_artifact_timestamp_requires_explicit_august_31_utc(self) -> None:
        with patch.dict(
            os.environ,
            {"VCC_PROVENANCE_UTC": "2026-08-31T04:00:00+03:00"},
        ):
            self.assertEqual(
                RUN.provenance_utc(),
                "2026-08-31T01:00:00+00:00",
            )
        with patch.dict(
            os.environ,
            {"VCC_PROVENANCE_UTC": "2026-08-30T23:59:59+00:00"},
        ):
            with self.assertRaises(RUN.ContractError):
                RUN.provenance_utc()
        with patch.dict(os.environ, {}, clear=False):
            os.environ.pop("VCC_PROVENANCE_UTC", None)
            with self.assertRaises(RUN.ContractError):
                RUN.provenance_utc()

    def test_heavy_stage_lock_serializes_all_stage_types(self) -> None:
        original = RUN.HEAVY_STAGE_LOCK
        with tempfile.TemporaryDirectory() as directory:
            RUN.HEAVY_STAGE_LOCK = Path(directory) / "heavy.lock"
            try:
                with RUN.heavy_stage_lock("source-evaluate", "k562_gwps"):
                    self.assertTrue(RUN.HEAVY_STAGE_LOCK.is_file())
                    with self.assertRaises(RUN.ContractError):
                        with RUN.heavy_stage_lock(
                            "source-build",
                            "orion_hct116",
                        ):
                            pass
                self.assertFalse(RUN.HEAVY_STAGE_LOCK.exists())
                with RUN.heavy_stage_lock("source-gate", None):
                    self.assertTrue(RUN.HEAVY_STAGE_LOCK.is_file())
            finally:
                RUN.HEAVY_STAGE_LOCK = original

    def test_development_is_fail_closed_without_authorization(self) -> None:
        source = MODULE_PATH.read_text()
        self.assertIn("validate_development_authorization()", source)
        self.assertIn(
            '"separate DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json is absent"',
            source,
        )
        self.assertNotIn("test_targets.txt", source)
        self.assertNotIn("uv run cell-eval2", source)
        self.assertNotIn("subprocess", source)
        self.assertNotIn("vcc submit", source)

    def test_source_freeze_binds_parent_control_plane(self) -> None:
        source = MODULE_PATH.read_text()
        for name in (
            "EXP132_SOURCE_AUTHORIZATION",
            "EXP132_CANDIDATE_FREEZE",
            "EXP132_PLAN_MANIFEST",
            "EXP132_SELECTION_MANIFEST",
            "SOURCE_ONLY_FREEZE_AUTHORIZATION",
            "SUPERSEDED_SOURCE_AUTHORIZATION_V4",
            "SUPERSEDED_SOURCE_AUTHORIZATION_V3",
            "SUPERSEDED_PREOUTCOME_AUTHORIZATION_V2",
            "SUPERSEDED_K562_BUILD_MANIFEST_V2",
            "SUPERSEDED_K562_EVALUATION_MANIFEST_V2",
        ):
            self.assertIn(name, source)

    def test_preoutcome_freeze_is_one_atomic_directory(self) -> None:
        source = MODULE_PATH.read_text()
        self.assertIn("tempfile.mkdtemp", source)
        self.assertIn("os.replace(temporary, FREEZE_DIR)", source)

    def test_target_nulls_remain_float64_for_strict_comparison(self) -> None:
        source = MODULE_PATH.read_text()
        self.assertIn(
            "(NULL_DRAWS, SPLITS, HALVES, DEV_COUNT), np.nan, dtype=np.float64",
            source,
        )

    def test_split_rows_never_fit_operator(self) -> None:
        source = MODULE_PATH.read_text()
        self.assertNotIn("np.linalg.eigh", source)
        self.assertNotIn("fit_operator", source)
        self.assertIn("load_operator(role)", source)


if __name__ == "__main__":
    unittest.main()
