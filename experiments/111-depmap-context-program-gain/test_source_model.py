from __future__ import annotations

import importlib.util
import json
import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "exp111_source_model", HERE / "source_model.py"
)
assert SPEC is not None and SPEC.loader is not None
MODEL = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODEL
SPEC.loader.exec_module(MODEL)


class SourceModelTests(unittest.TestCase):
    def test_hashed_target_split_is_fixed_and_complete(self) -> None:
        manifest = json.loads((HERE / "SOURCE_MANIFEST.json").read_text())
        controls = set(manifest["dataset_contract"]["control_targets"])
        targets = [
            target
            for target in manifest["dataset_contract"]["targets"]
            if target not in controls
        ]
        split = MODEL.hashed_target_split(targets)
        self.assertEqual(
            {key: len(value) for key, value in split.items()},
            {"fit": 54, "tune": 18, "test": 18},
        )
        self.assertTrue(set(split["fit"]).isdisjoint(split["tune"]))
        self.assertTrue(set(split["fit"]).isdisjoint(split["test"]))
        self.assertTrue(set(split["tune"]).isdisjoint(split["test"]))
        self.assertEqual(
            sorted(split["fit"] + split["tune"] + split["test"]),
            list(range(90)),
        )

    def test_context_weights_give_each_batch_equal_mass(self) -> None:
        batches = ["a", "a", "a", "b", "c", "c"]
        indices = list(range(len(batches)))
        weights = MODEL.context_weights(indices, batches)
        self.assertTrue(np.isclose(weights[:3].sum(), 1 / 3))
        self.assertTrue(np.isclose(weights[3], 1 / 3))
        self.assertTrue(np.isclose(weights[4:].sum(), 1 / 3))
        self.assertTrue(np.isclose(weights.sum(), 1.0))

    def test_pds_uses_exact_midranks_and_zero_rows(self) -> None:
        similarity = np.asarray(
            [
                [1.0, 0.0, 0.0],
                [0.2, 0.2, 0.1],
                [0.0, 0.0, 0.0],
            ]
        )
        observed = MODEL.pds_per_target(similarity)
        self.assertTrue(np.allclose(observed, [1.0, 0.75, 0.5]))

    def test_angular_cap_is_respected(self) -> None:
        q = np.asarray([[1.0, 0.0]])
        r = np.asarray([[0.0, 1.0]])
        scale, angle = MODEL.capped_residual_scale(
            np.sum(q * q, axis=1),
            np.sum(q * r, axis=1),
            np.sum(r * r, axis=1),
            5.0,
        )
        self.assertGreater(scale[0], 0.0)
        self.assertLess(scale[0], 1.0)
        self.assertLessEqual(angle[0], 5.0 + 1e-10)
        self.assertGreaterEqual(angle[0], 5.0 - 1e-8)

    def test_zero_beta_is_exact_identity_in_pds(self) -> None:
        q = np.eye(3)
        truth = np.eye(3)
        scores, angles = MODEL.candidate_context_scores(
            q_projection=np.zeros((3, 2)),
            q_norm2=np.ones(3),
            base_truth_dot=q @ truth.T,
            truth_projection=np.zeros((3, 2)),
            truth_norm=np.ones(3),
            predicted_gain=np.ones(2),
            beta=0.0,
            cap_degrees=1.0,
        )
        self.assertTrue(np.array_equal(scores, np.ones(3)))
        self.assertTrue(np.array_equal(angles, np.zeros(3)))

    def test_derangements_have_no_fixed_points(self) -> None:
        rng = np.random.default_rng(123)
        for size in (4, 8, 16, 18, 54):
            permutation = MODEL.random_derangement(size, rng)
            self.assertEqual(sorted(permutation.tolist()), list(range(size)))
            self.assertTrue(np.all(permutation != np.arange(size)))

    def test_candidate_grid_matches_preregistered_cartesian_product(self) -> None:
        protocol = json.loads((HERE / "preregistered_protocol.json").read_text())
        candidates = MODEL.candidate_grid(protocol)
        self.assertEqual(len(candidates), 144)
        self.assertEqual(len({candidate.key for candidate in candidates}), 144)

    def test_selection_tie_break_prefers_exact_identity_simplest_arm(self) -> None:
        candidates = [
            MODEL.Candidate(8, 4, 1.0, 0.0, 5.0),
            MODEL.Candidate(4, 2, 100.0, 0.0, 1.0),
            MODEL.Candidate(4, 2, 100.0, 0.25, 1.0),
        ]
        scores = {candidate.key: 0.0 for candidate in candidates}
        selected = MODEL.select_candidate(candidates, scores)
        self.assertEqual(selected, MODEL.Candidate(4, 2, 100.0, 0.0, 1.0))

    def test_synthetic_nested_pipeline_runs_all_six_outer_batches(self) -> None:
        rng = np.random.default_rng(7)
        contexts = [f"C{index}" for index in range(6)]
        batches = [f"B{index}" for index in range(6)]
        targets = [f"T{index:02d}" for index in range(90)]
        effects = rng.normal(size=(6, 90, 20)).astype(np.float32)
        controls = rng.normal(size=(6, 20)).astype(np.float32)
        flat = effects.reshape(6 * 90, 20)
        control_mean, control_components, control_scale, control_scores = (
            MODEL.fit_control_pca(controls, batches)
        )
        data = MODEL.SourceData(
            genes=[f"G{index}" for index in range(20)],
            contexts=contexts,
            batches=batches,
            batch_order=batches,
            targets=targets,
            controls=[],
            effects=effects,
            available=np.ones((6, 90), dtype=bool),
            control_profiles=controls,
            gram=np.asarray(flat @ flat.T, dtype=np.float64),
            control_mean=control_mean,
            control_components=control_components,
            control_score_scale=control_scale,
            control_scores=control_scores,
        )
        split = MODEL.hashed_target_split(targets)
        candidates = [
            MODEL.Candidate(4, 2, 100.0, 0.0, 1.0),
            MODEL.Candidate(4, 2, 10.0, 0.25, 3.0),
        ]
        result = MODEL.run_nested_selection(
            data,
            split,
            candidates,
            {},
            MODEL.identity_target_mapping(data),
            MODEL.identity_context_mapping(data),
            detailed=False,
        )
        self.assertEqual(len(result["outer_folds"]), 6)
        self.assertTrue(np.isfinite(result["macro_raw_pds_delta"]))


if __name__ == "__main__":
    unittest.main()
