#!/usr/bin/env python3
"""Synthetic and static tests for Experiment 144 development evaluation."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RUN_PATH = HERE / "evaluate.py"


def load_run():
    spec = importlib.util.spec_from_file_location(
        "experiment144_development_tested",
        RUN_PATH,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


M = load_run()


class Experiment144DevelopmentTests(unittest.TestCase):
    def test_protocol_freezes_one_nonidentity_candidate(self) -> None:
        protocol = M.protocol()
        self.assertEqual(protocol["candidate"]["name"], M.CANDIDATE_NAME)
        self.assertFalse(protocol["candidate"]["new_grid_or_selection"])
        self.assertEqual(protocol["candidate"]["blend_eta"], 0.25)
        self.assertEqual(protocol["target_alignment_null"]["draws"], 99)

    def test_command_surface_is_exact(self) -> None:
        self.assertEqual(M.COMMANDS, ("prepare", "evaluate", "audit"))
        parser = M.build_parser()
        for invalid in ("develop", "score", "build", "submit"):
            with self.subTest(invalid=invalid):
                with self.assertRaises(SystemExit):
                    parser.parse_args([invalid])

    def test_prepare_paths_do_not_include_truth_containers(self) -> None:
        normalized = [str(path.resolve()) for path in M.prepare_control_paths()]
        forbidden = (
            "data/donors/k562_essential/delta.npy",
            "data/donors/rpe1/delta.npy",
            "data/donors/nadig_jurkat/delta.npy",
        )
        self.assertFalse(
            any(path.endswith(value) for path in normalized for value in forbidden)
        )

    def test_old_whitener_inverse_roundtrip(self) -> None:
        rng = np.random.default_rng(144)
        matrix = rng.normal(size=(9, 12)).astype(np.float32)
        q, _ = np.linalg.qr(rng.normal(size=(12, 12)))
        evecs = q.astype(np.float32)
        evals = np.geomspace(0.1, 3.0, 12)
        whitened = M.CORE.apply_operator(matrix, evals, evecs, M.OLD_GAMMA)
        recovered = M.CORE.apply_operator(
            whitened,
            evals,
            evecs,
            -M.OLD_GAMMA,
        )
        np.testing.assert_allclose(recovered, matrix, rtol=2e-6, atol=2e-6)

    def test_target_nulls_are_fold_preserving_derangements(self) -> None:
        folds = np.repeat(np.arange(M.FOLD_COUNT), (123, 123, 123, 122, 122))
        permutations = M.make_null_permutations(folds)
        indices = np.arange(M.DEV_COUNT)
        self.assertEqual(permutations.shape, (M.NULL_DRAWS, M.DEV_COUNT))
        self.assertTrue(np.all(permutations != indices[None, :]))
        self.assertTrue(np.all(folds[permutations] == folds[None, :]))
        self.assertTrue(np.array_equal(permutations, M.make_null_permutations(folds)))

    def test_pds_rewards_identity_matching(self) -> None:
        similarity = np.eye(M.DEV_COUNT)
        pds = M.pds_from_square_similarity(similarity)
        np.testing.assert_array_equal(pds, np.ones(M.DEV_COUNT))

    def test_family_macro_weights_k562_once(self) -> None:
        score = {
            "k562_gwps": {
                "pds": np.ones(M.DEV_COUNT),
                "signed_credit": np.ones(M.DEV_COUNT),
            },
            "k562_essential": {
                "pds": np.zeros(M.DEV_COUNT),
                "signed_credit": np.zeros(M.DEV_COUNT),
            },
            "rpe1": {
                "pds": np.full(M.DEV_COUNT, 2.0),
                "signed_credit": np.full(M.DEV_COUNT, 2.0),
            },
            "nadig_jurkat": {
                "pds": np.full(M.DEV_COUNT, 3.0),
                "signed_credit": np.full(M.DEV_COUNT, 3.0),
            },
        }
        summary = M.summarize(score)
        self.assertEqual(summary["families"]["pds"]["k562"], 0.5)
        self.assertEqual(summary["macro"]["pds"], (0.5 + 2.0 + 3.0) / 3.0)

    def test_bootstrap_constant_signal_is_constant(self) -> None:
        blocks = np.arange(M.DEV_COUNT, dtype=np.int64) % M.BLOCK_COUNT
        draws = np.tile(np.arange(M.BLOCK_COUNT), (7, 1))
        score = {
            dataset: {
                "pds": np.full(M.DEV_COUNT, 0.01),
                "signed_credit": np.full(M.DEV_COUNT, 0.02),
            }
            for dataset in M.SCORED_DATASETS
        }
        summary, arrays = M.bootstrap_summary(score, blocks, draws)
        np.testing.assert_allclose(arrays["pds"], 0.01)
        np.testing.assert_allclose(arrays["signed_credit"], 0.02)
        self.assertAlmostEqual(summary["pds"]["lower_95"], 0.01)

    def test_pre_null_gate_requires_broad_positive_gain(self) -> None:
        summary = {
            "macro": {"pds": 0.004, "signed_credit": 0.001},
            "families": {
                "pds": {"k562": 0.004, "rpe1": 0.004, "jurkat": 0.004},
                "signed_credit": {
                    "k562": 0.001,
                    "rpe1": 0.001,
                    "jurkat": 0.001,
                },
            },
        }
        bootstrap = {
            "pds": {"lower_95": 0.001},
            "signed_credit": {"lower_95": 0.0001},
        }
        folds = {
            "pds": [0.001, 0.002, 0.003, 0.004, -0.0004],
            "signed_credit": [0.001] * 5,
        }
        concentration = {
            "pds": {"family": 0.4, "target_block": 0.2, "outer_fold": 0.4},
            "signed_credit": {
                "family": 0.4,
                "target_block": 0.2,
                "outer_fold": 0.4,
            },
        }
        self.assertTrue(
            all(
                M.pre_null_checks(
                    summary,
                    bootstrap,
                    folds,
                    concentration,
                ).values()
            )
        )
        summary["families"]["pds"]["rpe1"] = -1e-6
        self.assertFalse(
            M.pre_null_checks(
                summary,
                bootstrap,
                folds,
                concentration,
            )["every_family_pds_nonnegative"]
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
