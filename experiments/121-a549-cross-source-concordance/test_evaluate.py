from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


freeze_plans = load_module("freeze_plans", HERE / "freeze_plans.py")
evaluate = load_module("experiment_121_evaluate", HERE / "evaluate.py")


class ProtocolTests(unittest.TestCase):
    def test_preregistered_contract(self) -> None:
        protocol = json.loads(
            (HERE / "preregistered_protocol.json").read_text(encoding="utf-8")
        )
        authorization = json.loads(
            (HERE / "AUTHORIZATION_MANIFEST.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            protocol["status"],
            "frozen_before_either_source_expression_array_was_opened",
        )
        self.assertEqual(protocol["target_axis"]["count"], 20)
        self.assertTrue(protocol["target_axis"]["complete_gallery_required"])
        self.assertEqual(protocol["go_gate"]["macro_pds_minimum_inclusive"], 0.6)
        self.assertEqual(
            protocol["go_gate"]["threshold_adjustment"],
            "none",
        )
        self.assertTrue(
            authorization["created_before_source_expression_value_access"]
        )
        denied = {
            "open_vcc_held_perturbation_truth",
            "open_existing_outcome_directories",
            "cell_eval2",
            "count_generation",
            "submission_path_access",
            "submission",
            "commit",
            "push",
        }
        for key in denied:
            self.assertFalse(authorization["authorization"][key])


class PlanTests(unittest.TestCase):
    def test_derangements_and_bootstrap_are_deterministic(self) -> None:
        first = freeze_plans.generate_derangements(25, 20, 1212026)
        second = freeze_plans.generate_derangements(25, 20, 1212026)
        self.assertTrue(np.array_equal(first, second))
        self.assertTrue(np.all(first != np.arange(20)[None, :]))
        self.assertTrue(
            np.all(np.sort(first, axis=1) == np.arange(20)[None, :])
        )
        bootstrap = freeze_plans.generate_bootstrap(50, 20, 1212027)
        self.assertEqual(bootstrap.shape, (50, 20))
        self.assertTrue(np.all((bootstrap >= 0) & (bootstrap < 20)))

    def test_synthetic_freeze_refuses_overwrite(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as directory:
            output = Path(directory) / "freeze"
            freeze_plans.build(output)
            freeze_plans.verify(output)
            with self.assertRaises(freeze_plans.FreezeError):
                freeze_plans.build(output)


class MetricTests(unittest.TestCase):
    def test_midrank_and_zero_row_pds(self) -> None:
        similarity = np.asarray(
            [
                [0.8, 0.8, 0.1],
                [0.9, 0.2, 0.1],
                [0.0, 0.0, 0.0],
            ],
            dtype=np.float64,
        )
        scores = evaluate.pds_from_assignment(
            similarity,
            np.arange(3),
            np.asarray([False, False, True]),
            np.asarray([False, False, False]),
        )
        self.assertTrue(np.allclose(scores, [0.75, 0.5, 0.5]))

    def test_identity_retrieval_is_perfect(self) -> None:
        values = np.eye(4, dtype=np.float64)
        similarity, predictor_zero, reference_zero = evaluate.cosine_matrix(
            values,
            values,
        )
        scores = evaluate.pds_from_assignment(
            similarity,
            np.arange(4),
            predictor_zero,
            reference_zero,
        )
        self.assertTrue(np.array_equal(scores, np.ones(4)))

    def test_sign_metrics_use_predictor_magnitude(self) -> None:
        predictor = np.asarray([[10.0, -9.0, 1.0, 0.5, 0.1]])
        reference = np.asarray([[2.0, 3.0, 1.0, -1.0, -1.0]])
        precision, credit, selected = evaluate.sign_metrics(
            predictor,
            reference,
            fraction=0.4,
        )
        self.assertEqual(selected, 2)
        self.assertAlmostEqual(float(precision[0]), 0.5)
        self.assertAlmostEqual(float(credit[0]), 1.0 / 19.0)

    def test_response_mask_is_global_and_excludes_union(self) -> None:
        genes = ("A", "B", "C", "D", "E")
        altair = np.asarray(
            [
                [True, True, True, True, True],
                [True, True, False, True, True],
            ]
        )
        myllia = np.ones((2, 5), dtype=bool)
        mask = evaluate.build_response_mask(
            genes,
            altair,
            myllia,
            panel=("A",),
            targets=("D",),
        )
        self.assertTrue(np.array_equal(mask, [False, True, False, False, True]))

    def test_null_uses_inverse_mapping_in_reverse(self) -> None:
        similarity = np.eye(4, dtype=np.float64)
        permutation = np.asarray([[1, 2, 3, 0]], dtype=np.int16)
        zero = np.zeros(4, dtype=bool)
        null = evaluate.null_macro_pds(
            similarity,
            similarity.T,
            zero,
            zero,
            permutation,
        )
        self.assertEqual(null.shape, (1,))
        self.assertAlmostEqual(float(null[0]), 1.0 / 3.0)


if __name__ == "__main__":
    unittest.main()
