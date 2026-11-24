from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("experiment_141_run", HERE / "run.py")
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


class CalibrationTests(unittest.TestCase):
    def synthetic(self):
        rng = np.random.default_rng(7)
        n, g = 20, RUN.RESPONSE_COUNT
        base = {
            family: rng.normal(size=(n, g))
            for family in RUN.FAMILIES
        }
        residual_source = rng.normal(size=(n, g))
        pan = {
            family: base[family] + residual_source
            for family in RUN.FAMILIES
        }
        coefficient = np.linspace(0.0, 1.0, g)
        truth = {
            family: base["k562"] + residual_source * coefficient
            for family in RUN.FAMILIES
        }
        return base, pan, truth, coefficient

    def test_diagonal_fit_recovers_transfer(self):
        base, pan, truth, expected = self.synthetic()
        observed = RUN.fit_calibration(
            "k562",
            base,
            pan,
            truth,
            RUN.FAMILIES,
            np.arange(20),
            1.0,
        )
        self.assertTrue(np.allclose(observed, expected, atol=1e-10))

    def test_scalar_is_single_coefficient(self):
        base, pan, truth, _ = self.synthetic()
        observed = RUN.fit_calibration(
            "k562",
            base,
            pan,
            truth,
            RUN.FAMILIES,
            np.arange(20),
            0.0,
        )
        self.assertTrue(np.allclose(observed, observed[0]))
        self.assertGreaterEqual(observed[0], 0.0)
        self.assertLessEqual(observed[0], 1.0)

    def test_negative_coefficients_clip_to_zero(self):
        rng = np.random.default_rng(8)
        base = {
            family: np.zeros((4, RUN.RESPONSE_COUNT))
            for family in RUN.FAMILIES
        }
        pan = {
            family: rng.normal(size=(4, RUN.RESPONSE_COUNT))
            for family in RUN.FAMILIES
        }
        truth = {family: -pan["k562"] for family in RUN.FAMILIES}
        observed = RUN.fit_calibration(
            "k562",
            base,
            pan,
            truth,
            RUN.FAMILIES,
            np.arange(4),
            1.0,
        )
        self.assertTrue(np.array_equal(observed, np.zeros_like(observed)))

    def test_primary_is_fixed_before_outcomes(self):
        self.assertEqual(RUN.PRIMARY_CANDIDATE, "diag_rho050")
        self.assertEqual(RUN.CANDIDATES[RUN.PRIMARY_CANDIDATE], 0.5)
        self.assertEqual(RUN.CANDIDATES["scalar_control"], 0.0)

    def test_calibration_stays_in_destination_geometry(self):
        rng = np.random.default_rng(10)
        n, g = 8, RUN.RESPONSE_COUNT
        base = {
            family: np.zeros((n, g))
            for family in RUN.FAMILIES
        }
        destination_q = rng.normal(size=(n, g))
        pan = {
            "k562": destination_q,
            "rpe1": 100.0 * rng.normal(size=(n, g)),
            "jurkat": 100.0 * rng.normal(size=(n, g)),
        }
        expected = np.linspace(0.0, 1.0, g)
        truth = {
            family: destination_q * expected
            for family in RUN.FAMILIES
        }
        observed = RUN.fit_calibration(
            "k562",
            base,
            pan,
            truth,
            ("rpe1", "jurkat"),
            np.arange(n),
            1.0,
        )
        self.assertTrue(np.allclose(observed, expected, atol=1e-10))

    def test_teacher_policy_excludes_source_self_leakage(self):
        self.assertEqual(RUN.CALIBRATION_TEACHERS["k562"], ("rpe1", "jurkat"))
        self.assertEqual(RUN.CALIBRATION_TEACHERS["rpe1"], ("jurkat",))
        self.assertEqual(RUN.CALIBRATION_TEACHERS["jurkat"], ("rpe1",))

    def test_derangements_have_no_fixed_points(self):
        values = RUN.derangements(np.random.default_rng(9), 7, 20)
        self.assertEqual(values.shape, (7, 20))
        self.assertTrue(np.all(values != np.arange(20)[None, :]))

    def test_sealed_split_has_no_loader_or_command(self):
        source = (HERE / "run.py").read_text(encoding="utf-8")
        self.assertNotIn("load_test", source)
        self.assertNotIn('"evaluate"', source)
        self.assertNotIn('"confirm"', source)
        self.assertIn("sealed row in Pan center", source)
        self.assertIn("plan changed after authorization", source)
        self.assertIn("authorize-development", source)
        self.assertIn("DEVELOPMENT_RESULT_FREEZE", source)
        self.assertIn("INFERENCE_RESULT_FREEZE", source)


if __name__ == "__main__":
    unittest.main()
