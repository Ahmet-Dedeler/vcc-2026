from __future__ import annotations

import importlib.util
import json
import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


freeze_plans = load_module("experiment_125_freeze_plans", HERE / "freeze_plans.py")
sys.modules["freeze_plans"] = freeze_plans
run = load_module("experiment_125_run", HERE / "run.py")


def synthetic_prepared(seed: int = 125, n: int = 300):
    rng = np.random.default_rng(seed)
    latent = rng.normal(scale=0.7, size=(n, n))
    y0 = latent + rng.normal(scale=0.2, size=(n, n))
    y1 = latent + rng.normal(scale=0.8, size=(n, n))
    valid = np.ones((2, n, n), dtype=bool)
    diagonal = np.arange(n)
    y0[diagonal, diagonal] = 0.0
    y1[diagonal, diagonal] = 0.0
    valid[:, diagonal, diagonal] = False
    effects = np.stack((y0, y1))
    scales = np.ones((2, n), dtype=np.float64)
    return run.PreparedFamilies(
        effects=effects,
        valid=valid,
        normalized=effects.copy(),
        row_scales=scales,
        row_origin=np.arange(n),
    )


class FreezeTests(unittest.TestCase):
    def test_plan_is_frozen_and_valid(self) -> None:
        manifest = freeze_plans.verify()
        self.assertEqual(manifest["target_folds"], 5)
        self.assertEqual(manifest["null_draws_each"], 200)
        self.assertEqual(manifest["bootstrap_replicates"], 2000)

    def test_target_folds_are_deterministic_and_balanced(self) -> None:
        panel = tuple(f"G{i}" for i in range(300))
        first = freeze_plans.stable_target_folds(panel, 5)
        second = freeze_plans.stable_target_folds(panel, 5)
        self.assertTrue(np.array_equal(first, second))
        self.assertTrue(np.array_equal(np.bincount(first), np.full(5, 60)))


class PosteriorTests(unittest.TestCase):
    def test_covariance_decomposition_detects_noisier_family(self) -> None:
        prepared = synthetic_prepared()
        fit = np.ones(300, dtype=bool)
        forbidden = np.zeros(300, dtype=bool)
        moments = run.fit_raw_moments(prepared, fit, forbidden)
        self.assertGreater(moments.global_noise1, moments.global_noise0)
        self.assertGreater(moments.global_latent, 0.0)

    def test_forbidden_target_values_cannot_change_moments(self) -> None:
        prepared = synthetic_prepared()
        forbidden = np.zeros(300, dtype=bool)
        forbidden[:60] = True
        fit = ~forbidden
        first = run.fit_raw_moments(prepared, fit, forbidden)
        changed_effects = prepared.effects.copy()
        changed_effects[:, :60] = 1e9
        changed = run.PreparedFamilies(
            effects=changed_effects,
            valid=prepared.valid,
            normalized=changed_effects,
            row_scales=prepared.row_scales,
            row_origin=prepared.row_origin,
        )
        second = run.fit_raw_moments(changed, fit, forbidden)
        self.assertTrue(np.array_equal(first.prior_mean_raw, second.prior_mean_raw))
        self.assertTrue(np.array_equal(first.latent_raw, second.latent_raw))
        self.assertEqual(first.global_total, second.global_total)

    def test_permuted_biological_origin_is_excluded(self) -> None:
        prepared = synthetic_prepared()
        permutation = np.roll(np.arange(300), 60)
        permuted = run.PreparedFamilies(
            effects=prepared.effects[:, permutation],
            valid=prepared.valid[:, permutation],
            normalized=prepared.normalized[:, permutation],
            row_scales=prepared.row_scales[:, permutation],
            row_origin=permutation,
        )
        forbidden = np.zeros(300, dtype=bool)
        forbidden[:60] = True
        fit = np.ones(300, dtype=bool)
        moments = run.fit_raw_moments(permuted, fit, forbidden)
        self.assertEqual(moments.fit_biological_rows, 240)

    def test_effect_and_confidence_are_distinct(self) -> None:
        prepared = synthetic_prepared()
        moments = run.fit_raw_moments(
            prepared,
            np.ones(300, dtype=bool),
            np.zeros(300, dtype=bool),
        )
        setting = run.PosteriorSetting("x", 40.0, 0.05, 1.0)
        params = run.posterior_parameters(moments, setting)
        effect, confidence, support = run.posterior_prediction(
            prepared,
            (0, 1),
            params,
        )
        self.assertFalse(np.array_equal(np.abs(effect), confidence))
        self.assertTrue(np.all((confidence >= 0.0) & (confidence <= 1.0)))
        self.assertTrue(np.all(np.diag(effect) == 0.0))
        self.assertTrue(np.all(np.diag(support) == 0))


class TransportTests(unittest.TestCase):
    def test_repair_handles_sub_ulp_total_residual(self) -> None:
        values = np.linspace(0.0, 0.002257820922051579, 223)
        values *= 0.026444077993670705 / run.math.fsum(map(float, values))
        target = run.math.nextafter(
            run.math.fsum(map(float, values)),
            run.math.inf,
        )
        run.repair_to_sum(values, target, np.ones(values.size, dtype=bool))
        self.assertEqual(run.math.fsum(map(float, values)), target)
        self.assertTrue(np.all(values >= 0.0))

    def test_both_transports_preserve_required_invariants(self) -> None:
        q = np.linspace(1e-6, 4e-4, 300)
        q /= q.sum() * 5.0
        valid = np.ones(300, dtype=bool)
        residual = np.linspace(-2.0, 1.5, 300)
        support = np.ones(300, dtype=np.int16)
        for arm in run.TRANSPORT_NAMES:
            effect, output_support, audit = run.transport_row(
                q,
                valid,
                residual,
                support,
                17,
                arm,
            )
            self.assertTrue(all(audit["invariants"].values()))
            self.assertEqual(effect[17], 0.0)
            self.assertEqual(output_support[17], 0)
            self.assertTrue(np.all(np.isfinite(effect)))

    def test_production_log1p_clips_negative_pretransport_abundance(self) -> None:
        q = np.full(300, 1e-5)
        valid = np.ones(300, dtype=bool)
        residual = np.zeros(300)
        residual[8] = -20.0
        support = np.ones(300, dtype=np.int16)
        _, _, audit = run.transport_row(
            q,
            valid,
            residual,
            support,
            0,
            "production_log1p",
        )
        self.assertGreaterEqual(audit["preclip_negative"], 1)
        self.assertTrue(all(audit["invariants"].values()))


class ComparatorTests(unittest.TestCase):
    def test_equal_mean_agreement_matches_frozen_formula(self) -> None:
        prepared = synthetic_prepared(n=300)
        effect, confidence, support = run.equal_mean_agreement(prepared, (0, 1))
        expected = 0.5 * (
            prepared.normalized[0, 0, 1] + prepared.normalized[1, 0, 1]
        )
        agreement = abs(
            np.sign(prepared.normalized[0, 0, 1])
            + np.sign(prepared.normalized[1, 0, 1])
        ) / 2.0
        self.assertAlmostEqual(effect[0, 1], expected)
        self.assertAlmostEqual(confidence[0, 1], abs(expected) * agreement)
        self.assertEqual(support[0, 1], 2)


class ContractTests(unittest.TestCase):
    def test_protocol_grid_and_exclusions_are_exact(self) -> None:
        protocol = json.loads(
            (HERE / "preregistered_protocol.json").read_text(encoding="utf-8")
        )
        authorization = json.loads(
            (HERE / "SOURCE_ONLY_AUTHORIZATION.json").read_text(encoding="utf-8")
        )
        candidates = run.candidates_from_protocol(protocol)
        self.assertEqual(len(candidates), 8)
        self.assertEqual({candidate.transport for candidate in candidates}, set(run.TRANSPORT_NAMES))
        self.assertEqual(
            tuple(protocol["biological_family_contract"]["order"]),
            run.FAMILY_NAMES,
        )
        for key in (
            "read_pan",
            "read_dld1",
            "read_h1",
            "read_orion_held_harness_outcomes",
            "read_challenge_truth",
            "run_cell_eval2",
            "generate_counts",
            "build_submission",
            "submit",
            "commit",
            "push",
        ):
            self.assertFalse(authorization["authorization"][key])


if __name__ == "__main__":
    unittest.main()
