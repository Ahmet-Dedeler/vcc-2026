from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge


MODULE_PATH = Path(__file__).with_name("run.py")
SPEC = importlib.util.spec_from_file_location("experiment132_run", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


class Experiment132Tests(unittest.TestCase):
    def test_pds_identity_and_zero_prediction(self) -> None:
        truth = np.eye(5, dtype=np.float64)
        np.testing.assert_allclose(RUN.pds_per_target(truth, truth), 1.0)
        np.testing.assert_allclose(RUN.pds_per_target(np.zeros_like(truth), truth), 0.5)

    def test_signed_credit_rewards_correct_direction(self) -> None:
        truth = np.asarray([[2.0, -3.0, 0.0]])
        good = np.asarray([[4.0, -1.0, 0.0]])
        bad = -good
        self.assertEqual(RUN.signed_credit_per_target(good, truth)[0], 1.0)
        self.assertEqual(RUN.signed_credit_per_target(bad, truth)[0], -1.0)

    def test_collapsed_ridge_matches_repeated_upstream_fit(self) -> None:
        rng = np.random.default_rng(7)
        features = rng.normal(size=(31, 14)).astype(np.float32)
        sources = [rng.normal(size=(31, 9)).astype(np.float32) for _ in range(3)]
        response_mean = np.mean(np.stack(sources), axis=0)

        state = np.random.get_state()
        np.random.seed(17)
        try:
            repeated_features = np.tile(features, (3, 1))
            pca = PCA(7).fit(repeated_features)
        finally:
            np.random.set_state(state)
        upstream = Ridge(alpha=10).fit(
            pca.transform(repeated_features),
            np.vstack(sources),
        ).predict(pca.transform(features))

        collapsed = RUN.fit_projection(
            features,
            n_sources=3,
            n_components=7,
            alpha=10,
            seed=17,
        ).predict(response_mean)
        np.testing.assert_allclose(collapsed, upstream, rtol=2e-5, atol=2e-5)

    def test_permuted_projection_matches_direct_repeated_fit(self) -> None:
        rng = np.random.default_rng(11)
        features = rng.normal(size=(27, 12)).astype(np.float32)
        sources = [rng.normal(size=(27, 6)).astype(np.float32) for _ in range(4)]
        response_mean = np.mean(np.stack(sources), axis=0)
        permutation = rng.permutation(len(features))

        state = np.random.get_state()
        np.random.seed(23)
        try:
            assigned = features[permutation]
            repeated_features = np.tile(assigned, (4, 1))
            pca = PCA(6).fit(repeated_features)
        finally:
            np.random.set_state(state)
        direct = Ridge(alpha=10).fit(
            pca.transform(repeated_features),
            np.vstack(sources),
        ).predict(pca.transform(assigned))

        projected = RUN.fit_projection(
            features,
            n_sources=4,
            n_components=6,
            alpha=10,
            seed=23,
        ).predict(response_mean, permutation)
        np.testing.assert_allclose(projected, direct, rtol=3e-5, atol=3e-5)

    def test_outer_context_exclusion_has_five_sources(self) -> None:
        sample = {
            name: np.full((3, 2), index, dtype=np.float64)
            for index, name in enumerate(RUN.CONTEXTS)
        }
        means = RUN.source_means_by_held(sample)
        for held in RUN.HELD_GROUPS:
            expected = np.mean(
                np.stack([sample[name] for name in RUN.CONTEXTS if name != held]),
                axis=0,
            )
            np.testing.assert_array_equal(means[held], expected)

    def test_selection_fails_on_one_negative_context(self) -> None:
        contexts = {
            name: {"delta": {"pds": 0.01, "signed_credit": 0.01}}
            for name in RUN.SCORED_DATASETS
        }
        contexts["rpe1"]["delta"]["pds"] = -1e-6
        summary = {
            "contexts": contexts,
            "macro_delta": {"pds": 0.01, "signed_credit": 0.01},
        }
        nulls = {"pds": np.zeros(49), "signed_credit": np.zeros(49)}
        checks = RUN.selection_passes(summary, nulls)
        self.assertFalse(checks["pds_nonnegative_every_dataset"])
        self.assertFalse(all(checks.values()))

    def test_block_bootstrap_constant_is_exact(self) -> None:
        values = np.full(8, 0.125)
        labels = np.asarray([0, 0, 1, 2, 2, 3, 4, 4])
        draws = np.asarray([[0, 1, 2, 3, 4], [4, 4, 4, 4, 4]])
        np.testing.assert_allclose(RUN.block_bootstrap(values, labels, draws), 0.125)


if __name__ == "__main__":
    unittest.main()
