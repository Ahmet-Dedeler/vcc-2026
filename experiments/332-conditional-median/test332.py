import unittest
import numpy as np
import core332 as C


class DesignTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(332)
        self.r = rng.normal(size=(5, 20, 32))
        self.b = rng.uniform(0, 3, size=(5, 32))
        self.rms = np.sqrt(np.mean(self.r**2, axis=2))

    def test_outer_poison_does_not_change_features_or_labels(self):
        first = C.fold_features(self.r, self.b, self.rms, 2)
        self.r[2], self.rms[2] = 1e9, 1e9
        second = C.fold_features(self.r, self.b, self.rms, 2)
        for aa, bb in zip(first[0], second[0]):
            for a, b in zip(aa[1:], bb[1:]):
                np.testing.assert_array_equal(a, b)
        for a, b in zip(first[1:], second[1:]):
            np.testing.assert_array_equal(a, b)

    def test_source_sets_exclude_recipient_and_outer(self):
        self.assertEqual(C.source_indices(5, 2, 3), [0, 1, 4])

    def test_no_recipient_ablation_and_wrong_control(self):
        _, x, wrong, _ = C.fold_features(self.r, self.b, self.rms, 2)
        np.testing.assert_array_equal(x[..., :10], wrong[..., :10])
        self.assertGreater(np.max(np.abs(x[..., 10:] - wrong[..., 10:])), .1)

    def test_whole_target_split_and_balanced_reproducible_sample(self):
        names = [f"t{i}" for i in range(20)]
        a, b = C.target_blocks(names)
        self.assertFalse(set(a) & set(b))
        self.assertEqual(set(a) | set(b), set(range(20)))
        blocks = C.fold_features(self.r, self.b, self.rms, 2)[0]
        first = C.sample_training(blocks, a, 332, 64)
        second = C.sample_training(blocks, a, 332, 64)
        for x, y in zip(first, second):
            np.testing.assert_array_equal(x, y)
        provenance = first[-1]
        self.assertNotIn(2, provenance[:, 0])
        self.assertFalse(set(provenance[:, 1]) & set(b))
        self.assertEqual(len(set(map(tuple, provenance))), 64)
        self.assertEqual(np.unique(provenance[:, 0], return_counts=True)[1].tolist(), [16]*4)

    def test_identity_zero_and_ties(self):
        truth = np.eye(8)
        np.testing.assert_allclose(C.retrieval(truth, truth), 1)
        np.testing.assert_allclose(C.retrieval(truth*0, truth), .5)
        ties = np.ones((8, 8))
        np.testing.assert_allclose(C.retrieval(ties, ties), .5)
        out = C.endpoints(truth, truth, truth, np.ones(8))
        np.testing.assert_array_equal(out["absolute_error"], np.zeros(8))
        np.testing.assert_array_equal(out["squared_error"], np.zeros(8))
        np.testing.assert_array_equal(out["orthogonal_cosine"], np.zeros(8))

    def test_standardization_uses_training_only(self):
        x = np.arange(140).reshape(10, 14).astype(float)
        a, _ = C.standardize(x, np.zeros((2, 14)))
        b, _ = C.standardize(x, np.full((2, 14), 1e8))
        np.testing.assert_array_equal(a, b)
        np.testing.assert_allclose(a.mean(0), 0, atol=1e-7)

    def test_crossed_bootstrap_and_zero_contrast(self):
        result = C.paired_interval(np.zeros((4, 8)), n=100)
        self.assertEqual(result["target_ci95"], [0, 0])
        self.assertEqual(result["crossed_ci95"], [0, 0])
        self.assertEqual(result["minimum_leave_one_line_out"], 0)

    def test_predict_single_output_list_returns_array(self):
        pred = C.prediction_columns(np.arange(8), ["median"], 8)
        self.assertEqual(set(pred), {"median"})
        both = C.prediction_columns({"mean": np.arange(8), "median": np.arange(8)},
                                    ["mean", "median"], 8)
        self.assertEqual(set(both), {"median", "mean"})
        with self.assertRaises(ArithmeticError):
            C.prediction_columns(np.array([np.nan]), ["median"], 1)

    def test_skewed_fixture_distinguishes_mean_and_median_loss(self):
        values = np.array([0., 0., 0., 8.])
        med, mean = np.median(values), values.mean()
        self.assertLess(np.abs(values-med).mean(), np.abs(values-mean).mean())
        self.assertLess(np.square(values-mean).mean(), np.square(values-med).mean())

    def test_half_array_keeps_target_gene_axes(self):
        array = np.zeros((12, 82, 1000))
        selected = np.arange(512)
        self.assertEqual(array[2][:, selected].shape, (82, 512))


if __name__ == "__main__":
    unittest.main()
