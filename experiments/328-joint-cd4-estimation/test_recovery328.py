import unittest

import numpy as np

import recover328 as A


class ContinuationTests(unittest.TestCase):
    def test_converged_model_is_unchanged(self):
        ll = np.array([[0., -1.], [-1., 0.]])
        initial = A.M.fit(ll)
        self.assertTrue(initial["converged"])
        self.assertEqual(A.continue_em(ll, initial), initial)

    def test_continuation_is_same_em_trajectory(self):
        rng = np.random.default_rng(328)
        ll = rng.normal(size=(200, 5))
        initial = A.M.fit(ll, max_iter=2)
        self.assertFalse(initial["converged"])
        resumed = A.continue_em(ll, initial)
        straight = A.M.fit(ll, max_iter=2000)
        self.assertEqual(resumed["weights"], straight["weights"])
        self.assertEqual(resumed["mean_loglik_trace"], straight["mean_loglik_trace"])

    def test_wrong_endpoint_fails(self):
        ll = np.array([[0., -1.], [-1., 0.]])
        initial = A.M.fit(ll)
        initial["mean_loglik_trace"][-1] += .01
        with self.assertRaises(AssertionError):
            A.continue_em(ll, initial)


if __name__ == "__main__":
    unittest.main()
