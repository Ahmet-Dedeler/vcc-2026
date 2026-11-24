import numpy as np
from core341 import diagonal, fit, predict, sign_endpoints


def test_diagonal_known_answer():
    x = np.array([[-1., -1.], [1., 1.]])
    np.testing.assert_allclose(diagonal(x, 3*x), 2.)


def test_ridge_known_answer():
    x = np.array([[-1., 0.], [1., 0.]])
    model = fit(x, 3*x, gain=np.ones(2))
    np.testing.assert_allclose(model["beta"], [[1., 0.], [0., 0.]])
    np.testing.assert_allclose(predict(model, x), 2*x)


def test_identity_and_zero():
    x = np.random.default_rng(341).normal(size=(80, 12))
    model = fit(x, x, gain=np.ones(12))
    np.testing.assert_array_equal(model["beta"], 0.)
    np.testing.assert_array_equal(predict(model, x), x)
    np.testing.assert_array_equal(predict(model, np.zeros_like(x)), 0.)


def test_subset_linearity():
    rng = np.random.default_rng(342)
    x, y = rng.normal(size=(80, 12)), rng.normal(size=(80, 12))
    model = fit(x, y)
    sources = rng.normal(size=(5, 7, 12))
    expected = predict(model, sources.mean(0))
    measured = np.mean([predict(model, np.delete(sources, j, axis=0).mean(0))
                        for j in range(5)], axis=0)
    np.testing.assert_allclose(expected, measured, atol=1e-14, rtol=1e-14)


def test_sign_is_target_macro():
    truth = np.array([[1., -3.], [100., -100.]])
    prediction = np.array([[1., 0.], [100., -100.]])
    got = sign_endpoints(prediction, truth)
    np.testing.assert_allclose(got["weighted_sign"], [.25, 1.])
    assert got["weighted_sign"].mean() == .625


def test_top64_stable_tie():
    truth = np.ones((1, 65))
    prediction = np.ones((1, 65))
    prediction[0, -1] = -1
    got = sign_endpoints(prediction, truth)
    assert got["top64_weighted_sign"][0] == 1.
    assert got["weighted_sign"][0] == 64/65


def test_held_target_isolation():
    rng = np.random.default_rng(345)
    x, y = rng.normal(size=(80, 12)), rng.normal(size=(80, 12))
    first = fit(x[:40], y[:40])
    y[40:] = 1e9
    second = fit(x[:40], y[:40])
    for key in first:
        np.testing.assert_array_equal(first[key], second[key])


def test_inputs_unchanged():
    rng = np.random.default_rng(346)
    x, y = rng.normal(size=(80, 12)), rng.normal(size=(80, 12))
    xx, yy = x.copy(), y.copy()
    predict(fit(x, y), x)
    np.testing.assert_array_equal(x, xx)
    np.testing.assert_array_equal(y, yy)


def test_whole_target_derangement_is_consistent():
    from run341 import P
    permutation = P.derangement(20, 3410929)
    assert np.array_equal(np.sort(permutation), np.arange(20))
    assert np.all(permutation != np.arange(20))
    targets = np.arange(60).reshape(20, 3)
    labels = np.tile(targets[permutation], (5, 1))
    for view in range(5):
        np.testing.assert_array_equal(labels[20*view:20*(view+1)], targets[permutation])
