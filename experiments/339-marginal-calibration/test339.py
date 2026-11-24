import numpy as np
from core339 import fit_curves, predict, curve_map


def test_zero_anchor():
    x = np.arange(24).reshape(2, 6, 2)-8
    curves = fit_curves(x)
    for mode in ("marginal", "gaussian", "clipping"):
        np.testing.assert_array_equal(predict(curves, np.zeros((2, 3, 2)), mode), 0.)


def test_identical_curves_identity():
    x = np.tile(np.array([-2., -1., 0., 1., 2.])[None, :, None], (3, 1, 2))
    query = np.array([[[-1., .3], [.5, 1.5]]]*3)
    np.testing.assert_allclose(predict(fit_curves(x), query), query.mean(0), atol=1e-15)


def test_analytic_scale_transport():
    values = np.array([-2., -1., 0., 1., 2.])
    x = np.array([values, 3*values])[..., None]
    query = np.array([[[1.]], [[3.]]])
    for mode in ("marginal", "gaussian"):
        np.testing.assert_allclose(predict(fit_curves(x), query, mode), [[2.]])


def test_ties_are_averaged():
    result = curve_map(np.array([-1., 0., .5, 1.]),
                       np.array([-1., 0., 0., 1.]),
                       np.array([-2., 0., 2., 4.]))
    np.testing.assert_allclose(result, [-3., 0., 1.5, 3.])


def test_identical_curves_equals_clipping():
    x = np.tile(np.array([-2., -1., 0., 1., 2.])[None, :, None], (3, 1, 1))
    query = np.array([[[-4.], [-.5], [3.]]]*3)
    curves = fit_curves(x)
    np.testing.assert_allclose(predict(curves, query),
                               predict(curves, query, "clipping"), atol=1e-15)


def test_constant_curve_identity():
    x = np.ones((3, 5, 2))
    query = np.arange(12).reshape(3, 2, 2)-5
    for mode in ("marginal", "gaussian", "clipping"):
        np.testing.assert_allclose(predict(fit_curves(x), query, mode), query.mean(0))


def test_training_pairing_invariance():
    rng = np.random.default_rng(22)
    x = rng.normal(size=(3, 12, 4))
    shuffled = np.stack([row[rng.permutation(12)] for row in x])
    a, b = fit_curves(x), fit_curves(shuffled)
    np.testing.assert_array_equal(a["sorted"], b["sorted"])
    np.testing.assert_array_equal(a["barycenter"], b["barycenter"])
    np.testing.assert_allclose(a["sd"], b["sd"], atol=1e-15)


def test_query_does_not_refit():
    rng = np.random.default_rng(23)
    x, query = rng.normal(size=(3, 12, 4)), rng.normal(size=(3, 5, 4))
    curves = fit_curves(x)
    saved = {k: v.copy() for k, v in curves.items()}
    first = predict(curves, query)
    huge = np.concatenate((query, np.full((3, 5, 4), 1e8)), axis=1)
    np.testing.assert_array_equal(first, predict(curves, huge)[:5])
    for key in curves:
        np.testing.assert_array_equal(saved[key], curves[key])


def test_curve_gene_permutation():
    rng = np.random.default_rng(24)
    x, query = rng.normal(size=(3, 12, 4)), rng.normal(size=(3, 5, 4))
    permutation = np.array([2, 3, 0, 1])
    curves = fit_curves(x)
    manual = {k: v[..., permutation] for k, v in curves.items()}
    np.testing.assert_array_equal(predict(curves, query, gene_permutation=permutation),
                                  predict(manual, query))


def test_rejects_nonfinite():
    try:
        fit_curves(np.full((3, 3, 3), np.nan))
    except ValueError:
        return
    raise AssertionError("accepted nonfinite training")
