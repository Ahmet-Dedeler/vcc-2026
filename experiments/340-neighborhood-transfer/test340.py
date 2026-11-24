import numpy as np
from core340 import K, feature_columns, neighbors_from_similarity, random_neighbors, fit, predict


def fixture():
    rng = np.random.default_rng(340)
    x = rng.normal(size=(100, 12))
    graph = random_neighbors(np.full(12, K))
    return x, graph


def test_sparse_padding():
    x, graph = fixture()
    graph[0, 2:] = -1
    z = feature_columns(x, 0, graph)
    np.testing.assert_array_equal(z[:, 0], x[:, 0])
    np.testing.assert_array_equal(z[:, 1:3], x[:, graph[0, :2]])
    np.testing.assert_array_equal(z[:, 3:], 0.)


def test_graph_no_self_or_duplicates():
    x, _ = fixture()
    graph = neighbors_from_similarity(x.T@x)
    for i, row in enumerate(graph):
        assert i not in row and len(set(row)) == K


def test_missing_graph_rows():
    available = np.array([False]+[True]*11)
    graph = neighbors_from_similarity(np.eye(12), available=available, positive=True)
    np.testing.assert_array_equal(graph, -1)


def test_zero_residual_is_diagonal():
    x, graph = fixture()
    gain = np.arange(12)/10
    model = fit(x, x*gain, gain, graph)
    np.testing.assert_array_equal(model["beta"], 0.)
    np.testing.assert_array_equal(predict(model, x), x*gain)


def test_ridge_closed_form():
    x = np.array([[-1., 0.], [1., 0.]])
    graph = np.full((2, K), -1, dtype=int)
    model = fit(x, 3*x, np.ones(2), graph)
    np.testing.assert_allclose(model["beta"][0, 0], 1.)
    np.testing.assert_allclose(predict(model, x), 2*x)


def test_zero_query_and_linearity():
    x, graph = fixture()
    model = fit(x, x[:, ::-1], np.ones(12), graph)
    np.testing.assert_array_equal(predict(model, np.zeros_like(x)), 0.)
    rng = np.random.default_rng(341)
    five = rng.normal(size=(5, 13, 12))
    subset = np.mean([predict(model, np.delete(five, j, axis=0).mean(0))
                      for j in range(5)], axis=0)
    np.testing.assert_allclose(subset, predict(model, five.mean(0)), atol=1e-14, rtol=1e-14)


def test_held_data_not_part_of_fit():
    x, graph = fixture()
    storage = np.stack((x, x[:, ::-1], np.full_like(x, 1e6)))
    before = fit(storage[0], storage[1], np.ones(12), graph)
    storage[2] *= 1000
    after = fit(storage[0], storage[1], np.ones(12), graph)
    for key in before:
        np.testing.assert_array_equal(before[key], after[key])


def test_input_immutability():
    x, graph = fixture()
    saved, gsave = x.copy(), graph.copy()
    model = fit(x, x[:, ::-1], np.ones(12), graph)
    predict(model, x)
    np.testing.assert_array_equal(saved, x)
    np.testing.assert_array_equal(gsave, graph)


def test_neighbor_signal():
    rng = np.random.default_rng(344)
    x = rng.normal(size=(1000, 12))
    graph = np.full((12, K), -1, dtype=int)
    graph[0, 0] = 1
    y = x.copy()
    y[:, 0] += 2*x[:, 1]
    model = fit(x, y, np.ones(12), graph)
    assert np.mean((predict(model, x)[:, 0]-y[:, 0])**2) < np.mean((x[:, 0]-y[:, 0])**2)/2


def test_random_degree_match():
    degree = np.arange(12) % (K+1)
    graph = random_neighbors(degree)
    np.testing.assert_array_equal((graph >= 0).sum(1), degree)
    for gene, row in enumerate(graph):
        valid = row[row >= 0]
        assert gene not in valid and len(set(valid)) == len(valid)
