"""Small deterministic fixtures; no biological response data."""

import numpy as np
import core334 as C


def test_outer_isolation():
    rng = np.random.default_rng(1)
    a, rms = C.rms_direction(rng.normal(size=(6, 8, 20)))
    basal = rng.normal(size=(6, 20))
    train = np.arange(4)
    old = C.source_blocks(a, basal, rms, 0, train)
    query = C.query_blocks(a, basal, rms, 0, np.arange(4, 8))
    a[0], rms[0] = 1e9, 1e9
    new = C.source_blocks(a, basal, rms, 0, train)
    for (r,x,y),(r2,x2,y2) in zip(old,new):
        assert r == r2 and r != 0
        np.testing.assert_array_equal(x,x2)
        np.testing.assert_array_equal(y,y2)
    for x,x2 in zip(query,C.query_blocks(a,basal,rms,0,np.arange(4,8))):
        np.testing.assert_array_equal(x,x2)


def test_held_target_label_isolation():
    rng = np.random.default_rng(2)
    a, rms = C.rms_direction(rng.normal(size=(6, 8, 20)))
    basal = rng.normal(size=(6, 20))
    before = C.source_blocks(a, basal, rms, 0, np.arange(4))
    a[:, 4:] = 999
    rms[:, 4:] = 999
    after = C.source_blocks(a, basal, rms, 0, np.arange(4))
    for (_,x,y),(_,x2,y2) in zip(before,after):
        np.testing.assert_array_equal(x,x2)
        np.testing.assert_array_equal(y,y2)


def test_balanced_unique():
    ix = C.sample_indices(5, 20, 30, 2000, 1)
    assert len(set(map(tuple,ix))) == 2000
    assert np.ptp(np.bincount(ix[:,0])) <= 1
    assert ix[:,1].max() < 20 and ix[:,2].max() < 30


def test_subset_mean():
    rng = np.random.default_rng(3)
    a,r = C.rms_direction(rng.normal(size=(6, 8, 20)))
    b = rng.normal(size=(6,20))
    q = C.query_blocks(a,b,r,0,np.arange(8))
    np.testing.assert_allclose(np.mean([x[...,0] for x in q],0),a[1:].mean(0),atol=1e-14)


def test_derangement():
    ix = C.derangement(50,1)
    assert not (ix == np.arange(50)).any()
    np.testing.assert_array_equal(np.sort(ix),np.arange(50))


def test_metrics():
    a = np.eye(6)
    np.testing.assert_allclose(C.endpoints(a,a)["retrieval"],1)
    np.testing.assert_allclose(C.endpoints(np.zeros_like(a),a)["retrieval"],.5)
    np.testing.assert_allclose(C.endpoints(a,a)["squared_error"],0)


def test_diagonal():
    x = np.ones((4,3,len(C.SOURCE_FEATURES)))
    y = 3*np.ones((4,3))
    np.testing.assert_array_equal(C.diagonal([(0,x,y),(1,x,y)]),2)


def test_description():
    d = {"ABC":"Gene Symbol ABC a protein", "DEF":"Gene Symbol XYZ another",
         "GHI":"Gene Symbol GHI2 wrong boundary"}
    assert C.description(d,"ABC") == "a protein"
    assert C.description(d,"DEF") == C.description(d,"GHI") == ""
    assert C.description(d,"abc") == ""


def test_held_vocabulary():
    docs = [f"shared molecule mechanism{i%20} process{i%8}" for i in range(40)]
    docs.append("zebrazone nevertraining")
    e,v,s = C.text_representation(docs,np.arange(40),1)
    assert "zebrazone" not in v.vocabulary_
    assert not e[-1].any()
    assert np.isfinite(e).all()


def test_synthetic_learning():
    rng = np.random.default_rng(4)
    x = rng.uniform(-1,1,(4096,3))
    y = x[:,0]*x[:,1]
    m = C.model(4)
    m.fit(x[:3072],y[:3072])
    error = np.mean((m.predict(x[3072:])-y[3072:])**2)
    assert error < .04
