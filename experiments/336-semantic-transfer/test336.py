import numpy as np
from run336 import C,feature_projection


def test_projection_excludes_held():
    rng = np.random.default_rng(1)
    a = rng.normal(size=(80,32))
    have = np.ones(80,bool)
    train = np.arange(20)
    first,p,fit = feature_projection(a,have,train,40)
    a[20:40] = rng.normal(0,1000,(20,32))
    second,q,fit2 = feature_projection(a,have,train,40)
    np.testing.assert_array_equal(fit,fit2)
    np.testing.assert_allclose(p.components_,q.components_,atol=0,rtol=0)
    np.testing.assert_allclose(first[fit],second[fit],atol=0,rtol=0)


def test_missing_zero():
    a = np.random.default_rng(2).normal(size=(80,32))
    have = np.ones(80,bool)
    have[4] = False
    a[4] = 0
    x,p,fit = feature_projection(a,have,np.arange(20),40)
    assert not x[4].any()


def test_frozen_model():
    m = C.model(123).get_params()
    assert m["loss"] == "squared_error" and m["max_iter"] == 150
    assert m["max_leaf_nodes"] == 15 and m["max_depth"] == 4
    assert m["l2_regularization"] == 10 and not m["early_stopping"]


def test_derangements():
    for n,seed in ((1024,C.SEED+10),(512,C.SEED+11)):
        p = C.derangement(n,seed)
        assert not (p==np.arange(n)).any()
        np.testing.assert_array_equal(np.sort(p),np.arange(n))
