import numpy as np
from estimator337 import estimate


def fixture():
    source=np.array([[[1.,2.,3.,4.],[2.,1.,4.,3.]],
                     [[2.,3.,1.,4.],[1.,3.,2.,4.]],
                     [[3.,1.,2.,4.],[3.,2.,1.,4.]]])
    observed=np.ones_like(source,dtype=bool)
    return source.mean(0),source,observed,np.array([True,True,True,False])


def test_full_identity():
    b,s,o,m=fixture()
    out,_,d=estimate(b,s,o,m)
    for a in out.values():
        np.testing.assert_allclose(a,b,atol=1e-14)
    assert d["partial_support_coordinates"]==0


def test_missing_and_true_zero_differ():
    b,s,o,m=fixture()
    s[0,0,0]=0.
    out,_,d=estimate(b,s,o,m)
    assert d["partial_support_coordinates"]==0
    o[0,0,0]=False
    out,_,d=estimate(b,s,o,m)
    assert d["partial_support_coordinates"]==1
    assert not np.allclose(out["available"],b)


def test_one_source_not_directly_amplified():
    b,s,o,m=fixture()
    o[:2,:,:]=False
    s[:2,:,:]=0
    out,_,d=estimate(b,s,o,m)
    for a in out.values():
        np.testing.assert_allclose(a,b,atol=1e-14)


def test_uniform_weights_match_count():
    b,s,o,m=fixture()
    s[:]=1
    s[0,:,0]=0
    o[0,:,0]=False
    out,_,d=estimate(b,s,o,m)
    np.testing.assert_allclose(out["available"],out["count_only"],atol=1e-14)


def test_norm_and_outside():
    b,s,o,m=fixture()
    s[0,0,0]=0; o[0,0,0]=False
    out,eligible,d=estimate(b,s,o,m)
    for a in out.values():
        np.testing.assert_array_equal(a[~eligible],b[~eligible])
        np.testing.assert_allclose(np.linalg.norm(a[:,m],axis=1),np.linalg.norm(b[:,m],axis=1))


def test_deterministic():
    b,s,o,m=fixture()
    s[0,0,0]=0; o[0,0,0]=False
    first,_,_=estimate(b,s,o,m)
    second,_,_=estimate(b,s,o,m)
    np.testing.assert_array_equal(first["shuffled"],second["shuffled"])


def test_bound():
    s=np.ones((3,3,4))
    o=np.ones_like(s,dtype=bool)
    s[1:,0,1]=1000
    s[0,0,0]=0; o[0,0,0]=False
    b=s.mean(0)
    m=np.array([True,True,True,False])
    out,_,d=estimate(b,s,o,m)
    assert d["factor_min_median_max"][-1]==4
    assert d["factor_cap_coordinates"]>0


def test_zero_row():
    b,s,o,m=fixture()
    b[0]=0
    out,_,_=estimate(b,s,o,m)
    for a in out.values():
        assert not a[0].any()
