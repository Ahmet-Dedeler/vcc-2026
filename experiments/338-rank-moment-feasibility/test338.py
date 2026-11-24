import numpy as np
from core338 import bounds,phi,u_known


def test_constant_depth_mean():
    c=np.zeros(8)
    x=np.array([10,10,10,10])
    r=np.full(4,100)
    got=bounds(c,x,r,r)
    assert got["minimum"]["status"]==got["maximum"]["status"]==0
    np.testing.assert_allclose([got["minimum"]["mean_cpm"],got["maximum"]["mean_cpm"]],1e5)


def test_variable_known_room():
    c=np.zeros(8)
    x=np.array([10,10,10,10])
    got=bounds(c,x,np.full(4,100),np.array([50,50,150,150]))
    assert got["minimum"]["status"]==got["maximum"]["status"]==0
    np.testing.assert_allclose(got["minimum"]["mean_cpm"],(2/50+38/150)/4*1e6)
    np.testing.assert_allclose(got["maximum"]["mean_cpm"],(38/50+2/150)/4*1e6)


def test_zero_stays_zero():
    c=np.zeros(8)
    got=bounds(c,np.zeros(4,dtype=int),np.full(4,100),np.array([50,50,150,150]))
    assert got["minimum"]["mean_cpm"]==got["maximum"]["mean_cpm"]==0.


def test_midranks_with_ties():
    c=np.log1p(np.array([0.,0.,100.,100.,1000.]))
    x=np.array([0,1,1,2])
    d=np.array([10000,10000,10000,20000])
    direct,lookup,z=u_known(c,x,d)
    assert direct==lookup and np.isfinite(z)
    assert phi(c,np.array([0]),10000)[0]==.2


def test_lp_residuals():
    got=bounds(np.zeros(8),np.array([1,2,3,4]),np.full(4,100),np.array([75,75,125,125]))
    for key in ("minimum","maximum"):
        assert got[key]["status"]==0
        assert got[key]["equality_max_abs"]<1e-8
        assert got[key]["minimum_variable"]>=-1e-9
        assert got[key]["count_total_absolute_error"]<1e-8
        assert abs(got[key]["duality_gap_scaled"])<1e-8


def test_saved_histogram_constraints():
    control=np.log1p(np.array([0.,0.,1000.,10000.]))
    counts=np.array([0,1,2,3])
    old=np.full(4,1000)
    new=np.array([750,750,1250,1250])
    got=bounds(control,counts,old,new)
    for key in ("minimum","maximum"):
        assert got[key]["status"]==0
        hist=got[key]["histogram"]
        assert abs(sum(h["count"]*h["multiplicity"] for h in hist)-6)<1e-7
        assert abs(sum(h["multiplicity"] for h in hist if h["count"]==0)-1)<1e-7
        for depth in (750,1250):
            assert abs(sum(h["multiplicity"] for h in hist if h["depth"]==depth)-2)<1e-7
        rank_sum=sum(float(phi(control,h["count"],h["depth"]))*h["multiplicity"] for h in hist)
        assert abs(rank_sum-float(phi(control,counts,old).sum()))<1e-7
        mean=sum(h["count"]/(h["depth"]/1e6)*h["multiplicity"] for h in hist)/4
        assert abs(mean-got[key]["mean_cpm"])<1e-7


def test_rejects_fractional_counts():
    try:
        bounds(np.zeros(8),np.array([1.,1.5]),np.array([100,100]),np.array([75,125]))
    except ValueError:
        return
    raise AssertionError("fractional counts silently truncated")


def test_unconstrained_contains_constrained():
    control=np.log1p(np.array([0.,0.,1000.,10000.]))
    args=(control,np.array([0,1,2,3]),np.full(4,1000),np.array([750,750,1250,1250]))
    strict=bounds(*args)
    free=bounds(*args,rank_constraint=False)
    assert free["minimum"]["mean_cpm"]<=strict["minimum"]["mean_cpm"]+1e-7
    assert free["maximum"]["mean_cpm"]>=strict["maximum"]["mean_cpm"]-1e-7
