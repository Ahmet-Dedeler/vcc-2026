"""Availability is a mask, not a test for a numerically nonzero effect."""

import numpy as np

NAMES = ("cap15","available","count_only","shuffled")


def estimate(base,sources,observed,operator_mask,seed=3370909):
    base = np.asarray(base,dtype=float)
    source = np.asarray(sources,dtype=float)
    observed = np.asarray(observed,dtype=bool)
    op = np.asarray(operator_mask,dtype=bool)
    assert source.shape == observed.shape and source.shape[1:] == base.shape
    assert source.shape[0] == 3 and op.shape == (base.shape[1],)
    assert np.isfinite(base).all() and np.isfinite(source).all()
    have = observed.any(2)
    norms = np.linalg.norm(source,axis=2)
    weights = np.zeros_like(norms)
    for d in range(3):
        median = float(np.median(norms[d,have[d]])) if have[d].any() else 0.
        weights[d,have[d]] = 1.
        good = have[d] & (norms[d]>0) & (median>0)
        weights[d,good] = np.maximum(1,norms[d,good]/median)**-1.5
    total = weights.sum(0)
    normalized = np.divide(weights,total[None,:],out=np.zeros_like(weights),
                           where=total[None,:]>0)
    mass = np.einsum("dt,dtg->tg",normalized,observed)
    count = observed.sum(0)
    carrying = have.sum(0)
    partial = op[None,:] & (count>=2) & (count<carrying[:,None])
    factor = np.ones_like(base)
    factor[partial] = np.minimum(4,1/np.maximum(mass[partial],1e-12))
    count_factor = np.ones_like(base)
    numerator = np.broadcast_to(carrying[:,None],base.shape)
    count_factor[partial] = numerator[partial]/count[partial]
    shuffled = factor.copy()
    rng = np.random.default_rng(seed)
    eligible_shuffle = op[None,:] & (count>=2)
    for t in range(len(base)):
        indices = np.flatnonzero(eligible_shuffle[t])
        if len(indices)>1:
            cycle = rng.permutation(indices)
            shuffled[t,cycle] = factor[t,np.roll(cycle,1)]
    eligible = np.broadcast_to(op,base.shape).copy()
    baseline_norm = np.linalg.norm(base[:,op],axis=1)
    outputs = {"cap15":base.copy()}
    for name,f in (("available",factor),("count_only",count_factor),("shuffled",shuffled)):
        out = base*f
        norm = np.linalg.norm(out[:,op],axis=1)
        for t in range(len(base)):
            if baseline_norm[t]>0 and norm[t]>0:
                out[t,op] *= baseline_norm[t]/norm[t]
            else:
                out[t] = base[t]
        np.testing.assert_array_equal(out[~eligible],base[~eligible])
        np.testing.assert_allclose(np.linalg.norm(out[:,op],axis=1),baseline_norm,
                                   rtol=1e-13,atol=1e-13)
        outputs[name] = out
    details = {"policy":"available cap-weight mass; max factor 4; operator norm restored",
               "partial_support_coordinates":int(partial.sum()),
               "partial_support_targets":int(partial.any(1).sum()),
               "factor_min_median_max":np.quantile(factor[partial],[0,.5,1]).tolist()
                                        if partial.any() else [1.,1.,1.],
               "factor_cap_coordinates":int((factor[partial]==4).sum()),
               "normalized_donor_weights":normalized.T.tolist(),
               "before_restoration_full_support_identity":True,
               "after_restoration_full_support_identity_claimed":False,
               "gene_counts_at_least_two":int(eligible_shuffle.sum())}
    return outputs,eligible,details
