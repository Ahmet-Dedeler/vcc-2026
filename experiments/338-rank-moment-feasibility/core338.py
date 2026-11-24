"""Continuous per-gene histogram bounds, not a joint integer emitter."""

import numpy as np
from scipy.optimize import linprog
from scipy.stats import rankdata


def phi(control_sorted,counts,depth):
    values=np.log1p(np.asarray(counts,dtype=float)/(np.asarray(depth,dtype=float)/1e6))
    return (np.searchsorted(control_sorted,values,side="left")+
            np.searchsorted(control_sorted,values,side="right"))/(2*len(control_sorted))


def u_known(control_sorted,counts,depths):
    values=np.log1p(np.asarray(counts,dtype=float)/(np.asarray(depths,dtype=float)/1e6))
    m,n=len(control_sorted),len(values)
    ranks=rankdata(np.r_[control_sorted,values],method="average")
    direct=ranks[m:].sum()-n*(n+1)/2
    lookup=phi(control_sorted,counts,depths).sum()*m
    _,tie_counts=np.unique(np.r_[control_sorted,values],return_counts=True)
    N=m+n
    tie=float(np.sum(tie_counts.astype(float)**3-tie_counts))
    variance=n*m*(N+1)/12*(1-tie/(N**3-N))
    z=(direct-n*m/2)/np.sqrt(variance) if variance>0 else 0.
    return float(direct),float(lookup),float(z)


def bounds(control_sorted,original_counts,original_depths,new_depths,rank_constraint=True):
    original_counts=np.asarray(original_counts)
    original_depths=np.asarray(original_depths)
    new_depths=np.asarray(new_depths)
    control_sorted=np.asarray(control_sorted)
    if (control_sorted.ndim!=1 or not len(control_sorted) or
        not np.all(np.isfinite(control_sorted)) or np.any(np.diff(control_sorted)<0)):
        raise ValueError("finite nonempty sorted control required")
    for values in (original_counts,original_depths,new_depths):
        if (values.ndim!=1 or not np.all(np.isfinite(values)) or
            np.any(values!=np.rint(values))):
            raise ValueError("finite integer vector required")
    if (original_counts.shape!=original_depths.shape or
        original_counts.shape!=new_depths.shape or not len(original_counts) or
        np.any(original_counts<0) or np.any(original_depths<=0) or
        np.any(new_depths<=0) or np.any(original_counts>original_depths)):
        raise ValueError("invalid count/depth geometry")
    counts=np.asarray(original_counts,dtype=np.int64)
    rows=np.asarray(new_depths,dtype=np.int64)
    assert counts.ndim==rows.ndim==1 and len(counts)==len(rows)
    assert np.all(counts>=0) and np.all(rows>0)
    n=len(rows)
    total=int(counts.sum())
    zero=int((counts==0).sum())
    target_phi=float(phi(control_sorted,counts,original_depths).sum())
    depths,groups=np.unique(rows,return_counts=True)
    K=min(total,int(rows.min()),max(64,4*int(counts.max(initial=0))))
    levels=np.arange(K+1)
    group_id=np.repeat(np.arange(len(groups)),K+1)
    k=np.tile(levels,len(groups))
    lib=depths[group_id]
    matrix=[(group_id==j).astype(float) for j in range(len(groups))]
    rhs=list(map(float,groups))
    # Normalize the count equation to avoid a large coefficient range.
    count_scale=max(1,K)
    matrix.append(k/count_scale); rhs.append(total/count_scale)
    matrix.append((k==0).astype(float)); rhs.append(float(zero))
    f=phi(control_sorted,k,lib)
    if rank_constraint:
        matrix.append(f); rhs.append(target_phi)
    A,b=np.array(matrix),np.array(rhs)
    objective=(k/(lib/1e6))/n
    scale=max(float(np.mean(counts/(np.asarray(original_depths)/1e6))),1.)
    results=[]
    for sign in (1,-1):
        result=linprog(sign*objective/scale,A_eq=A,b_eq=b,bounds=(0,None),
                       method="highs",options={"time_limit":2.,
                       "primal_feasibility_tolerance":1e-9,
                       "dual_feasibility_tolerance":1e-9})
        record={"status":int(result.status),"message":str(result.message)}
        if result.success:
            residual=float(np.max(np.abs(A@result.x-b)))
            minimum=float(np.min(result.x))
            count_error=float(abs(k@result.x-total))
            phi_error=float(abs(f@result.x-target_phi))
            dual=np.asarray(result.eqlin.marginals)
            reduced=sign*objective/scale-A.T@dual
            duality_gap=float(sign*objective@result.x/scale-b@dual)
            if (residual>1e-6 or minimum< -1e-7 or count_error>1e-4 or
                (rank_constraint and phi_error>1e-6) or
                np.min(reduced)< -1e-7 or abs(duality_gap)>1e-6):
                raise ArithmeticError(f"LP residual {residual}, minimum {minimum}")
            active=np.flatnonzero(result.x>1e-8)
            record.update(mean_cpm=float(objective@result.x),
                          equality_max_abs=residual,minimum_variable=minimum,
                          count_total_absolute_error=count_error,
                          phi_sum_absolute_error=phi_error,
                          dual_equality_marginals=dual.tolist(),
                          dual_objective_scaled=float(b@dual),
                          minimum_reduced_cost=float(np.min(reduced)),
                          duality_gap_scaled=duality_gap,
                          histogram=[{"depth":int(lib[i]),"count":int(k[i]),
                                      "multiplicity":float(result.x[i])} for i in active])
        results.append(record)
    return {"minimum":results[0],"maximum":results[1],"count_support_max":K,
            "depth_groups":[{"depth":int(d),"cells":int(c)} for d,c in zip(depths,groups)],
            "baseline_mean_cpm":float(np.mean(counts/(np.asarray(original_depths)/1e6))),
            "baseline_gene_total":total,"baseline_zero_cells":zero,
            "baseline_phi_sum":target_phi,"rank_sum_constrained":rank_constraint,
            "continuous_multiplicities":True}
