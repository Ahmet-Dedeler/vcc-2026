"""Actual-count, all-six-member paired analysis for the frozen policy."""

import json
import numpy as np
import polars as pl
import run337 as A
import report323 as B
from cell_eval2.prep import bulk_lognorm_means
from cell_eval2.metrics.delta import _across_pert_budget


def report():
    A.verify()
    if (A.WORK/"results.json").exists():
        raise FileExistsError("terminal report exists")
    assert json.loads((A.WORK/"known-answer.json").read_text())["pass"]
    with np.load(A.BASE/"work/reference.npz") as z:
        r={k:z[k] for k in z.files}
    targets,genes=r["targets"],r["genes"]
    logmeans=bulk_lognorm_means(r["sums"].astype(float),A.E.TS)
    truth,control=logmeans[:-1],logmeans[-1]
    lookup={g:i for i,g in enumerate(genes)}
    own=np.array([lookup[t] for t in targets])
    n=len(targets)
    rng=np.random.default_rng(20260909)
    draw=rng.integers(n,size=(10000,n))
    counts=np.array([np.bincount(row,minlength=n) for row in draw])
    counts=np.vstack((np.ones(n,dtype=int),counts))
    members=list(A.E.MEMBERS)
    assert members==["pds","mse","nmae","fid","reach","jac"]
    anchors=json.loads((A.BASE/"work/h1-anchors.json").read_text())
    base=np.array([anchors["baseline"][m] for m in members])
    replicate=np.array([anchors["replicate"][m] for m in members])
    raw_values,totals,rows,checks={},{},[],{}
    emitted_geometry={}
    base_pb=None
    scored=~np.isin(genes,targets)
    for name in A.NAMES:
        metric=json.loads((A.WORK/f"metrics-{name}.json").read_text())
        path=A.WORK/f"stats-{name}.npz"
        assert A.R.sha(path)==metric["statistics_sha256"]
        with np.load(path) as z:
            np.testing.assert_array_equal(z["targets"],targets)
            pm=bulk_lognorm_means(z["sums"].astype(float),A.E.TS)
            jk=z["jk"]
        emitted=(pm-control)[:,scored]
        if base_pb is None:
            base_pb=emitted
        nb=np.linalg.norm(base_pb,axis=1)
        ratio=np.divide(np.linalg.norm(emitted,axis=1),nb,
                        out=np.ones_like(nb),where=nb>0)
        emitted_geometry[name]={
            "row_norm_ratio_min_median_max":np.quantile(ratio,[0,.5,1]).tolist(),
            "relative_frobenius_change":float(np.linalg.norm(emitted-base_pb)/
                                              max(np.linalg.norm(base_pb),1e-12)),
            "nonzero_baseline_sign_changes":int(((np.sign(emitted)!=np.sign(base_pb))&
                                                  (base_pb!=0)).sum())}
        mse,budget,rho,den=B.mse_bootstrap(counts,pm,truth,control,jk,r["jk"][:-1],r["jk"][-1],own)
        expected_budget=_across_pert_budget(pm,dict(enumerate(map(int,own))),
                                           np.full(n,1/(len(genes)-1)))
        assert abs(mse[0]-metric["raw"]["mse"])<1e-10
        assert abs(budget[0]-expected_budget)<1e-10 and np.all(den>0)
        table=pl.read_csv(A.WORK/f"per-target-{name}.csv")
        raw=np.empty((len(counts),6))
        for j,m in enumerate(members):
            if m=="mse":
                raw[:,j]=mse
            else:
                filtered=table.filter(pl.col("metric")==A.E.MEMBERS[m])
                mapping=dict(zip(filtered["perturbation"],filtered["value"]))
                values=np.array([mapping.get(str(t),np.nan) for t in targets])
                valid=np.isfinite(values)
                denominator=counts[:,valid].sum(1)
                assert np.all(denominator>0)
                raw[:,j]=(counts[:,valid]@values[valid])/denominator
            assert abs(raw[0,j]-metric["raw"][m])<1e-12,(name,m)
        assert np.isfinite(raw).all()
        scaled=B.scaled(raw,base,replicate)
        totals[name]=scaled.mean(1)
        raw_values[name]=raw
        rows.append({"arm":name,**metric["raw"],"local_h1_overall":float(totals[name][0]),
                     "local_scaled_mse":float(scaled[0,1])})
        checks[name]={"mse_identity_abs":float(abs(mse[0]-metric["raw"]["mse"])),
                      "budget_identity_abs":float(abs(budget[0]-expected_budget)),
                      "rho_range":[float(rho.min()),float(rho.max())]}
    comparisons={}
    for other in ("cap15","count_only","shuffled"):
        delta=totals["available"]-totals[other]
        raw=raw_values["available"]-raw_values[other]
        comparisons[other]={"overall":float(delta[0]),"overall_ci":B.ci(delta[1:]),
            "raw":{m:{"point":float(raw[0,j]),**B.ci(raw[1:,j])} for j,m in enumerate(members)}}
    primary=comparisons["cap15"]
    gates={"overall_positive":primary["overall_ci"]["lo"]>0,
           "beats_shuffled":comparisons["shuffled"]["overall_ci"]["lo"]>0,
           "nmae_point_improves":primary["raw"]["nmae"]["point"]<0,
           "nmae_upper_bound":primary["raw"]["nmae"]["hi"]<=.005,
           "count_control_not_better":comparisons["count_only"]["overall"]>=0}
    result={"status":"EARNS_FURTHER_DEVELOPMENT_ONLY" if all(gates.values()) else
                      "DOES_NOT_EARN_FURTHER_DEVELOPMENT",
            "gates":gates,"arms":rows,"comparisons":comparisons,"identity_checks":checks,
            "emitted_non_target_geometry":emitted_geometry,
            "scope":"85-target reused H1 development; not independent A/B/C validation",
            "bootstrap":{"draws":10000,"seed":20260909,"unit":"target","fixed_gallery":True,
                         "panel_mse_and_claim_cap_recomputed":True},
            "input_sha256":A.R.sha(A.HERE/"inputs.json"),
            "effects_sha256":A.R.sha(A.WORK/"effects.npz"),
            "reserve_read":False,"submission":False}
    A.R.write_json(A.WORK/"results.json",result)
    print(json.dumps({"status":result["status"],"gates":gates,"comparisons":comparisons},indent=2))
