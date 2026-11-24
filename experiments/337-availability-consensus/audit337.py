"""Rescore sufficient statistics; independently form the dense weighted MSE budget."""

import csv
import json
import numpy as np
from threadpoolctl import threadpool_limits
import run337 as A
from cell_eval2.metrics.delta import PRED_TRACE_CAP_K
from cell_eval2.prep import bulk_lognorm_means


def dense_mse(weights,p,t,ctrl,pjk,rjk,cjk,own):
    """Direct gene-wise weighted sums, no bootstrap Gram/target-patch shortcut."""
    n,g=p.shape
    a=weights.astype(float)/(g-1)
    total=a.sum(1)
    sums=a@p
    denominators=np.broadcast_to(total[:,None],sums.shape).copy()
    for i,col in enumerate(own):
        sums[:,col]-=a[:,i]*p[i,col]
        denominators[:,col]-=a[:,i]
    centering=np.divide(sums*sums,denominators,out=np.zeros_like(sums),
                        where=denominators>0).sum(1)
    square=p*p
    square[np.arange(n),own]=0
    budget=np.maximum(0,a@square.sum(1)-centering)
    claim=a@np.minimum(pjk,PRED_TRACE_CAP_K*rjk)
    rho=np.minimum(1,np.divide(budget,claim,out=np.ones_like(claim),where=claim>0))
    error=(p-t)**2
    error[np.arange(n),own]=0
    energy=(t-ctrl)**2
    energy[np.arange(n),own]=0
    num=a@(error.sum(1)-rjk)-rho*claim
    den=a@(energy.sum(1)-rjk-cjk)
    assert np.all(den>0)
    return num/den


def main():
    A.verify()
    out=A.WORK/"independent-audit.json"
    if out.exists():
        raise FileExistsError(out)
    result=json.loads((A.WORK/"results.json").read_text())
    composition=json.loads((A.WORK/"compose.json").read_text())
    assert A.R.sha(A.WORK/"effects.npz")==composition["effects_sha256"]
    r,c,emit_control=A.E.load_reference()
    reference=bulk_lognorm_means(r["sums"].astype(float),A.E.TS)
    truth,control=reference[:-1],reference[-1]
    lookup={g:i for i,g in enumerate(r["genes"])}
    own=np.array([lookup[t] for t in r["targets"]])
    targets=list(r["targets"])
    members=list(A.E.MEMBERS)
    anchors=json.loads((A.BASE/"work/h1-anchors.json").read_text())
    base=np.array([anchors["baseline"][m] for m in members])
    replicate=np.array([anchors["replicate"][m] for m in members])
    rng=np.random.default_rng(20260909)
    indices=rng.integers(len(targets),size=(10000,len(targets)))
    weights=np.array([np.bincount(row,minlength=len(targets)) for row in indices])
    weights=np.vstack((np.ones(len(targets),dtype=int),weights))
    raw_by_arm,overall={},{}
    maximum_raw,maximum_rows={},{}
    for name in A.NAMES:
        record=json.loads((A.WORK/f"metrics-{name}.json").read_text())
        path=A.WORK/f"stats-{name}.npz"
        assert A.R.sha(path)==record["statistics_sha256"]
        with np.load(path) as z:
            stats={k:z[k] for k in z.files}
        recomputed,rows,legs=A.E.score(r,c,stats["sums"],stats["jk"],
                                      stats["p"],stats["q"],stats["lfc"])
        maximum_raw[name]=max(abs(recomputed[m]-record["raw"][m]) for m in members)
        assert maximum_raw[name]<1e-12
        saved=list(csv.DictReader((A.WORK/f"per-target-{name}.csv").open()))
        index={(x["perturbation"],x["metric"]):float(x["value"]) for x in saved}
        differences=[]
        for row in rows:
            old=index[row["perturbation"],row["metric"]]
            new=row["value"]
            if np.isnan(old) and np.isnan(new):
                continue
            differences.append(abs(old-new))
        maximum_rows[name]=max(differences)
        assert maximum_rows[name]<1e-12
        p=bulk_lognorm_means(stats["sums"].astype(float),A.E.TS)
        raw=np.empty((len(weights),len(members)))
        for start in range(0,len(weights),64):
            raw[start:start+64,1]=dense_mse(weights[start:start+64],p,truth,control,
                                           stats["jk"],r["jk"][:-1],r["jk"][-1],own)
        for j,m in enumerate(members):
            if m=="mse":
                continue
            values=np.array([index.get((str(t),A.E.MEMBERS[m]),np.nan) for t in targets])
            valid=np.isfinite(values)
            raw[:,j]=(weights[:,valid]@values[valid])/weights[:,valid].sum(1)
        assert np.max(np.abs(raw[0]-[record["raw"][m] for m in members]))<1e-9
        scaled=(raw-base)/(replicate-base)
        scaled[:,1]=np.clip(scaled[:,1],0,1)
        scaled[:,2]=np.maximum(scaled[:,2],-6)
        raw_by_arm[name]=raw
        overall[name]=scaled.mean(1)
        print(f"audited {name}",flush=True)
    differences={}
    for name in ("cap15","count_only","shuffled"):
        expected=result["comparisons"][name]
        delta=overall["available"]-overall[name]
        vals=[abs(float(delta[0])-expected["overall"]),
              *np.abs(np.quantile(delta[1:],[.025,.975])-
                       [expected["overall_ci"]["lo"],expected["overall_ci"]["hi"]])]
        for j,m in enumerate(members):
            d=raw_by_arm["available"][:,j]-raw_by_arm[name][:,j]
            e=expected["raw"][m]
            vals.extend([abs(d[0]-e["point"]),
                         *np.abs(np.quantile(d[1:],[.025,.975])-[e["lo"],e["hi"]])])
        differences[name]=float(max(vals))
        assert differences[name]<1e-9,(name,differences[name])
    A.R.write_json(out,{"status":"PASS","six_metric_replay_max_abs":maximum_raw,
                        "per_target_replay_max_abs":maximum_rows,
                        "independent_bootstrap_max_abs":differences,
                        "result_sha256":A.R.sha(A.WORK/"results.json"),
                        "scope":"sufficient-statistic scorer replay plus independent dense weighted MSE bootstrap",
                        "not_reproduced":"raw count regeneration or independent source-factor reconstruction"})


if __name__=="__main__":
    with threadpool_limits(limits=6):
        main()
