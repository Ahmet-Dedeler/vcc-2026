"""Independent prediction readback and paired interval verification."""

import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE/"work"
OLD = ROOT/"experiments/334-text-conditioned-transfer/recovery-work"


def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(8<<20),b""):
            h.update(block)
    return h.hexdigest()


def main():
    output = WORK/"independent-audit.json"
    if output.exists():
        raise FileExistsError(output)
    result = json.loads((WORK/"results.json").read_text())
    frozen = json.loads((HERE/"inputs.json").read_text())
    for rel,expected in frozen["sha256"].items():
        assert sha(ROOT/rel)==expected,rel
    rows = list(csv.DictReader((WORK/"per-target.csv").open()))
    keyed = {(r["line"],r["target"],r["arm"]):r for r in rows}
    assert len(keyed)==len(rows)==6*1024*11
    oldrows = list(csv.DictReader((OLD/"per-target.csv").open()))
    metrics = ("absolute_error","squared_error","cosine","retrieval","decoy_cosine")
    for r in oldrows:
        saved = keyed[r["line"],r["target"],r["arm"]]
        for m in metrics:
            assert float(saved[m])==float(r[m])
    newarms = ("semantic","semantic_shuffled","semantic_labels_shuffled")
    lines = ["k562_gwps","orion_hct116","orion_hek293t","rpe1","nadig_jurkat","nadig_hepg2"]
    maxerr = {m:0. for m in metrics}
    values = {a:[] for a in newarms}
    for line in lines:
        with np.load(OLD/line/"evaluation.npz") as z:
            truth,targets = z["truth"],list(z["targets"])
        outputs = {a:np.full_like(truth,np.nan) for a in newarms}
        for b in (0,1):
            receipt = json.loads((WORK/line/f"fit-{b}.json").read_text())
            for name,ext in (("prediction","npz"),("models","pickle")):
                assert sha(WORK/line/f"{name}-{b}.{ext}")==receipt[f"{name}_sha256"]
            assert receipt["source_only_replay_max_abs"]==0.
            with np.load(WORK/line/f"prediction-{b}.npz") as z:
                for a in newarms:
                    outputs[a][z["test_indices"]]=z[a]
        tn = truth/np.maximum(np.linalg.norm(truth,axis=1)[:,None],1e-12)
        off = ~np.eye(len(truth),dtype=bool)
        for a,p in outputs.items():
            pn = p/np.maximum(np.linalg.norm(p,axis=1)[:,None],1e-12)
            s = np.einsum("ig,jg->ij",pn,tn,optimize=True)
            own = s.diagonal().copy()
            ret = (((s<own[:,None])&off).sum(1)+.5*((s==own[:,None])&off).sum(1))/(len(truth)-1)
            measured = (np.abs(p-truth).mean(1),np.einsum("tg,tg->t",p-truth,p-truth)/512,
                        own,ret,np.where(off,s,0).sum(1)/(len(truth)-1))
            values[a].append(ret)
            for m,v in zip(metrics,measured):
                error = float(np.max(np.abs(v-[float(keyed[line,t,a][m]) for t in targets])))
                maxerr[m]=max(maxerr[m],error)
                assert error<5e-12,(a,m,error)
    interval_errors={}
    for comparator in ("source_mean","source_only","text","semantic_shuffled","diagonal"):
        other=np.array([[float(keyed[line,t,comparator]["retrieval"]) for t in targets] for line in lines])
        delta=np.array(values["semantic"])-other
        rng=np.random.default_rng(3340909)
        target_draws,crossed_draws=[],[]
        for _ in range(100):
            ti=rng.integers(1024,size=(100,1024))
            li=rng.integers(6,size=(100,6))
            tw=np.array([np.bincount(i,minlength=1024) for i in ti])/1024
            lw=np.array([np.bincount(i,minlength=6) for i in li])/6
            byline=tw@delta.T
            target_draws.extend(byline.mean(1))
            crossed_draws.extend(np.einsum("bl,bl->b",byline,lw))
        saved=result["comparisons"][comparator]["retrieval"]
        error=max(float(np.max(np.abs(np.quantile(v,[.025,.975])-saved[k]))) for v,k in
                  ((target_draws,"target_ci95"),(crossed_draws,"crossed_ci95")))
        assert error<5e-12
        interval_errors[comparator]=error
    record={"status":"PASS","endpoint_rows":len(rows),"unchanged_prior_rows":len(oldrows),
            "new_endpoint_rows":len(rows)-len(oldrows),"maximum_metric_differences":maxerr,
            "retrieval_interval_differences":interval_errors,
            "source_replay_maximum_difference":0.,"result_sha256":sha(WORK/"results.json"),
            "scope":"independent saved-output metrics, hashes, prior-row identity and retrieval intervals",
            "not_reproduced":"semantic tree refitting or independent PCA/model prediction replay"}
    output.write_text(json.dumps(record,indent=2,sort_keys=True)+"\n")
    print(json.dumps(record,indent=2))


if __name__=="__main__":
    with threadpool_limits(limits=6):
        main()
