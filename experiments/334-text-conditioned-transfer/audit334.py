"""Independent readback: no import of pilot feature, fit or metric functions."""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE/"recovery-work"
SEED = 3340909
BANKS = ("data/donors/k562_gwps", "data/donors/orion_hct116",
         "data/donors/orion_hek293t", "data/donors/rpe1",
         "data/donors/nadig_jurkat",
         "experiments/166-donor-acquisition/donors/nadig_hepg2")


def digest(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read(p):
    return json.loads(p.read_text())


def main():
    if (WORK/"independent-audit.json").exists():
        raise FileExistsError("independent audit exists")
    selection, result = read(HERE/"selection.json"), read(WORK/"results.json")
    attestation = read(WORK/"recovery-attestation.json")
    assert digest(WORK/"results.json") == attestation["result_sha256"]
    assert digest(HERE/"inputs.json") == attestation["original_inputs_sha256"]
    assert digest(HERE/"recovery-inputs.json") == attestation["recovery_inputs_sha256"]
    for manifest in ("inputs.json", "recovery-inputs.json"):
        for rel, expected in read(HERE/manifest)["sha256"].items():
            assert digest(ROOT/rel) == expected, rel
    assert digest(WORK/"per-target.csv") == result["per_target_sha256"]
    records = list(csv.DictReader((WORK/"per-target.csv").open()))
    keyed = {(r["line"],r["target"],r["arm"]):r for r in records}
    assert len(keyed) == len(records) == 6*1024*8
    targets, genes = selection["targets"], selection["genes"]
    raw, basal = [], []
    for bank in BANKS:
        p = ROOT/bank
        ts = {s:i for i,s in enumerate((p/"perts.txt").read_text().splitlines())}
        gs = {s:i for i,s in enumerate((p/"genes.txt").read_text().splitlines())}
        d = np.load(p/"delta.npy", mmap_mode="r")
        gi = [gs[g] for g in genes]
        raw.append(np.stack([np.asarray(d[ts[t], gi], dtype=np.float64) for t in targets]))
        basal.append(np.asarray(np.load(p/"b_ctrl.npy")[gi], dtype=np.float64))
    raw, basal = np.array(raw), np.array(basal)
    rms = np.sqrt(np.einsum("ltg,ltg->lt", raw, raw)/len(genes))
    direction = raw/np.maximum(rms[...,None],1e-12)
    maximum = {m:0. for m in ("absolute_error","squared_error","cosine",
                              "retrieval","decoy_cosine")}
    max_training_difference, max_baseline_difference = 0., 0.
    metric_values = {arm:{m:[] for m in maximum} for arm in result["means"]}
    receipts = 0
    for outer, bank in enumerate(BANKS):
        line = Path(bank).name
        path = WORK/line
        recipients = [i for i in range(6) if i != outer]
        with np.load(path/"evaluation.npz") as z:
            truth, baseline = z["truth"],z["baseline"]
            np.testing.assert_allclose(truth,direction[outer],atol=2e-14,rtol=2e-14)
            difference = float(np.max(np.abs(baseline-direction[recipients].mean(0))))
            assert difference < 2e-14
            max_baseline_difference = max(max_baseline_difference,difference)
        predictions = {a:np.full_like(truth,np.nan) for a in result["means"]}
        coverage = np.zeros(len(targets),bool)
        for block in (0,1):
            receipt = read(path/f"fit-{block}.json")
            for kind, suffix in (("prediction","npz"),("training","npz"),("models","pickle")):
                assert digest(path/f"{kind}-{block}.{suffix}") == receipt[f"{kind}_sha256"]
            assert digest(path/"evaluation.npz") == receipt["evaluation_sha256"]
            with np.load(path/f"training-{block}.npz") as z:
                x,y,ys = z["x"],z["y"],z["y_shuffled"]
                ix,train,test = z["indices"],z["train_indices"],z["test_indices"]
                assert not set(train)&set(test)
                assert len(set(map(tuple,ix))) == len(ix) == 32768
                assert np.ptp(np.bincount(ix[:,0])) <= 1
                rng = np.random.default_rng(SEED+block+20)
                cycle = rng.permutation(len(train))
                perm = np.empty(len(train),dtype=int)
                perm[cycle] = np.roll(cycle,1)
                rebuilt = np.empty_like(x)
                yy,yys = np.empty_like(y),np.empty_like(ys)
                for j,r in enumerate(recipients):
                    donors = [i for i in recipients if i != r]
                    assert outer not in donors and r not in donors
                    mask = ix[:,0] == j
                    t,g = train[ix[mask,1]],ix[mask,2]
                    values = direction[donors][:,t,g]
                    b = basal[donors][:,g]
                    drms = rms[donors][:,t]
                    consensus = direction[donors].mean(0)
                    consensus_rms = np.sqrt(np.einsum("tg,tg->t",consensus,consensus)/len(genes))
                    rebuilt[mask] = np.column_stack((
                        values.mean(0),np.median(values,axis=0),values.std(0),
                        (values>0).mean(0),np.abs(values).mean(0),values.min(0),values.max(0),
                        b.mean(0),b.std(0),drms.mean(0),drms.std(0),consensus_rms[t]))
                    yy[mask] = direction[r,t,g]
                    yys[mask] = direction[r,train[perm[ix[mask,1]]],g]
                error = max(float(np.max(np.abs(x-rebuilt))),
                            float(np.max(np.abs(y-yy))),float(np.max(np.abs(ys-yys))))
                assert error < 5e-13, (line,block,error)
                max_training_difference = max(max_training_difference,error)
            with np.load(path/f"prediction-{block}.npz") as z:
                indices = z["test_indices"]
                np.testing.assert_array_equal(indices,test)
                assert not coverage[indices].any()
                coverage[indices] = True
                for arm in predictions:
                    predictions[arm][indices] = z[arm]
            receipts += 1
        assert coverage.all()
        norm_truth = truth/np.maximum(np.sqrt(np.einsum("tg,tg->t",truth,truth))[:,None],1e-12)
        for arm,p in predictions.items():
            norm_pred = p/np.maximum(np.sqrt(np.einsum("tg,tg->t",p,p))[:,None],1e-12)
            similarities = np.einsum("ig,jg->ij",norm_pred,norm_truth,optimize=True)
            own = similarities.diagonal().copy()
            off_diagonal = ~np.eye(len(targets),dtype=bool)
            retrieval = ((similarities<own[:,None])&off_diagonal).sum(1)
            retrieval = (retrieval+.5*((similarities==own[:,None])&off_diagonal).sum(1))/(len(targets)-1)
            metrics = {"absolute_error":np.abs(p-truth).mean(1),
                       "squared_error":np.einsum("tg,tg->t",p-truth,p-truth)/len(genes),
                       "cosine":own,"retrieval":retrieval,
                       "decoy_cosine":np.where(off_diagonal,similarities,0).sum(1)/(len(targets)-1)}
            for metric,v in metrics.items():
                saved = np.array([float(keyed[line,t,arm][metric]) for t in targets])
                error = float(np.max(np.abs(v-saved)))
                assert error < 5e-12,(line,arm,metric,error)
                maximum[metric] = max(maximum[metric],error)
                metric_values[arm][metric].append(v)
    interval_errors = {}
    for comparator in ("source_mean","source_only","string","text_shuffled","diagonal"):
        for metric in ("retrieval","cosine","squared_error"):
            delta = np.array(metric_values["text"][metric])-np.array(metric_values[comparator][metric])
            if metric == "squared_error":
                delta = -delta
            rng = np.random.default_rng(SEED)
            target_draws,crossed_draws = [],[]
            for _ in range(100):
                ti = rng.integers(1024,size=(100,1024))
                li = rng.integers(6,size=(100,6))
                tw = np.array([np.bincount(i,minlength=1024) for i in ti])/1024
                lw = np.array([np.bincount(i,minlength=6) for i in li])/6
                perline = tw@delta.T
                target_draws.extend(perline.sum(1)/6)
                crossed_draws.extend(np.einsum("bl,bl->b",perline,lw))
            saved = result["comparisons"][comparator][metric]
            errors = [np.max(np.abs(np.quantile(draws,[.025,.975])-saved[key]))
                      for draws,key in ((target_draws,"target_ci95"),(crossed_draws,"crossed_ci95"))]
            interval_errors[f"{comparator}/{metric}"] = float(max(errors))
            assert max(errors) < 5e-12
    output = {"status":"PASS","endpoint_rows":len(records),"prediction_receipts":receipts,
              "maximum_endpoint_differences":maximum,
              "maximum_training_feature_or_label_difference":max_training_difference,
              "maximum_baseline_difference":max_baseline_difference,
              "interval_differences":interval_errors,
              "source_inputs_and_recovery_hashes":"verified",
              "scope":"independent raw-source training/label reconstruction and saved-prediction metrics",
              "not_reproduced":"refitting trees, full model-prediction replay, text-SVD numerical refit",
              "result_sha256":digest(WORK/"results.json")}
    (WORK/"independent-audit.json").write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
    print(json.dumps(output,indent=2))


if __name__ == "__main__":
    with threadpool_limits(limits=6):
        main()
