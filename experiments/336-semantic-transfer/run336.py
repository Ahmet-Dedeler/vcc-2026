"""Change one representation, replay all source-only known answers."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import pickle
import resource
import subprocess
import sys
import time
import numpy as np
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD = ROOT/"experiments/334-text-conditioned-transfer"
FEATURE = ROOT/"experiments/335-biomedical-text-features"
WORK = HERE/"work"
sys.path.insert(0,str(OLD))
import core334 as C
import run334 as R
NEW_ARMS = ("semantic","semantic_shuffled","semantic_labels_shuffled")
ALL_ARMS = C.ARMS+NEW_ARMS


def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(8<<20),b""):
            h.update(block)
    return h.hexdigest()


def read(p):
    return json.loads(p.read_text())


def write(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    q = p.with_suffix(p.suffix+".tmp")
    q.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")
    q.replace(p)


def paths():
    own = [HERE/n for n in ("PREREGISTRATION.md","run336.py","test336.py")]
    prior = [OLD/n for n in ("inputs.json","recovery-inputs.json","selection.json",
                             "core334.py","run334.py","recovery-work/results.json",
                             "recovery-work/independent-audit.json")]
    prior += [FEATURE/n for n in ("inputs.json","work/results.json","work/features.npz")]
    for line in R.LINES:
        prior.append(OLD/"recovery-work"/line/"evaluation.npz")
        for b in (0,1):
            prior += [OLD/"recovery-work"/line/f"{stem}-{b}.{ext}" for stem,ext in
                      (("training","npz"),("prediction","npz"),("fit","json"))]
    return own+prior+[ROOT/"uv.lock"]


def freeze():
    if (HERE/"inputs.json").exists():
        raise FileExistsError("already frozen")
    assert read(FEATURE/"work/results.json")["status"] == "ENGINEERING_FEASIBLE_ONLY"
    assert read(OLD/"recovery-work/independent-audit.json")["status"] == "PASS"
    write(HERE/"inputs.json",{"sha256":{str(p.relative_to(ROOT)):sha(p) for p in paths()},
                             "prior_development_outcomes_known":True})


def verify():
    m = read(HERE/"inputs.json")
    assert set(m["sha256"]) == {str(p.relative_to(ROOT)) for p in paths()}
    for name,expected in m["sha256"].items():
        assert sha(ROOT/name) == expected,name
    # Also recheck raw sources pinned transitively by 334, without changing its tree.
    for manifest in ("inputs.json","recovery-inputs.json"):
        for name,expected in read(OLD/manifest)["sha256"].items():
            assert sha(ROOT/name) == expected,name
    own = [str((HERE/n).relative_to(ROOT)) for n in
           ("PREREGISTRATION.md","run336.py","test336.py","inputs.json")]
    for name in own:
        subprocess.run(["git","ls-files","--error-unmatch",name],cwd=ROOT,check=True,
                       stdout=subprocess.DEVNULL)
    assert not subprocess.check_output(["git","status","--porcelain","--",*own],cwd=ROOT)


def feature_projection(values,have,train,n_targets):
    fit = np.r_[train,np.arange(n_targets,len(values))]
    vectors,pca = C.string_representation(values,have,fit,C.SEED)
    return vectors,pca,fit


def predict(fitted,queries,pairs,n_targets,n_genes):
    total = np.zeros((n_targets,n_genes))
    for q in queries:
        flat = q.reshape(-1,len(C.SOURCE_FEATURES))
        out = np.empty(len(flat))
        for start in range(0,len(flat),8192):
            stop = start+8192
            out[start:stop] = fitted.predict(np.column_stack((flat[start:stop],
                                pairs[start:stop])))+flat[start:stop,0]
        total += out.reshape(n_targets,n_genes)/len(queries)
    return total


def run():
    verify()
    if (WORK/"execution.json").exists() or list(WORK.glob("*/prediction-*.npz")):
        raise FileExistsError("prior run exists")
    start = time.monotonic()
    start_utc = datetime.now(timezone.utc).isoformat()
    selection = read(OLD/"selection.json")
    targets,genes = selection["targets"],selection["genes"]
    fresult = read(FEATURE/"work/results.json")
    assert sha(FEATURE/"work/features.npz") == fresult["features_sha256"]
    old_text = read(OLD/"recovery-work/features.json")["0"]
    assert fresult["description_sha256"] == old_text["description_sha256"]
    with np.load(FEATURE/"work/features.npz") as z:
        assert list(z["genes"]) == targets+genes
        matrix,have = np.asarray(z["matrix"],dtype=np.float64),z["available"]
        assert matrix.shape == (1536,768) and np.isfinite(matrix).all()
        assert have[:1024].sum() == 1021 and have[1024:].sum() == 507
    loaded = [R.load(b,targets,genes) for b in R.BANKS]
    raw,basal = np.array([a for a,b in loaded]),np.array([b for a,b in loaded])
    direction,rms = C.rms_direction(raw)
    del loaded,raw
    receipts = []
    for outer,line in enumerate(R.LINES):
        folder = WORK/line
        folder.mkdir(parents=True,exist_ok=True)
        for block in (0,1):
            if time.monotonic()-start > 2700:
                raise TimeoutError("45-minute budget")
            split_start = time.monotonic()
            previous = OLD/"recovery-work"/line
            prior_receipt = read(previous/f"fit-{block}.json")
            training = previous/f"training-{block}.npz"
            assert sha(training) == prior_receipt["training_sha256"]
            with np.load(training) as z:
                x,y,ys = z["x"],z["y"],z["y_shuffled"]
                ix,train,test = z["indices"],z["train_indices"],z["test_indices"]
            assert not set(train)&set(test)
            vector,pca,fit = feature_projection(matrix,have,train,len(targets))
            assert not set(test)&set(fit)
            nt,ng = len(targets),len(genes)
            tvec,gvec = vector[:nt],vector[nt:]
            pt,pg = C.derangement(nt,C.SEED+10),C.derangement(ng,C.SEED+11)
            variations = {"semantic":(tvec,gvec),"semantic_shuffled":(tvec[pt],gvec[pg]),
                          "semantic_labels_shuffled":(tvec,gvec),
                          "source_replay":(np.zeros_like(tvec),np.zeros_like(gvec))}
            query = C.query_blocks(direction,basal,rms,outer,test)
            predictions,models = {},{}
            for arm,(tv,gv) in variations.items():
                extra = C.pair_features(tv,gv,train[ix[:,1]],ix[:,2])
                response = ys if arm == "semantic_labels_shuffled" else y
                fitted = C.model(C.SEED+100*outer+block)
                fitted.fit(np.column_stack((x,extra)),response-x[:,0])
                ti,gi = np.repeat(test,ng),np.tile(np.arange(ng),len(test))
                pairs = C.pair_features(tv,gv,ti,gi)
                predictions[arm] = predict(fitted,query,pairs,len(test),ng)
                models[arm] = fitted
            with np.load(previous/f"prediction-{block}.npz") as z:
                replay_error = float(np.max(np.abs(predictions["source_replay"]-z["source_only"])))
                assert replay_error <= 1e-12,replay_error
                np.testing.assert_array_equal(test,z["test_indices"])
            out = folder/f"prediction-{block}.npz"
            np.savez_compressed(out,test_indices=test,**predictions)
            model_path = folder/f"models-{block}.pickle"
            with model_path.open("wb") as f:
                pickle.dump({"models":models,"pca":pca,"vectors":vector,"fit_rows":fit},
                            f,protocol=5)
            record = {"line":line,"block":block,"seconds":time.monotonic()-split_start,
                      "source_only_replay_max_abs":replay_error,
                      "prediction_sha256":sha(out),"models_sha256":sha(model_path),
                      "prior_training_sha256":sha(training),"inputs_sha256":sha(HERE/"inputs.json")}
            write(folder/f"fit-{block}.json",record)
            receipts.append(record)
            print(f"{line} block {block} complete; replay error {replay_error:.2g}",flush=True)
    verify()
    write(WORK/"execution.json",{"start_utc":start_utc,
          "finish_utc":datetime.now(timezone.utc).isoformat(),
          "seconds":time.monotonic()-start,"fits":receipts,
          "max_rss_bytes_macos":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    report()


def report():
    if (WORK/"results.json").exists():
        raise FileExistsError("terminal report exists")
    selection = read(OLD/"selection.json")
    table,values = [],{a:{} for a in ALL_ARMS}
    for line in R.LINES:
        old = OLD/"recovery-work"/line
        with np.load(old/"evaluation.npz") as z:
            truth = z["truth"]
        predictions = {a:np.full_like(truth,np.nan) for a in ALL_ARMS}
        coverage = np.zeros(len(truth),bool)
        for b in (0,1):
            receipt = read(WORK/line/f"fit-{b}.json")
            assert sha(WORK/line/f"prediction-{b}.npz") == receipt["prediction_sha256"]
            with np.load(old/f"prediction-{b}.npz") as z:
                ix = z["test_indices"]
                for a in C.ARMS:
                    predictions[a][ix] = z[a]
            with np.load(WORK/line/f"prediction-{b}.npz") as z:
                np.testing.assert_array_equal(ix,z["test_indices"])
                assert not coverage[ix].any()
                coverage[ix] = True
                for a in NEW_ARMS:
                    predictions[a][ix] = z[a]
        assert coverage.all()
        for a,p in predictions.items():
            metric = C.endpoints(p,truth)
            for m,v in metric.items():
                values[a].setdefault(m,[]).append(v)
            for i,t in enumerate(selection["targets"]):
                table.append({"line":line,"target":t,"arm":a,
                              **{m:float(v[i]) for m,v in metric.items()}})
    values = {a:{m:np.array(v) for m,v in ms.items()} for a,ms in values.items()}
    for a,ms in read(OLD/"recovery-work/results.json")["means"].items():
        for m,expected in ms.items():
            assert abs(values[a][m].mean()-expected) < 1e-12,(a,m)
    path = WORK/"per-target.csv"
    with path.open("w",newline="") as f:
        writer = csv.DictWriter(f,fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)
    comparisons = {}
    for a in ALL_ARMS:
        if a == "semantic":
            continue
        comparisons[a] = {}
        for m,v in values["semantic"].items():
            delta = v-values[a][m]
            if m in ("absolute_error","squared_error"):
                delta = -delta
            comparisons[a][m] = C.interval(delta)
    required = ("source_mean","source_only","diagonal","string","text","semantic_shuffled")
    gates = {f"retrieval_beats_{a}":min(comparisons[a]["retrieval"]["target_ci95"][0],
                                     comparisons[a]["retrieval"]["crossed_ci95"][0]) > 0
             for a in required}
    gates["retrieval_point"] = comparisons["source_mean"]["retrieval"]["mean"] >= .005
    gates["cosine_noninferiority"] = all(comparisons[a]["cosine"]["mean"] >= -.005
                                        for a in ("source_mean","diagonal"))
    gates["squared_error"] = comparisons["source_mean"]["squared_error"]["mean"] >= 0
    contrast = comparisons["source_only"]["retrieval"]
    gates["line_robustness"] = (contrast["positive_lines"] >= 5 and
                               contrast["minimum_leave_one_line_out"] > 0)
    depth = {k:float(np.array(contrast["line_means"])[selection[k]].mean())
             for k in ("lower_cell_lines","higher_cell_lines")}
    gates["cell_number_sensitivity"] = all(v>0 for v in depth.values())
    gates["source_replay"] = all(r["source_only_replay_max_abs"] <= 1e-12 for
                                r in read(WORK/"execution.json")["fits"])
    result = {"status":"EARNS_SEPARATE_COUNT_TEST" if all(gates.values()) else
                       "DOES_NOT_EARN_FURTHER_DEVELOPMENT",
              "gates":gates,"comparisons":comparisons,"cell_number_sensitivity":depth,
              "means":{a:{m:float(v.mean()) for m,v in ms.items()} for a,ms in values.items()},
              "per_target_sha256":sha(path),"inputs_sha256":sha(HERE/"inputs.json"),
              "execution_sha256":sha(WORK/"execution.json"),
              "all_six_competition_members_unpriced":True,"reserve_read":False,"submission":False}
    write(WORK/"results.json",result)
    print(json.dumps({"status":result["status"],"gates":gates,
                     "new_arm_means":{a:result["means"][a] for a in NEW_ARMS}},indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action",choices=("freeze","verify","run","report"))
    with threadpool_limits(limits=6):
        globals()[parser.parse_args().action]()
