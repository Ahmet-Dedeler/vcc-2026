"""Control-only relaxed feasibility; no reference perturbation response is loaded."""

import argparse
from datetime import datetime,timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import subprocess
import sys
import time
import numpy as np
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
WORK=HERE/"work"
BASE=ROOT/"experiments/323-current-family-validation"
sys.path.insert(0,str(BASE))
import run323 as R
import evaluate323 as E
from core338 import bounds,u_known

CONTROL_FIELDS=("genes","targets","b_ctrl","libraries","keep","cs","tie","ref_lin",
                "emit_cs","emit_tie","emit_ref_lin","counts")


def order(names,salt):
    return sorted(range(len(names)),key=lambda i:hashlib.sha256((salt+str(names[i])).encode()).hexdigest())


def inputs():
    return [HERE/n for n in ("PREREGISTRATION.md","core338.py","run338.py","test338.py",
                             "selection.json")]+[
        BASE/n for n in ("inputs.json","evaluation-inputs.json","run323.py","evaluate323.py",
                         "exact_de323.py","work/reference.npz","work/reference.json",
                         "work/effects.npz","work/compose.json","work/metrics-cap15.json",
                         "work/stats-cap15.npz")]+[
        ROOT/"experiments/147-dual-moment-emitter/dual_moment.py",
        ROOT/"data/challenge/pert_counts.csv",ROOT/"uv.lock"]


def freeze():
    if (HERE/"inputs.json").exists() or (HERE/"selection.json").exists():
        raise FileExistsError("already frozen")
    import csv
    with np.load(BASE/"work/reference.npz") as z:
        genes=z["genes"].astype(str)
        targets=z["targets"].astype(str)
        keep=z["keep"]
        cpm=z["ref_lin"]
    exclude=set(targets)
    exclude.update(row["target_gene"] for row in csv.DictReader(
        (ROOT/"data/challenge/pert_counts.csv").open()))
    kept=np.flatnonzero(keep)
    eligible=np.array([j for j,i in enumerate(kept) if genes[i] not in exclude])
    q=np.quantile(cpm[eligible],[.25,.5,.75])
    bins=np.searchsorted(q,cpm[eligible],side="right")
    selected=[]
    for quartile in range(4):
        choices=eligible[bins==quartile]
        assert len(choices)>=8
        names=[genes[kept[j]] for j in choices]
        selected += [int(choices[i]) for i in order(names,f"338-gene-{quartile}|")[:8]]
    tids=order(targets,"338-target|")[:4]
    R.write_json(HERE/"selection.json",{
        "target_indices":tids,"targets":[str(targets[i]) for i in tids],
        "control_gene_indices":selected,
        "full_gene_indices":[int(kept[j]) for j in selected],
        "genes":[str(genes[kept[j]]) for j in selected],
        "control_mean_cpm":[float(cpm[j]) for j in selected],
        "control_quartile_boundaries":q.tolist(),"control_fields_loaded":["genes","targets","keep","ref_lin"],
        "reference_perturbation_values_loaded":False})
    R.write_json(HERE/"inputs.json",{"sha256":{str(p.relative_to(ROOT)):R.sha(p) for p in inputs()}})


def verify():
    m=json.loads((HERE/"inputs.json").read_text())
    assert set(m["sha256"])=={str(p.relative_to(ROOT)) for p in inputs()}
    for name,expected in m["sha256"].items():
        assert R.sha(ROOT/name)==expected,name
    own=[str(p.relative_to(ROOT)) for p in inputs() if p.parent==HERE]
    own.append(str((HERE/"inputs.json").relative_to(ROOT)))
    for name in own:
        subprocess.run(["git","ls-files","--error-unmatch",name],cwd=ROOT,check=True,
                       stdout=subprocess.DEVNULL)
    assert not subprocess.check_output(["git","status","--porcelain","--",*own],cwd=ROOT)
    E.verify()


def control_only():
    # NPZ members are loaded lazily. Do not iterate through z.files / materialize truths.
    with np.load(BASE/"work/reference.npz") as z:
        data={name:z[name] for name in CONTROL_FIELDS}
    c=E.rc_emit.RcControlBlock(keep=data["keep"],cs=data["cs"],tie=data["tie"],
                               ref_lin=data["ref_lin"],n_controls=int(data["counts"][-1]),min_cpm=5)
    ec=E.rc_emit.RcControlBlock(keep=data["keep"],cs=data["emit_cs"],tie=data["emit_tie"],
                                ref_lin=data["emit_ref_lin"],n_controls=c.n_controls,min_cpm=5)
    return data,c,ec


def baseline(r,ec,selection):
    with np.load(BASE/"work/effects.npz") as z:
        np.testing.assert_array_equal(z["targets"],r["targets"])
        effects=R.map_axis(z["cap15"],z["genes"],r["genes"])
    with np.load(BASE/"work/stats-cap15.npz") as z:
        expected_sums=z["sums"]
    expected=json.loads((BASE/"work/metrics-cap15.json").read_text())
    b=R.builder()
    pick_rng,emit_rng=b.make_rngs(0)
    captured={}
    digest=hashlib.sha256()
    for i,t in enumerate(r["targets"]):
        profile=b.profile_from_effect(effects[i],r["b_ctrl"],min_cpm=10.)
        picks=pick_rng.choice(len(r["libraries"]),400,replace=False)
        spread=b.emission_spread_vector(effects[i],r["b_ctrl"],min_cpm=10.,
            tail_spread=1.,keep_top=8000,keep_top_frac=None,keep_spread=.25,
            rank_score=None,spread_map="block",spread_lo=0.,spread_hi=.60)
        gate=b.emission_rank_axis(effects[i],r["b_ctrl"],min_cpm=10.)[0]
        block,_=E.rc_emit.rc_hybrid_block(profile,r["libraries"][picks],.1,emit_rng,
            spread=spread,overdispersion=4.,control=ec,key=np.where(gate,effects[i],0.),
            jitter_seed=1000+i,call_set_first=True,call_set=E.rc_emit.CallSet(400,2200,1.9),
            emit_counts=E.emit_counts)
        np.testing.assert_array_equal(block.sum(0),expected_sums[i])
        digest.update(np.asarray(block,dtype="<u4",order="C").tobytes())
        if i in selection["target_indices"]:
            captured[int(i)]=block.copy()
        if (i+1)%20==0:
            print(f"baseline reconstruction {i+1}/85",flush=True)
    assert digest.hexdigest()==expected["emitted_count_digest_sha256"]
    assert set(captured)==set(selection["target_indices"])
    np.savez_compressed(WORK/"selected-blocks.npz",
        targets=np.array(selection["targets"]),
        blocks=np.stack([captured[i] for i in selection["target_indices"]]),genes=r["genes"])
    R.write_json(WORK/"baseline.json",{
        "all_85_column_sums_exact":True,"whole_count_digest_exact":True,
        "emitted_count_digest_sha256":digest.hexdigest(),
        "selected_blocks_sha256":R.sha(WORK/"selected-blocks.npz"),
        "reference_fields_loaded":CONTROL_FIELDS,
        "reference_perturbation_values_loaded":False})
    return captured


def run():
    verify()
    if (WORK/"selected-blocks.npz").exists():
        raise FileExistsError("prior run exists")
    WORK.mkdir(parents=True,exist_ok=True)
    start=time.monotonic()
    started=datetime.now(timezone.utc).isoformat()
    selection=json.loads((HERE/"selection.json").read_text())
    r,c,ec=control_only()
    captured=baseline(r,ec,selection)
    dm=R.load_module("dual_moment338",ROOT/"experiments/147-dual-moment-emitter/dual_moment.py")
    results,checks=[],[]
    for tindex,target in zip(selection["target_indices"],selection["targets"]):
        block=captured[tindex]
        original_depths=block.sum(1)
        new_depths,_=dm.two_depth_row_margins(int(original_depths.sum()),400,.25,seed=3380909)
        for j,g,full in zip(selection["control_gene_indices"],selection["genes"],
                            selection["full_gene_indices"]):
            if time.monotonic()-start>2700:
                raise TimeoutError("45-minute total budget")
            counts=block[:,full]
            direct,looked,z=u_known(c.cs[j],counts,original_depths)
            logged=np.log1p(counts.astype(float)/(original_depths.astype(float)/1e6))
            native=E.rc_emit.wilcoxon_z(np.ascontiguousarray(c.cs[j:j+1]),
                                       c.tie[j:j+1],np.ascontiguousarray(logged[None,:]))[0]
            assert abs(direct-looked)<1e-7 and abs(z-native)<1e-12,(target,g,z,native)
            same=bounds(c.cs[j],counts,original_depths,original_depths)
            variable=bounds(c.cs[j],counts,original_depths,new_depths)
            no_rank=bounds(c.cs[j],counts,original_depths,new_depths,rank_constraint=False)
            known=all(same[k]["status"]==0 for k in ("minimum","maximum"))
            if known:
                v=same["baseline_mean_cpm"]
                tol=1e-5*max(1,v)
                known=(same["minimum"]["mean_cpm"]<=v+tol and
                       same["maximum"]["mean_cpm"]>=v-tol and
                       same["maximum"]["mean_cpm"]-same["minimum"]["mean_cpm"]<=
                       (1e6*counts.sum()/400)*(1/original_depths.min()-1/original_depths.max())+2*tol)
            feasible=all(variable[k]["status"]==0 for k in ("minimum","maximum"))
            room=0.
            if feasible:
                original=variable["baseline_mean_cpm"]
                room=min(original-variable["minimum"]["mean_cpm"],
                         variable["maximum"]["mean_cpm"]-original)/max(original,1.)
            results.append({"target":target,"gene":g,"control_gene_index":j,
                            "same_depth":same,"variable_depth":variable,"without_rank":no_rank,
                            "constrained_extrema_solved":feasible,"symmetric_relative_room":float(room),
                            "room_at_least_five_percent":bool(feasible and room>=.05)})
            checks.append({"target":target,"gene":g,"baseline_lp_pass":bool(known),
                           "U_absolute_difference":float(abs(direct-looked)),
                           "z_absolute_difference":float(abs(z-native))})
        R.write_json(WORK/"partial.json",{"cases":results,"checks":checks,
                                        "new_accuracy_evaluation":False,
                                        "reference_perturbation_values_loaded":False})
        print(f"{target}: all 32 feasibility cases completed",flush=True)
    fractions={t:sum(r["room_at_least_five_percent"] for r in results if r["target"]==t)/32
               for t in selection["targets"]}
    gates={"baseline_and_rank_instruments":all(c["baseline_lp_pass"] for c in checks),
           "at_least_95_percent_extrema_solved":
               sum(r["constrained_extrema_solved"] for r in results)>=.95*128,
           "each_target_at_least_quarter_with_five_percent_room":all(v>=.25 for v in fractions.values())}
    verify()
    R.write_json(WORK/"results.json",{
        "status":"EARNS_INTEGER_ROW_COUPLING_FEASIBILITY_ONLY" if all(gates.values()) else
                 "DOES_NOT_EARN_INTEGER_CONSTRUCTION",
        "gates":gates,"target_room_fractions":fractions,"cases":results,"checks":checks,
        "interpretation":"continuous per-gene fixed-U bounds under finite count support; not an integer emitter",
        "not_established":["joint row feasibility","tie-corrected z invariance","BH call sets",
                           "prediction improvement","competition score"],
        "reference_perturbation_values_loaded":False,"all_six_members_unpriced":True,
        "start_utc":started,"finish_utc":datetime.now(timezone.utc).isoformat(),
        "seconds":time.monotonic()-start,
        "max_rss_bytes_macos":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "inputs_sha256":R.sha(HERE/"inputs.json"),"baseline_sha256":R.sha(WORK/"baseline.json"),
        "submission":False})
    print(json.dumps({"gates":gates,"target_room_fractions":fractions},indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("action",choices=("freeze","verify","run"))
    with threadpool_limits(limits=6):
        globals()[p.parse_args().action]()
