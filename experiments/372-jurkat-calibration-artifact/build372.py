#!/usr/bin/env python3
"""372 step 1 — axes, fit, poison, shuffled, held-out sign check, artifact.

Mirrors 346 fit346.py/effects346.py verbatim (imports core341/core334, same
loader, same rms_direction, same ridge, same blend). Deviations documented in
PLAN.md and RESULTS.md. Writes only into this experiment's work/.
"""
from __future__ import annotations
import csv, hashlib, json, os, resource, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path
os.environ.setdefault("NUMBA_NUM_THREADS","4")
os.environ.setdefault("OMP_NUM_THREADS","4")
os.environ.setdefault("OPENBLAS_NUM_THREADS","4")
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
for _p in (str(ROOT/"experiments/334-text-conditioned-transfer"), str(ROOT/"experiments/341-known-line-calibration")):
    if _p not in sys.path: sys.path.insert(0,_p)
import core334 as P
import core341 as C
SEED=3460910
HELD_SEED=3720910
DONORS=("data/donors/k562_gwps","data/donors/orion_hct116","data/donors/orion_hek293t")
RECIP="data/donors/nadig_jurkat"
LINES=("k562_gwps","orion_hct116","orion_hek293t","nadig_jurkat")

def names(p:Path): return p.read_text().splitlines()
def sha(p:Path):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def bank(path:Path,targets,genes):
    perts,gene_names=names(path/"perts.txt"),names(path/"genes.txt")
    raw=np.load(path/"delta.npy",mmap_mode="r")
    if raw.shape!=(len(perts),len(gene_names)): raise ValueError(f"{path}: delta shape {raw.shape}")
    pi={p:i for i,p in enumerate(perts)}; gi={g:i for i,g in enumerate(gene_names)}
    rows=[pi[t] for t in targets]; cols=[gi[g] for g in genes]
    out=np.array(raw[np.ix_(rows,cols)],dtype=np.float64)
    basal=np.array(np.load(path/"b_ctrl.npy")[cols],dtype=np.float64)
    if not (np.isfinite(out).all() and np.isfinite(basal).all()): raise ValueError(f"{path}: non-finite")
    np.testing.assert_array_equal(out[:3,:3],[[raw[pi[t],gi[g]] for g in genes[:3]] for t in targets[:3]])
    return out,basal

def build_jurkat_axes():
    panel=[r["target_gene"] for r in csv.DictReader((ROOT/"data/challenge/pert_counts.csv").open())]
    panel_order=panel  # pert_counts.csv order
    canonical=[r["gene_name"] for r in csv.DictReader((ROOT/"data/challenge/gene_names.csv").open())]
    donor_perts=[set(names(ROOT/d/"perts.txt")) for d in DONORS]
    donor_genes=[set(names(ROOT/d/"genes.txt")) for d in DONORS]
    j_perts=names(ROOT/RECIP/"perts.txt"); j_genes=names(ROOT/RECIP/"genes.txt")
    j_pert_set=set(j_perts)
    panel_set=set(panel)
    overlap_j_panel=j_pert_set & panel_set
    common4=set.intersection(*donor_genes)&set(j_genes)&set(canonical)
    training=sorted((set.intersection(*donor_perts)&j_pert_set)-panel_set)
    # covered: common minus training symbols minus panel symbols (both as gene symbols)
    excluded_genes=set(training)|panel_set
    covered=sorted(common4-excluded_genes)
    # query = panel targets with full 3-donor consensus, in panel order
    test=[t for t in panel_order if all(t in s for s in donor_perts)]
    missing=[t for t in panel_order if t not in set(test)]
    # how many panel/training symbols intersected the pre-exclusion axis
    panel_on_axis=sorted(common4 & panel_set)
    train_on_axis=sorted(common4 & set(training))
    return dict(panel=panel,panel_order=panel_order,canonical=canonical,training=training,
        covered=covered,test=test,missing=missing,panel_on_axis=panel_on_axis,
        train_on_axis=train_on_axis,common4=common4,j_perts=j_perts,j_genes=j_genes,
        overlap_j_panel=sorted(overlap_j_panel))

def assemble(axes,poison:bool):
    order=axes["training"]+axes["test"]
    genes=axes["covered"]
    raws,basals=[],[]
    for d in DONORS:
        a,b=bank(ROOT/d,order,genes); raws.append(a); basals.append(b)
    # recipient rows: training targets exist; test (panel) targets do not (asserted)
    j_perts_set=set(axes["j_perts"])
    have_train=[t for t in axes["training"] if t in j_perts_set]
    assert len(have_train)==len(axes["training"]), "all training targets must be in Jurkat bank"
    jur_train,_=bank(ROOT/RECIP,axes["training"],genes)
    jur_basal=np.array(np.load(ROOT/RECIP/"b_ctrl.npy")[[names(ROOT/RECIP/'genes.txt').index(g) for g in genes]],dtype=np.float64)
    jur=np.zeros((len(order),len(genes)),dtype=np.float64)
    jur[:len(axes["training"])]=jur_train
    # test rows stay 0 (no Jurkat rows); poison replaces them with noise
    if poison:
        rng=np.random.default_rng(SEED+77)
        jur[len(axes["training"]):]=1e6*rng.standard_normal((len(axes["test"]),len(genes)))
    raws.append(jur); basals.append(jur_basal)
    raw=np.asarray(raws); del raws
    direction,rms=P.rms_direction(raw)
    if poison: rms[3,len(axes["training"]):]=1e6
    return direction,np.asarray(basals),rms,raw

def fit_all(axes,direction,basal,rms):
    n_train=len(axes["training"])
    train=np.arange(n_train); test=np.arange(n_train,direction.shape[1])
    outer=3; donors=[0,1,2]
    baseline=direction[donors].mean(0)
    inner=P.source_blocks(direction,basal,rms,outer,train)
    cross_gain=P.diagonal(inner)
    x=np.concatenate([f[...,0] for _,f,_ in inner])
    cross_y=np.concatenate([r for _,_,r in inner])
    same_y=np.tile(direction[outer,train],(len(donors),1))
    permutation=P.derangement(n_train,SEED)
    shuffled_y=np.tile(direction[outer,train[permutation]],(len(donors),1))
    models={"crossline_full":C.fit(x,cross_y,gain=cross_gain),
            "same_line_full":C.fit(x,same_y),
            "labels_shuffled_full":C.fit(x,shuffled_y)}
    del cross_y,same_y,shuffled_y
    query=baseline[test]
    predictions={"source_mean":query.copy(),
                 "crossline_diag":query*cross_gain,
                 "same_line_diag":query*models["same_line_full"]["gain"]}
    predictions.update({a:C.predict(m,query) for a,m in models.items()})
    subset=np.mean([C.predict(models["same_line_full"],direction[[d for d in donors if d!=omit]][:,test].mean(0)) for omit in donors],axis=0)
    subset_error=float(np.abs(subset-predictions["same_line_full"]).max())
    sc=models["same_line_full"]["scale"]
    beta=models["same_line_full"]["beta"]
    residual=float(np.abs((x/sc).T@((x/sc)@beta)+len(x)*beta-(x/sc).T@(np.tile(direction[outer,train],(len(donors),1))-x*models["same_line_full"]["gain"])).max())
    return models,predictions,dict(cross_gain=cross_gain,permutation=permutation,
        subset_linear_max_difference=subset_error,normal_equation_max_residual=residual,
        n_rows=int(len(x)),n_features=int(x.shape[1]))

def blend_rows(d,p,beta):
    d=np.asarray(d,dtype=np.float64); nd=np.linalg.norm(d)
    if nd<=0: raise ValueError("zero consensus row")
    if beta==0.0: return d.copy()
    p=np.asarray(p,dtype=np.float64); npn=np.linalg.norm(p)
    if npn<=0: raise ValueError("zero prediction")
    m=(1.0-beta)*(d/nd)+beta*(p/npn); nm=np.linalg.norm(m)
    if nm<=1e-9: raise ValueError("degenerate blend")
    return m*(nd/nm)

def heldout(axes,direction,basal,rms):
    """Within-bank held-out: split training targets into halves, fit on half1,
    evaluate weighted sign on half2's own Jurkat truth. Returns gains vs diagonal."""
    n=len(axes["training"])
    idx=np.arange(n)
    # deterministic even/odd split on sorted order
    halves={"evenodd":(idx[0::2],idx[1::2])}
    rng=np.random.default_rng(HELD_SEED)
    perm=rng.permutation(n)
    halves["random"]=(np.sort(perm[:n//2]),np.sort(perm[n//2:]))
    out={}
    for name,(ia,ib) in halves.items():
        donors=[0,1,2]; outer=3
        # build blocks on half A only: need source_blocks restricted to those targets
        trainA=ia
        inner=P.source_blocks(direction,basal,rms,outer,trainA)
        cross_gain=P.diagonal(inner)
        x=np.concatenate([f[...,0] for _,f,_ in inner])
        same_y=np.tile(direction[outer,trainA],(len(donors),1))
        permA=P.derangement(len(trainA),SEED)
        shuf_y=np.tile(direction[outer,trainA[permA]],(len(donors),1))
        m_same=C.fit(x,same_y); m_shuf=C.fit(x,shuf_y)
        # query half B consensus
        baseB=direction[donors][:,ib].mean(0)
        pred_same=C.predict(m_same,baseB)
        pred_diag=baseB*cross_gain
        pred_shuf=C.predict(m_shuf,baseB)
        truth=direction[outer,ib]
        s_same=C.sign_endpoints(pred_same,truth)
        s_diag=C.sign_endpoints(pred_diag,truth)
        s_shuf=C.sign_endpoints(pred_shuf,truth)
        d_same=s_same["weighted_sign"]-s_diag["weighted_sign"]
        d_shuf=s_shuf["weighted_sign"]-s_diag["weighted_sign"]
        # paired target bootstrap CI for the gain
        rng2=np.random.default_rng(HELD_SEED+1)
        draws=10000; bs=[]
        for s in range(0,draws,100):
            nn=min(100,draws-s); ti=rng2.integers(len(ib),size=(nn,len(ib)))
            bs.extend(d_same[ti].mean(1))
        out[name]=dict(n_fit=int(len(ia)),n_eval=int(len(ib)),
            mean_same=float(s_same["weighted_sign"].mean()),mean_diag=float(s_diag["weighted_sign"].mean()),
            mean_shuf=float(s_shuf["weighted_sign"].mean()),
            gain_same_minus_diag=float(d_same.mean()),
            ci95=np.quantile(bs,[.025,.975]).tolist(),
            frac_improved=float((d_same>0).mean()),
            gain_shuf_minus_diag=float(d_shuf.mean()),
            top64_gain=float((s_same["top64_weighted_sign"]-s_diag["top64_weighted_sign"]).mean()))
    return out

def main():
    started,clock=datetime.now(timezone.utc).isoformat(),time.monotonic()
    WORK.mkdir(parents=True,exist_ok=True)
    axes=build_jurkat_axes()
    print(json.dumps({k:(len(v) if isinstance(v,list) else v) for k,v in axes.items() if k not in ("common4",)},indent=2,sort_keys=True),flush=True)
    assert len(axes["overlap_j_panel"])==0, f"Jurkat bank carries panel targets: {axes['overlap_j_panel']}"
    direction,basal,rms,raw=assemble(axes,poison=False)
    print(f"assemble done {time.monotonic()-clock:.1f}s dir{str(direction.shape)}",flush=True)
    models,predictions,receipt=fit_all(axes,direction,basal,rms)
    print(f"fit done subset={receipt['subset_linear_max_difference']:.3e} resid={receipt['normal_equation_max_residual']:.3e} rows={receipt['n_rows']} feat={receipt['n_features']}",flush=True)
    # poison replay
    pd, pb, pr, _=assemble(axes,poison=True)
    np.testing.assert_array_equal(direction[:3],pd[:3]); np.testing.assert_array_equal(rms[:3],pr[:3]); np.testing.assert_array_equal(basal,pb)
    pm, pp, prec=fit_all(axes,pd,pb,pr)
    for arm in models:
        for key in models[arm]: np.testing.assert_array_equal(models[arm][key],pm[arm][key])
    for arm in predictions: np.testing.assert_array_equal(predictions[arm],pp[arm])
    np.testing.assert_array_equal(receipt["cross_gain"],prec["cross_gain"])
    print("poison PASS byte-identical",flush=True)
    del pd,pb,pr,pm,pp
    ho=heldout(axes,direction,basal,rms)
    print(json.dumps(ho,indent=2,sort_keys=True),flush=True)
    # build artifact at beta 0.5
    beta=0.5
    cons=predictions["source_mean"]  # donor consensus query, N x G
    praw=predictions["same_line_full"]
    N,G=cons.shape
    mat=np.empty_like(cons)
    cnorm=np.linalg.norm(cons,axis=1)
    for i in range(N): mat[i]=blend_rows(cons[i],praw[i],beta)
    # sanity
    assert np.isfinite(mat).all() and np.isfinite(praw).all() and np.isfinite(cons).all()
    un=lambda a:a/np.maximum(np.linalg.norm(a,axis=1,keepdims=True),1e-12)
    cos=(un(mat)*un(cons)).sum(1)
    cos_raw=(un(praw)*un(cons)).sum(1)
    print(f"cos blended median {float(np.median(cos)):.4f} min {float(cos.min()):.4f} max {float(cos.max()):.4f}",flush=True)
    print(f"cos raw median {float(np.median(cos_raw)):.4f}",flush=True)
    # save
    (WORK/"fitted-training-targets.txt").write_text("\n".join(axes["training"])+"\n")
    (WORK/"fitted-covered-genes.txt").write_text("\n".join(axes["covered"])+"\n")
    (WORK/"fitted-query-targets.txt").write_text("\n".join(axes["test"])+"\n")
    np.savez_compressed(WORK/"predictions372.npz",genes=np.array(axes["covered"]),test_targets=np.array(axes["test"]),
        training_targets=np.array(axes["training"]),cross_gain=receipt["cross_gain"],
        **{a:v for a,v in predictions.items()})
    for arm,model in models.items(): np.savez_compressed(WORK/f"model-{arm}.npz",**model)
    np.savez(WORK/"ctx_a_calibrated.npz",targets=np.array(axes["test"]),genes=np.array(axes["covered"]),
        matrix=mat,prediction_raw=praw,consensus=cons,consensus_norm=cnorm,
        beta=np.array(0.5),mode=np.array("blend_rows"),context=np.array("A"),norm_rule=np.array("fixed_norm_on_covered"))
    # sidecar
    bank_shas={}
    for d in list(DONORS)+[RECIP]:
        bank_shas[d]=dict(delta=sha(ROOT/d/"delta.npy"),perts=sha(ROOT/d/"perts.txt"),genes=sha(ROOT/d/"genes.txt"),b_ctrl=sha(ROOT/d/"b_ctrl.npy"))
    sidecar=dict(context="A",mode="blend_rows",beta=0.5,norm_rule="fixed_norm_on_covered",
        donor_banks=list(DONORS),recipient_bank=RECIP,
        n_training_targets=len(axes["training"]),n_covered_genes=len(axes["covered"]),
        n_query_targets=len(axes["test"]),n_panel_total=len(axes["panel"]),
        n_panel_missing_no_consensus=len(axes["missing"]),missing_targets=axes["missing"],
        jurkat_panel_overlap=axes["overlap_j_panel"],
        n_panel_genes_dropped_from_axis=len(axes["panel_on_axis"]),panel_genes_dropped=axes["panel_on_axis"],
        n_training_genes_dropped_from_axis=len(axes["train_on_axis"]),
        common4_canonical=int(len(axes["common4"])),
        ridge_rows=receipt["n_rows"],ridge_features=receipt["n_features"],
        ridge_lambda_rule="core341.fit: lambda = number of rows, on RMS-scaled features",
        ridge_lambda=int(receipt["n_rows"]),
        held_target_poison_pass=True,
        subset_linear_max_difference=receipt["subset_linear_max_difference"],
        normal_equation_max_residual=receipt["normal_equation_max_residual"],
        heldout_within_bank=ho,
        cosine_blended_vs_consensus=dict(median=float(np.median(cos)),min=float(cos.min()),max=float(cos.max()),mean=float(cos.mean())),
        cosine_raw_vs_consensus=dict(median=float(np.median(cos_raw)),min=float(cos_raw.min()),max=float(cos_raw.max())),
        covered_fraction_of_canonical=float(len(axes["covered"])/len(axes["canonical"])),
        bank_shas=bank_shas,
        started_utc=started,finished_utc=datetime.now(timezone.utc).isoformat(),
        seconds=time.monotonic()-clock,
        max_rss_bytes_macos=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        deviations=["assignment says 100 A targets; panel is 300 shared, 259 with full consensus, 0 in Jurkat bank — shipped all 259 in panel order",
            "no measured_mask for Jurkat; axis = 4-bank∩canonical − training − panel",
            "matrix = beta=0.5 fixed-norm blend (346 blend verbatim, design row = consensus row); raw/consensus/norm also stored"])
    (WORK/"ctx_a_calibrated.json").write_text(json.dumps(sidecar,indent=2,sort_keys=True)+"\n")
    (WORK/"fit372.json").write_text(json.dumps({k:v for k,v in sidecar.items() if k!="bank_shas"},indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:sidecar[k] for k in ("n_training_targets","n_covered_genes","n_query_targets","ridge_lambda","held_target_poison_pass")},indent=2),flush=True)

if __name__=="__main__": main()
