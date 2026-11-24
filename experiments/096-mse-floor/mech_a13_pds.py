"""pds_cosine under the shipped shrinkage. The one member the mse solve cannot see.

pds ranks each predicted effect against every measured effect by cosine, with ALL 300 panel
target genes removed from the feature space (`exclusion_scope: panel`) and midrank ties.
Scale-invariant, so effect_scale cannot touch it -- but the per-gene reweighting CHANGES
DIRECTION, and that is exactly what pds reads. Computed exactly on both harnesses.
"""
import importlib.util as ilu, sys
import numpy as np
ROOT="~/Code/virtual-cell-challenge-2026"
sys.path.insert(0,f"{ROOT}/src"); sys.path.insert(0,f"{ROOT}/experiments/096-mse-floor")
_s=ilu.spec_from_file_location("mkpm",f"{ROOT}/scripts/make_submission_pm_cd4.py")
MK=ilu.module_from_spec(_s); sys.modules['mkpm']=MK; _s.loader.exec_module(MK)
from vcc2026.multi_donor import MultiDonorModel
from vcc2026.donor_transfer import Donor
from decompose import terms
from msekit import group_stats
PROD=MK.PRODUCTION_DONORS
B_PDS,R_PDS=0.50020,0.95228     # board anchors

def pds(dhat,dreal):
    U=dhat/np.maximum(np.linalg.norm(dhat,axis=1,keepdims=True),1e-300)
    V=dreal/np.maximum(np.linalg.norm(dreal,axis=1,keepdims=True),1e-300)
    Dm=1.0-U@V.T; n=Dm.shape[1]; out=[]
    for p in range(Dm.shape[0]):
        own=Dm[p,p]
        k=int((Dm[p]<own).sum())+0.5*(int((Dm[p]==own).sum())-1)
        out.append(1.0-k/(n-1))
    return float(np.mean(out))

for nm,real,pred,names in (
  ("h1-test (3 donors, applicable)",f"{ROOT}/data/harness/h1-test.h5ad",
   f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",PROD),
  ("orion (2 donors, NOT applicable)",f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
   f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
   tuple(n for n in PROD if n!="orion_hct116"))):
    rp,rs,rc,rj,rg=group_stats(real); pp,ps,pc,pj,pg=group_stats(pred)
    t=terms(pp,ps,pc,pj,rp,rs,rc,rj,rg); dh,d=t["dh"],t["d"]
    tg=[str(x) for x in pp[t["keep"]]]; hp={g:i for i,g in enumerate(rg)}
    panel=np.array([hp[x] for x in tg if x in hp])
    mask=np.ones(rg.size,bool); mask[panel]=False        # exclusion_scope: panel
    mdl=MultiDonorModel({n:Donor.load(f"{ROOT}/data/donors/{n}") for n in names},rg,gene_norm=False,eb=False)
    KK=np.vstack([mdl.donor_agreement(x) for x in tg])
    zs,carried=MK.donor_cosine_key(mdl,tg)
    print(f"\n{nm}: {len(tg)} perturbations")
    print(f"  {'gamma_g':>8}{'g_pert':>8}{'raw pds':>10}{'scaled':>9}{'d scaled':>10}")
    base=None
    for gg,gp in ((0.0,0.0),(0.25,0.0),(0.5,0.0),(1.0,0.0),(2.0,0.0),(3.0,0.0),(0.0,1.0),(1.0,1.0)):
        al=np.vstack([MK.agreement_shrink_vector(KK[i],gg,protect=hp.get(tg[i])) for i in range(len(tg))])
        zp=MK.pert_shrink_vector(zs,carried,gp) if gp>0 else np.ones(len(tg))
        v=pds((zp[:,None]*al*dh)[:,mask], d[:,mask])
        sc=(v-B_PDS)/(R_PDS-B_PDS)
        if base is None: base=sc
        print(f"  {gg:>8.2f}{gp:>8.2f}{v:>10.4f}{sc:>9.4f}{sc-base:>+10.4f}")
