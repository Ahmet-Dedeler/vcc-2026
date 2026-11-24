"""Is the pds disagreement about the CELL LINE or about the DONOR COUNT?

h1-test (3-donor key) says shrinkage HELPS pds at low gamma; orion (forced 2-donor) says it
hurts at every gamma. Those are different estimators, so the comparison is confounded. Run
h1-test with the SAME 2-donor set orion is forced to use. If h1-test flips, the driver is the
donor count and the challenge -- which ships 3 donors -- follows the 3-donor result. If h1-test
still helps, it is a genuine cell-line disagreement and neither surface settles it.
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
PROD=MK.PRODUCTION_DONORS; TWO=tuple(n for n in PROD if n!="orion_hct116")
B_PDS,R_PDS=0.50020,0.95228
def pds(dhat,dreal):
    U=dhat/np.maximum(np.linalg.norm(dhat,axis=1,keepdims=True),1e-300)
    V=dreal/np.maximum(np.linalg.norm(dreal,axis=1,keepdims=True),1e-300)
    Dm=1.0-U@V.T; n=Dm.shape[1]
    return float(np.mean([1.0-((Dm[p]<Dm[p,p]).sum()+0.5*((Dm[p]==Dm[p,p]).sum()-1))/(n-1)
                          for p in range(Dm.shape[0])]))
rp,rs,rc,rj,rg=group_stats(f"{ROOT}/data/harness/h1-test.h5ad")
pp,ps,pc,pj,pg=group_stats(f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad")
t=terms(pp,ps,pc,pj,rp,rs,rc,rj,rg); dh,d,w=t["dh"],t["d"],t["w"]
tg=[str(x) for x in pp[t["keep"]]]; hp={g:i for i,g in enumerate(rg)}
panel=np.array([hp[x] for x in tg if x in hp]); mask=np.ones(rg.size,bool); mask[panel]=False
X0=float((np.einsum("ij,ij->i",dh,d)*w).sum())
for names,lbl in ((PROD,"h1-test, 3 donors (= the shipped key)"),
                  (TWO ,"h1-test, 2 donors (= orion's forced key)")):
    mdl=MultiDonorModel({n:Donor.load(f"{ROOT}/data/donors/{n}") for n in names},rg,gene_norm=False,eb=False)
    KK=np.vstack([mdl.donor_agreement(x) for x in tg])
    print(f"\n{lbl}: median non-zero genes/target {int(np.median((KK>0).sum(1)))}")
    print(f"  {'gamma':>7}{'raw pds':>10}{'d scaled':>10}{'X''/X':>9}")
    base=None
    for g in (0.0,0.25,0.5,1.0,2.0,3.0):
        al=np.vstack([MK.agreement_shrink_vector(KK[i],g,protect=hp.get(tg[i])) for i in range(len(tg))])
        v=pds((al*dh)[:,mask],d[:,mask]); sc=(v-B_PDS)/(R_PDS-B_PDS)
        if base is None: base=sc
        m=float((np.einsum("ij,ij->i",al*dh,d)*w).sum())/X0
        print(f"  {g:>7.2f}{v:>10.4f}{sc-base:>+10.4f}{m:>9.3f}")
