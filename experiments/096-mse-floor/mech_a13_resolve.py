"""Re-solve the shrinkage operating point on the SHIPPED estimator, post tie-break fix.

Imports `agreement_shrink_vector`, `donor_cosine_key`, `pert_shrink_vector` and
`MultiDonorModel` from the builder at HEAD rather than reimplementing them -- the previous
solve used `argsort(argsort(-key))`, which ranked the ~16.5k zero-agreement genes by ARRAY
POSITION and is a different estimator from the one that will ship.

APPLICABILITY, checked before transferring anything. Production keys on THREE donors
(k562_gwps, orion_hct116, orion_hek293t). `orion_hct116` is the same cell line as the
orion-hct116-thinned harness, so that surface can only ever run a TWO-donor key -- it is
structurally incapable of exhibiting the shipped estimator. h1-test is H1 hESC and none of the
three donors is H1, so it can run the real thing.
"""
import importlib.util as ilu, sys
import numpy as np
ROOT="~/Code/virtual-cell-challenge-2026"
sys.path.insert(0,f"{ROOT}/src"); sys.path.insert(0,f"{ROOT}/experiments/096-mse-floor")
_s=ilu.spec_from_file_location("mkpm",f"{ROOT}/scripts/make_submission_pm_cd4.py")
MK=ilu.module_from_spec(_s); sys.modules['mkpm']=MK; _s.loader.exec_module(MK)
from vcc2026.multi_donor import MultiDonorModel
from vcc2026.donor_transfer import Donor
from cell_eval2.metrics.delta import _across_pert_budget, _row_weights
from cell_eval2.prep import bulk_lognorm_means
from decompose import terms
from msekit import TS, group_stats

PROD=MK.PRODUCTION_DONORS
print(f"shipped: PRODUCTION_DONORS={PROD}  AGREEMENT_SHRINK_RANK_SCALE={MK.AGREEMENT_SHRINK_RANK_SCALE}")

def model_for(genes, names):
    return MultiDonorModel({n: Donor.load(f"{ROOT}/data/donors/{n}") for n in names},
                           np.asarray(genes), gene_norm=False, eb=False)

def agree_rows(genes, targets, names):
    m=model_for(genes,names)
    return m, np.vstack([m.donor_agreement(str(t)) for t in targets])

def alphas(A, targets, gpos, gamma):
    out=np.empty_like(A)
    for i,t in enumerate(targets):
        out[i]=MK.agreement_shrink_vector(A[i], gamma, protect=gpos.get(str(t)))
    return out

GAM=[0.0,0.25,0.5,0.75,1.0,1.5,2.0,3.0]
# ---------------- applicability ----------------
z=np.load(f"{ROOT}/experiments/096-mse-floor/work/a10-stats.npz",allow_pickle=True)
genes=z['genes'].astype(str); G=genes.size; gpos={g:i for i,g in enumerate(genes)}
perts=z['perts'].astype(str)
selA=np.flatnonzero(np.char.startswith(perts,"A|")); tgtA=[p.split("|",1)[1] for p in perts[selA]]
_,KA=agree_rows(genes,tgtA,PROD)
print(f"\nAPPLICABILITY -- non-zero genes per target in donor_agreement")
print(f"  challenge (3 donors)      median {int(np.median((KA>0).sum(1))):6d}  "
      f"rows with <2 distinct {int((np.array([np.unique(r).size for r in KA])<2).sum())}/{len(tgtA)}")
for nm,path,names in (("h1-test  (3 donors, leak-free)",f"{ROOT}/data/harness/h1-test.h5ad",PROD),
                      ("orion    (2 donors, forced)",f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
                       tuple(n for n in PROD if n!="orion_hct116"))):
    rp,_,_,_,rg=group_stats(path)
    tg=[str(p) for p in rp if str(p)!="non-targeting"]
    _,Kh=agree_rows(rg,tg,names)
    print(f"  {nm:26s} median {int(np.median((Kh>0).sum(1))):6d}  "
          f"rows with <2 distinct {int((np.array([np.unique(r).size for r in Kh])<2).sum())}/{len(tg)}")

# ---------------- challenge: A'(g), B'(g) exact ----------------
ctrl=np.load(f"{ROOT}/experiments/096-mse-floor/work/ctrl-stats.npz",allow_pickle=True)
Ag={g:[] for g in GAM}; Bg={g:[] for g in GAM}; SV={g:[] for g in GAM}
for c in "ABC":
    sel=np.flatnonzero(np.char.startswith(perts,f"{c}|"))
    tgt=[p.split("|",1)[1] for p in perts[sel]]
    bctrl=bulk_lognorm_means(ctrl[f'sum_{c}'][None,:].astype(np.float64),TS)[0]
    cols={k:gpos[t] for k,t in enumerate(tgt) if t in gpos}
    w=_row_weights(len(tgt),G,cols)
    dh=bulk_lognorm_means(z['sums'][sel].astype(np.float64),TS)-bctrl[None,:]
    for k,col in cols.items(): dh[k,col]=0.0
    _,KK=agree_rows(genes,tgt,PROD)
    for g in GAM:
        a=alphas(KK,tgt,gpos,g); D=a*dh
        Ag[g].append(float((np.einsum("ij,ij->i",D,D)*w).sum()))
        Bg[g].append(_across_pert_budget(bctrl[None,:]+D,cols,w))
        SV[g].append(float(np.mean(np.linalg.norm(D,axis=1)/np.maximum(np.linalg.norm(dh,axis=1),1e-12))))
A={g:float(np.mean(Ag[g])) for g in GAM}; B={g:float(np.mean(Bg[g])) for g in GAM}
S={g:float(np.mean(SV[g])) for g in GAM}

# ---------------- harness: X'/X ----------------
mult={g:{} for g in GAM}
for nm,real,pred,names in (("h1",f"{ROOT}/data/harness/h1-test.h5ad",
   f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",PROD),
  ("orion",f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
   f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
   tuple(n for n in PROD if n!="orion_hct116"))):
    rp,rs,rc,rj,rg=group_stats(real); pp,ps,pc,pj,pg=group_stats(pred)
    t=terms(pp,ps,pc,pj,rp,rs,rc,rj,rg); dh,d,w=t["dh"],t["d"],t["w"]
    for k,c2 in t["cols"].items(): dh[k,c2]=0.0; d[k,c2]=0.0
    tg=[str(x) for x in pp[t["keep"]]]; hp={g2:i for i,g2 in enumerate(rg)}
    _,KK=agree_rows(rg,tg,names)
    X0=float((np.einsum("ij,ij->i",dh,d)*w).sum())
    for g in GAM:
        a=alphas(KK,tg,hp,g)
        mult[g][nm]=float((np.einsum("ij,ij->i",a*dh,d)*w).sum())/X0
np.savez(f"{ROOT}/experiments/096-mse-floor/work/a13_resolve.npz",
         gam=np.array(GAM), A=np.array([A[g] for g in GAM]), B=np.array([B[g] for g in GAM]),
         surv=np.array([S[g] for g in GAM]),
         mh1=np.array([mult[g]['h1'] for g in GAM]), mor=np.array([mult[g]['orion'] for g in GAM]))
bs=np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve2.npz")
X0b=float(bs['X']); sD=float(bs['sumD']); Kc=float(bs['Kc']); K=float(bs['K'])
BA,SP=0.987862879,0.950162566
print(f"\n{'gamma':>6}{'surv':>8}{'A''':>10}{'B''':>10}{'(A-B)/A':>9}{'X/X h1':>9}{'orion':>8}"
      f"{'s*=escale':>11}{'raw mse':>10}{'score':>8}{'overall':>9}")
for g in GAM:
    m=mult[g]['h1']                      # h1-test is the only applicable surface (see above)
    Xp=X0b*m; ss=np.linspace(0.05,8.0,6000)
    v=1+(ss**2*A[g]-2*ss*Xp+Kc-np.minimum(K,ss**2*B[g]))/sD
    i=int(np.argmin(v)); sc=min(max((BA-v[i])/SP,0),1)
    print(f"{g:>6.2f}{S[g]:>8.4f}{A[g]:>10.5f}{B[g]:>10.5f}{(A[g]-B[g])/A[g]:>9.4f}"
          f"{mult[g]['h1']:>9.3f}{mult[g]['orion']:>8.3f}{ss[i]:>11.4f}{v[i]:>10.4f}{sc:>8.4f}{sc/6:>+9.5f}")
