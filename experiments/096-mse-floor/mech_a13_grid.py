"""Joint (gamma_gene, gamma_pert) solve: exact mse on the challenge + pds on the applicable surface."""
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
bs=np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve2.npz")
X0b=float(bs['X']); sD=float(bs['sumD']); Kc=float(bs['Kc']); K=float(bs['K'])
BA,SP,B_PDS,R_PDS=0.987862879,0.950162566,0.50020,0.95228
GG=[0.0,0.125,0.25,0.5,1.0]; GP=[0.0,0.5,1.0,2.0]

# challenge A'(gg,gp), B'(gg,gp)
z=np.load(f"{ROOT}/experiments/096-mse-floor/work/a10-stats.npz",allow_pickle=True)
genes=z['genes'].astype(str); G=genes.size; gpos={g:i for i,g in enumerate(genes)}
perts=z['perts'].astype(str)
ctrl=np.load(f"{ROOT}/experiments/096-mse-floor/work/ctrl-stats.npz",allow_pickle=True)
CH={}
for c in "ABC":
    sel=np.flatnonzero(np.char.startswith(perts,f"{c}|")); tgt=[p.split("|",1)[1] for p in perts[sel]]
    bctrl=bulk_lognorm_means(ctrl[f'sum_{c}'][None,:].astype(np.float64),TS)[0]
    cols={k:gpos[t] for k,t in enumerate(tgt) if t in gpos}; w=_row_weights(len(tgt),G,cols)
    dh=bulk_lognorm_means(z['sums'][sel].astype(np.float64),TS)-bctrl[None,:]
    for k,col in cols.items(): dh[k,col]=0.0
    mdl=MultiDonorModel({n:Donor.load(f"{ROOT}/data/donors/{n}") for n in PROD},genes,gene_norm=False,eb=False)
    KK=np.vstack([mdl.donor_agreement(str(t)) for t in tgt])
    zs,car=MK.donor_cosine_key(mdl,[str(t) for t in tgt])
    CH[c]=(dh,cols,w,bctrl,KK,zs,car,tgt)
# h1-test multipliers + pds
rp,rs,rc,rj,rg=group_stats(f"{ROOT}/data/harness/h1-test.h5ad")
pp,ps,pc,pj,pg=group_stats(f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad")
t=terms(pp,ps,pc,pj,rp,rs,rc,rj,rg); dhh,dd,wh=t["dh"],t["d"],t["w"]
tgh=[str(x) for x in pp[t["keep"]]]; hp={g:i for i,g in enumerate(rg)}
pan=np.array([hp[x] for x in tgh if x in hp]); msk=np.ones(rg.size,bool); msk[pan]=False
mh=MultiDonorModel({n:Donor.load(f"{ROOT}/data/donors/{n}") for n in PROD},rg,gene_norm=False,eb=False)
KH=np.vstack([mh.donor_agreement(x) for x in tgh]); zsh,carh=MK.donor_cosine_key(mh,tgh)
X0h=float((np.einsum("ij,ij->i",dhh,dd)*wh).sum())
def pds_of(D):
    U=D[:,msk]/np.maximum(np.linalg.norm(D[:,msk],axis=1,keepdims=True),1e-300)
    V=dd[:,msk]/np.maximum(np.linalg.norm(dd[:,msk],axis=1,keepdims=True),1e-300)
    Dm=1.0-U@V.T; n=Dm.shape[1]
    return float(np.mean([1.0-((Dm[p]<Dm[p,p]).sum()+0.5*((Dm[p]==Dm[p,p]).sum()-1))/(n-1) for p in range(Dm.shape[0])]))
base_pds=(pds_of(dhh)-B_PDS)/(R_PDS-B_PDS)
print(f"{'g_gene':>7}{'g_pert':>7}{'escale':>9}{'raw mse':>9}{'d_mse':>9}{'d_pds':>9}{'SUM':>9}")
best=None
for gg in GG:
    for gp in GP:
        As,Bs=[],[]
        for c in "ABC":
            dh,cols,w,bctrl,KK,zs,car,tgt=CH[c]
            al=np.vstack([MK.agreement_shrink_vector(KK[i],gg,protect=gpos.get(str(tgt[i]))) for i in range(len(tgt))])
            zp=MK.pert_shrink_vector(zs,car,gp) if gp>0 else np.ones(len(tgt))
            D=(zp[:,None]*al)*dh
            As.append(float((np.einsum("ij,ij->i",D,D)*w).sum())); Bs.append(_across_pert_budget(bctrl[None,:]+D,cols,w))
        Ax,Bx=float(np.mean(As)),float(np.mean(Bs))
        alh=np.vstack([MK.agreement_shrink_vector(KH[i],gg,protect=hp.get(tgh[i])) for i in range(len(tgh))])
        zph=MK.pert_shrink_vector(zsh,carh,gp) if gp>0 else np.ones(len(tgh))
        Dh=(zph[:,None]*alh)*dhh
        m=float((np.einsum("ij,ij->i",Dh,dd)*wh).sum())/X0h
        ss=np.linspace(0.05,9.0,7000)
        v=1+(ss**2*Ax-2*ss*X0b*m+Kc-np.minimum(K,ss**2*Bx))/sD
        i=int(np.argmin(v)); sc=min(max((BA-v[i])/SP,0),1)
        dmse=sc/6-0.1049/6
        dpds=((pds_of(Dh)-B_PDS)/(R_PDS-B_PDS)-base_pds)/6
        tot=dmse+dpds
        print(f"{gg:>7.3f}{gp:>7.2f}{ss[i]:>9.4f}{v[i]:>9.4f}{dmse:>+9.5f}{dpds:>+9.5f}{tot:>+9.5f}")
        if best is None or tot>best[0]: best=(tot,gg,gp,ss[i],v[i],dmse,dpds)
print(f"\nBEST: gamma_gene={best[1]}, gamma_pert={best[2]}, effect_scale={best[3]:.4f}, "
      f"raw mse {best[4]:.4f}, d_mse {best[5]:+.5f}, d_pds {best[6]:+.5f}, sum {best[0]:+.5f}")
