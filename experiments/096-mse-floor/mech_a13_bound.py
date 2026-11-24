"""What bounds gamma, now that the mse surface is monotone in it?

The corrected estimator makes raw mse fall monotonically out to gamma=3 -- an unbounded gain is
a sign of extrapolating past where the model is valid, not of a free lunch. Three candidate
bounds, all measurable: profile realisability (clipping), the call set, and pds-side direction.
Also: does the per-perturbation multiplier COMPOSE with the per-gene one, or duplicate it?
"""
import importlib.util as ilu, sys
import numpy as np, anndata as ad
from scipy.sparse import csr_matrix
from scipy.stats import false_discovery_control, norm as _n
ROOT="~/Code/virtual-cell-challenge-2026"
sys.path.insert(0,f"{ROOT}/src"); sys.path.insert(0,f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0,f"{ROOT}/experiments/092-call-set-truncation")
_s=ilu.spec_from_file_location("mkpm",f"{ROOT}/scripts/make_submission_pm_cd4.py")
MK=ilu.module_from_spec(_s); sys.modules['mkpm']=MK; _s.loader.exec_module(MK)
from vcc2026.multi_donor import MultiDonorModel
from vcc2026.donor_transfer import Donor, emit_counts
from de_kernel import ctrl_tiesums, wilcoxon_z
from cell_eval2.metrics.delta import _across_pert_budget, _row_weights
from cell_eval2.prep import bulk_lognorm_means
from decompose import terms
from msekit import TS, group_stats
PROD=MK.PRODUCTION_DONORS
def nl(x):
    s=x.sum(1); s[s==0]=1.0
    return np.log1p(x*(1e6/s[:,None]))

z=np.load(f"{ROOT}/experiments/096-mse-floor/work/a10-stats.npz",allow_pickle=True)
genes=z['genes'].astype(str); G=genes.size; gpos={g:i for i,g in enumerate(genes)}
perts=z['perts'].astype(str)
ctrl=np.load(f"{ROOT}/experiments/096-mse-floor/work/ctrl-stats.npz",allow_pickle=True)
ctx="A"; sel=np.flatnonzero(np.char.startswith(perts,f"{ctx}|"))
tgt=[p.split("|",1)[1] for p in perts[sel]]
bctrl=bulk_lognorm_means(ctrl[f'sum_{ctx}'][None,:].astype(np.float64),TS)[0]
dh=bulk_lognorm_means(z['sums'][sel].astype(np.float64),TS)-bctrl[None,:]
mdl=MultiDonorModel({n:Donor.load(f"{ROOT}/data/donors/{n}") for n in PROD},genes,gene_norm=False,eb=False)
KK=np.vstack([mdl.donor_agreement(str(t)) for t in tgt])
zs,carried=MK.donor_cosine_key(mdl,[str(t) for t in tgt])
print(f"per-pert key: carried {int(carried.sum())}/{len(tgt)}")

A=ad.read_h5ad(f"{ROOT}/data/challenge/context_{ctx}.h5ad")
C=nl(np.asarray(A.X.todense(),dtype=np.float64)); cpm=np.expm1(C).mean(0); keep=cpm>5.0
cs=np.ascontiguousarray(np.sort(C[:,keep].T,axis=1)); tie=ctrl_tiesums(cs); del C,A
idx=np.linspace(0,len(tgt)-1,12).astype(int)
GS=[(0.0,0.9366,0.0),(0.5,1.5674,0.0),(1.0,1.7397,0.0),(2.0,1.9822,0.0),(1.0,None,1.0)]
print(f"\n{'gamma_g':>8}{'g_pert':>8}{'escale':>9}{'clip%':>8}{'n_pred':>9}{'|eff| p99 x prod':>18}")
for gg,esc,gp in GS:
    if esc is None:  # composed arm: recompute s* with both weights in
        w_=_row_weights(len(tgt),G,{k:gpos[t] for k,t in enumerate(tgt) if t in gpos})
        zp=MK.pert_shrink_vector(zs,carried,gp)
        al=np.vstack([MK.agreement_shrink_vector(KK[i],gg,protect=gpos.get(str(t))) for i,t in enumerate(tgt)])
        D=(zp[:,None]*al)*dh
        cols={k:gpos[t] for k,t in enumerate(tgt) if t in gpos}
        Bp=_across_pert_budget(bctrl[None,:]+D,cols,w_)
        esc=float(np.sqrt(0.246529/Bp)); print(f"   [composed: recomputed s* = {esc:.4f}]")
    clip=[]; npred=[]; amp=[]
    for j in idx:
        al=MK.agreement_shrink_vector(KK[j],gg,protect=gpos.get(str(tgt[j])))
        zp=MK.pert_shrink_vector(zs,carried,gp)[j] if gp>0 else 1.0
        eff=esc*zp*al*dh[j]
        raw=np.expm1(bctrl+eff)
        clip.append(float((raw<0).mean())); amp.append(float(np.percentile(np.abs(eff),99)/max(np.percentile(np.abs(dh[j]),99),1e-12)))
        prof=np.clip(raw,0,None); prof/=prof.sum()
        top=np.argpartition(-np.abs(eff),8000)[:8000]
        sp=np.full(G,1.0); ph=np.full(G,4.0); sp[top]=0.25
        v,col,ip=emit_counts(prof,np.full(400,z['sums'][sel][j].sum()/400),0.0,
                             np.random.default_rng(11000+j),mode="quasi",spread=sp,overdispersion=ph)
        M=csr_matrix((v.astype(np.float64),col.astype(np.int64),ip.astype(np.int64)),shape=(400,G))
        P=nl(np.asarray(M.todense(),dtype=np.float64))
        zz=wilcoxon_z(cs,tie,np.ascontiguousarray(P[:,keep].T))
        padj=false_discovery_control(np.nan_to_num(2.0*_n.sf(np.abs(zz)),nan=1.0),method="bh")
        npred.append(int((padj<0.05).sum()))
    print(f"{gg:>8.2f}{gp:>8.2f}{esc:>9.4f}{100*np.mean(clip):>8.3f}{np.mean(npred):>9.0f}{np.mean(amp):>18.2f}")
