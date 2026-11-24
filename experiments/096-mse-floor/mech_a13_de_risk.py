"""Does the shrinkage move amplitude INTO or OUT OF the genes we actually call?

The mse gain is settled. The open risk is the DE members: shrinkage multiplies each gene by
alpha = exp(-gamma*donor_rank/2000) and then effect_scale multiplies everything back up. A called
gene whose donor rank is high gets alpha*scale < 1 and may drop out; one at rank 0 gets
alpha*scale > 1 and is reinforced. Which dominates is a property of where our CALLED genes sit in
the donor-rank order, and that is measurable against the real challenge controls.
"""
import sys, time
import numpy as np, anndata as ad
from scipy.sparse import csr_matrix
from scipy.stats import false_discovery_control, norm as _n
ROOT="~/Code/virtual-cell-challenge-2026"
sys.path.insert(0,f"{ROOT}/src"); sys.path.insert(0,f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0,f"{ROOT}/experiments/092-call-set-truncation")
from de_kernel import ctrl_tiesums, wilcoxon_z
from cell_eval2.prep import bulk_lognorm_means
from msekit import TS
from mech_miller import donor_key
from vcc2026.donor_transfer import emit_counts
MIN_CPM,ALPHA,NC=5.0,0.05,400
def nl(x):
    s=x.sum(1); s[s==0]=1.0
    return np.log1p(x*(1e6/s[:,None]))
ctx="A"; NP=12
A=ad.read_h5ad(f"{ROOT}/data/challenge/context_{ctx}.h5ad")
Xc=np.asarray(A.X.todense(),dtype=np.float64); C=nl(Xc)
cpm=np.expm1(C).mean(0); keep=cpm>MIN_CPM
cs=np.ascontiguousarray(np.sort(C[:,keep].T,axis=1)); tie=ctrl_tiesums(cs); del C,Xc
z=np.load(f"{ROOT}/experiments/096-mse-floor/work/a10-stats.npz",allow_pickle=True)
genes=z['genes'].astype(str); G=genes.size
ctrl=np.load(f"{ROOT}/experiments/096-mse-floor/work/ctrl-stats.npz",allow_pickle=True)
bctrl=bulk_lognorm_means(ctrl[f'sum_{ctx}'][None,:].astype(np.float64),TS)[0]
perts=z['perts'].astype(str); sel=np.flatnonzero(np.char.startswith(perts,f"{ctx}|"))
tgt=np.asarray([p.split("|",1)[1] for p in perts[sel]])
pm=bulk_lognorm_means(z['sums'][sel].astype(np.float64),TS); DH=pm-bctrl[None,:]
idx=np.linspace(0,sel.size-1,NP).astype(int)
dk=donor_key(genes,list(tgt[idx])); blank=(dk.sum(1)==0)
rank=np.argsort(np.argsort(-dk,axis=1),axis=1)
print(f"context {ctx}: {int(keep.sum())}/{G} genes past the DE gate, {NP} perturbations, "
      f"{int(blank.sum())} blank donor rows")
for gam,scale in ((0.0,0.9368),(0.5,1.4025),(1.0,1.5423)):
    npred=[]; net_called=[]; frac_up=[]
    for j,i in enumerate(idx):
        e=DH[i].copy()
        a=np.exp(-gam*rank[j]/2000.0)
        if blank[j]: a=np.ones(G)
        eff=scale*a*e
        prof=np.expm1(bctrl+eff); np.clip(prof,0,None,out=prof); prof/=prof.sum()
        top=np.argpartition(-np.abs(eff),8000)[:8000]
        sp=np.full(G,1.0); ph=np.full(G,4.0); sp[top]=0.25; ph[top]=4.0
        libs=np.full(NC,z['sums'][sel][i].sum()/NC)
        v,col,ip=emit_counts(prof,libs,0.0,np.random.default_rng(9000+i),mode="quasi",
                             spread=sp,overdispersion=ph)
        M=csr_matrix((v.astype(np.float64),col.astype(np.int64),ip.astype(np.int64)),shape=(NC,G))
        P=nl(np.asarray(M.todense(),dtype=np.float64))
        zz=wilcoxon_z(cs,tie,np.ascontiguousarray(P[:,keep].T))
        padj=false_discovery_control(np.nan_to_num(2.0*_n.sf(np.abs(zz)),nan=1.0),method="bh")
        sig=padj<ALPHA; npred.append(int(sig.sum()))
        gnet=(scale*a)[keep][sig]
        net_called.append(float(np.median(gnet))); frac_up.append(float((gnet>=1.0).mean()))
    print(f"  gamma={gam:<4} effect_scale={scale:<7.4f} n_pred {np.mean(npred):7.0f}   "
          f"median (scale*alpha) on CALLED genes {np.mean(net_called):.3f}   "
          f"fraction of called genes reinforced {np.mean(frac_up):.3f}")
