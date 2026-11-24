"""Why did `jac` fall under over-dispersion when both harnesses said it would rise?

a8 delivered jac -0.003439 scaled against a predicted +0.0022 -- most of the shortfall.
`jac = |R n R^| / |R u R^|`, so a smaller call set shrinks the union (good) and can shrink
the intersection (bad). Both harnesses showed jac RISING with phi, so whatever happened on
the board is something they do not reproduce. This decomposes the ratio into its parts.
"""
import os, sys
import numpy as np, anndata as ad
from scipy.stats import false_discovery_control, norm
ROOT="~/Code/virtual-cell-challenge-2026"
sys.path.insert(0,f"{ROOT}/src"); sys.path.insert(0,f"{ROOT}/experiments/092-call-set-truncation")
from vcc2026.donor_transfer import emit_counts
from de_kernel import ctrl_tiesums, wilcoxon_z
REAL=os.environ.get("REAL",f"{ROOT}/data/harness/orion-hct116-thinned.h5ad")
PRED=os.environ.get("PRED",f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad")
if not os.path.isabs(REAL): REAL=os.path.join(ROOT,REAL)
if not os.path.isabs(PRED): PRED=os.path.join(ROOT,PRED)
CTRL,COL,ALPHA,MIN_CPM,N,L,SP="non-targeting","target_gene",0.05,5.0,400,100_000,0.25
def nl(x):
    s=x.sum(1); s[s==0]=1.0
    return np.log1p(x*(1e6/s[:,None]))
R=ad.read_h5ad(REAL); lab=np.asarray(R.obs[COL].astype(str).values,dtype="U64")
X=R.X; X=np.asarray(X.todense() if hasattr(X,"todense") else X,dtype=np.float64)
ci=np.flatnonzero(lab==CTRL); C=nl(X[ci]); cpm=np.expm1(C).mean(0); keep=cpm>MIN_CPM
cs=np.ascontiguousarray(np.sort(C[:,keep].T,axis=1)); tie=ctrl_tiesums(cs)
ref=np.expm1(C[:,keep]).mean(0); del C
G_full=X.shape[1]; gk=np.asarray([str(g) for g in R.var_names])[keep]
gpos={g:i for i,g in enumerate(gk)}
P=ad.read_h5ad(PRED); pl=np.asarray(P.obs[COL].astype(str).values,dtype="U64")
tg=[t for t in sorted(set(lab)-{CTRL}) if (pl==t).sum()>0][:int(os.environ.get("N_PERTS","30"))]
res={}
for phi in (0.0,4.0):
    inter,union,npred,nreal,jac=[],[],[],[],[]
    for i,t in enumerate(tg):
        rr=np.flatnonzero(lab==t); B=nl(X[rr])
        z=wilcoxon_z(cs,tie,np.ascontiguousarray(B[:,keep].T)); del B
        rq=false_discovery_control(2*norm.sf(np.abs(z)),method="bh")
        rs=rq<ALPHA
        j=gpos.get(t)
        if j is not None: rs[j]=False
        pr=P[np.flatnonzero(pl==t)].to_memory().X
        tot=np.asarray((pr.todense() if hasattr(pr,"todense") else pr)).sum(0).astype(np.int64); del pr
        prof=(tot+0.5)/(N*L)
        rng=np.random.default_rng(90_000+i)
        v,c,ip=emit_counts(prof,np.full(N,L),0.1,rng,mode="quasi",
                           spread=np.full(G_full,SP),overdispersion=phi)
        blk=np.zeros((N,G_full))
        for r in range(N):
            sl=slice(ip[r],ip[r+1]); blk[r,c[sl]]=v[sl]
        nb=nl(blk); del blk
        z2=wilcoxon_z(cs,tie,np.ascontiguousarray(nb[:,keep].T)); del nb
        pq2=false_discovery_control(2*norm.sf(np.abs(z2)),method="bh")
        ps=pq2<ALPHA
        if j is not None: ps[j]=False
        I=int((rs&ps).sum()); U=int((rs|ps).sum())
        inter.append(I); union.append(U); npred.append(int(ps.sum())); nreal.append(int(rs.sum()))
        jac.append(1.0 if U==0 else I/U)
    res[phi]=(np.mean(inter),np.mean(union),np.mean(npred),np.mean(nreal),np.mean(jac))
print(f"{len(tg)} perturbations, harness {os.path.basename(REAL)}\n")
print(f"{'phi':>5}{'n_real':>9}{'n_pred':>9}{'|R n Rhat|':>12}{'|R u Rhat|':>12}{'jac':>9}{'recall':>9}{'precision':>11}")
for phi,(I,U,NP_,NR,J) in res.items():
    print(f"{phi:>5.1f}{NR:>9.0f}{NP_:>9.0f}{I:>12.1f}{U:>12.1f}{J:>9.4f}{I/max(NR,1):>9.3f}{I/max(NP_,1):>11.3f}")
a,b=res[0.0],res[4.0]
print(f"\n  intersection {a[0]:.1f} -> {b[0]:.1f}  ({(b[0]-a[0])/a[0]*100:+.1f}%)")
print(f"  union        {a[1]:.1f} -> {b[1]:.1f}  ({(b[1]-a[1])/a[1]*100:+.1f}%)")
print("\n  jac rises only if the union shrinks faster than the intersection.")
