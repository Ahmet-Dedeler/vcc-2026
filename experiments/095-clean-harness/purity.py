"""Does over-dispersion put BETTER genes first, or just shorten the list?

`reach` is the deepest prefix of the reference's DE pool, ordered by our confidence, whose
sign purity clears 0.9. It can rise for two very different reasons: the head genuinely gets
purer, or the pool gets reordered so that a lucky run of correct signs sits at the front.
The first transfers; the second is noise. This measures the purity profile directly.
"""
import os, sys
import numpy as np, anndata as ad
from scipy.stats import false_discovery_control, norm
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src"); sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
from vcc2026.donor_transfer import emit_counts
from de_kernel import ctrl_tiesums, wilcoxon_z
REAL=f"{ROOT}/data/harness/orion-hct116-thinned.h5ad"
PRED=f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"
CTRL,COL,ALPHA,MIN_CPM,N,L,SP = "non-targeting","target_gene",0.05,5.0,400,100_000,0.25
def nl(x):
    s=x.sum(1); s[s==0]=1.0
    return np.log1p(x*(1e6/s[:,None]))
R=ad.read_h5ad(REAL); rl_=np.asarray(R.obs[COL].astype(str).values,dtype="U64")
X=R.X; X=np.asarray(X.todense() if hasattr(X,"todense") else X,dtype=np.float64)
ci=np.flatnonzero(rl_==CTRL); C=nl(X[ci]); cpm=np.expm1(C).mean(0); keep=cpm>MIN_CPM
cs=np.ascontiguousarray(np.sort(C[:,keep].T,axis=1)); tie=ctrl_tiesums(cs)
ref=np.expm1(C[:,keep]).mean(0); del C
G_full=X.shape[1]; gk=np.asarray([str(g) for g in R.var_names])[keep]
fr=np.empty(gk.size,np.int64); fr[np.argsort(gk)]=np.arange(gk.size)
gpos={g:i for i,g in enumerate(gk)}
P=ad.read_h5ad(PRED); pl_=np.asarray(P.obs[COL].astype(str).values,dtype="U64")
tg=[t for t in sorted(set(rl_)-{CTRL}) if (pl_==t).sum()>0][:int(os.environ.get("N_PERTS","30"))]
DEPTHS=[1,2,5,10,20,50,100]
prof={k:{d:[] for d in DEPTHS} for k in ("production","phi4")}
for i,t in enumerate(tg):
    rr=np.flatnonzero(rl_==t); B=nl(X[rr])
    z=wilcoxon_z(cs,tie,np.ascontiguousarray(B[:,keep].T))
    rq=false_discovery_control(2*norm.sf(np.abs(z)),method="bh")
    rlv_all=np.log2((np.expm1(B[:,keep]).mean(0)+1e-9)/(ref+1e-9)); del B
    rs=rq<ALPHA
    j=gpos.get(t)
    if j is not None: rs[j]=False
    gs=np.flatnonzero(rs)
    if gs.size<20: continue
    pr=P[np.flatnonzero(pl_==t)].to_memory().X
    tot=np.asarray((pr.todense() if hasattr(pr,"todense") else pr)).sum(0).astype(np.int64); del pr
    prof_p=(tot+0.5)/(N*L)
    for name,phi in (("production",0.0),("phi4",4.0)):
        rng=np.random.default_rng(90_000+i)
        v,c,ip=emit_counts(prof_p,np.full(N,L),0.1,rng,mode="quasi",
                           spread=np.full(G_full,SP),overdispersion=phi)
        blk=np.zeros((N,G_full))
        for r in range(N):
            sl=slice(ip[r],ip[r+1]); blk[r,c[sl]]=v[sl]
        nb=nl(blk); del blk
        z2=wilcoxon_z(cs,tie,np.ascontiguousarray(nb[:,keep].T))
        pp2=2*norm.sf(np.abs(z2)); pq2=false_discovery_control(pp2,method="bh")
        plv_all=np.log2((np.expm1(nb[:,keep]).mean(0)+1e-9)/(ref+1e-9)); del nb
        ps=pq2<ALPHA
        rlv,plv=rlv_all[gs],plv_all[gs]
        idn=(~np.isnan(rlv))&(rlv!=0)
        mt=idn&(~np.isnan(plv))&(plv!=0)&(np.sign(plv)==np.sign(rlv))
        absl=np.abs(np.nan_to_num(plv,nan=-np.inf))
        o=np.lexsort((fr[gs],-absl,pp2[gs],pq2[gs],~ps[gs]))
        nd=np.cumsum(idn[o]); nmt=np.cumsum(mt[o])
        for d in DEPTHS:
            k=np.searchsorted(nd,d)
            if k<nd.size and nd[k]>0: prof[name][d].append(nmt[k]/nd[k])
print(f"{len(tg)} perturbations\n")
print(f"{'depth into the reference DE pool':<34}{'production':>12}{'phi=4':>10}{'delta':>9}")
for d in DEPTHS:
    a=np.mean(prof['production'][d]); b=np.mean(prof['phi4'][d])
    print(f"  first {d:>3} adjudicable genes{'':<12}{a:>12.4f}{b:>10.4f}{b-a:>+9.4f}")
print("\n0.9 is the purity pass mark reach requires.")
