"""Does over-dispersion make the emitted p-order track our own predicted effect?

On the real challenge artefact, Spearman(p_adj rank, |lfc| rank) is +0.015 -- our confidence
order is essentially independent of our model's own magnitudes, because Poisson noise is too
weak to compete with each gene's control-distribution skewness. If over-dispersion works for
the reason section 9 claims, that correlation should RISE.
"""
import os, sys
import numpy as np, anndata as ad
from scipy.stats import false_discovery_control, norm, spearmanr
ROOT="~/Code/virtual-cell-challenge-2026"
sys.path.insert(0,f"{ROOT}/src"); sys.path.insert(0,f"{ROOT}/experiments/092-call-set-truncation")
from vcc2026.donor_transfer import emit_counts
from de_kernel import ctrl_tiesums, wilcoxon_z
REAL=f"{ROOT}/data/harness/orion-hct116-thinned.h5ad"
PRED=f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"
CTRL,COL,ALPHA,MIN_CPM,N,L,SP="non-targeting","target_gene",0.05,5.0,400,100_000,0.25
def nl(x):
    s=x.sum(1); s[s==0]=1.0
    return np.log1p(x*(1e6/s[:,None]))
R=ad.read_h5ad(REAL); lab=np.asarray(R.obs[COL].astype(str).values,dtype="U64")
X=R.X; X=np.asarray(X.todense() if hasattr(X,"todense") else X,dtype=np.float64)
ci=np.flatnonzero(lab==CTRL); C=nl(X[ci]); cpm=np.expm1(C).mean(0); keep=cpm>MIN_CPM
cs=np.ascontiguousarray(np.sort(C[:,keep].T,axis=1)); tie=ctrl_tiesums(cs)
ref=np.expm1(C[:,keep]).mean(0); del C
G_full=X.shape[1]
P=ad.read_h5ad(PRED); pl=np.asarray(P.obs[COL].astype(str).values,dtype="U64")
tg=[t for t in sorted(set(lab)-{CTRL}) if (pl==t).sum()>0][:int(os.environ.get("N_PERTS","25"))]
out={0.0:[],2.0:[],4.0:[],8.0:[]}
for i,t in enumerate(tg):
    pr=P[np.flatnonzero(pl==t)].to_memory().X
    tot=np.asarray((pr.todense() if hasattr(pr,"todense") else pr)).sum(0).astype(np.int64); del pr
    prof=(tot+0.5)/(N*L)
    base=None
    for phi in out:
        rng=np.random.default_rng(90_000+i)
        v,c,ip=emit_counts(prof,np.full(N,L),0.1,rng,mode="quasi",
                           spread=np.full(G_full,SP),overdispersion=phi)
        blk=np.zeros((N,G_full))
        for r in range(N):
            sl=slice(ip[r],ip[r+1]); blk[r,c[sl]]=v[sl]
        nb=nl(blk); del blk
        z=wilcoxon_z(cs,tie,np.ascontiguousarray(nb[:,keep].T))
        q=false_discovery_control(2*norm.sf(np.abs(z)),method="bh")
        lfc=np.log2((np.expm1(nb[:,keep]).mean(0)+1e-9)/(ref+1e-9)); del nb
        if phi==0.0: base=np.abs(lfc)      # |effect| from the PRODUCTION emission, fixed
        m=np.isfinite(q)&np.isfinite(base)
        out[phi].append(spearmanr(q[m], -base[m]).statistic)
print(f"{len(tg)} perturbations\n")
print(f"{'phi':>6}{'Spearman(p_adj, |predicted effect|)':>40}")
for phi in out:
    v=np.array(out[phi]); print(f"{phi:>6.1f}{np.nanmean(v):>40.4f}")
print("\nHigher = the emitted p-value order tracks our own model's magnitudes more closely.")
