"""174b -- condition the recipient co-expression on cell state (remove top PCs) and re-test."""
import numpy as np, sys, time, h5py
sys.path.insert(0, "experiments/153-de-bench"); sys.path.insert(0, "experiments/168-sharpening")
from debench import Surface
from h1_tiebreak import labels
t0=time.time()
S = Surface("experiments/164-callset-artefact/de164.a15.npz")
E = np.load("experiments/169-sharpening-six-member/work/h1test-effect.npz", allow_pickle=True)
F = np.load("experiments/172-sign-confidence-key/work/h1test-feats.npz", allow_pickle=True)
eg = {str(g): i for i, g in enumerate(E["genes"])}; et = {str(t): i for i, t in enumerate(E["targets"])}
gi = np.array([eg[g] for g in S.genes]); pi = np.array([et[p] for p in S.perts])
H = E["H"][pi][:, gi].astype(float); cons = F["cons"].astype(float)
with h5py.File("data/harness/h1-test.h5ad","r") as h:
    genes=labels(h,"var/_index"); lab=labels(h,"obs/target_gene")
    ci=np.flatnonzero(lab=="non-targeting"); n,G=(int(x) for x in h["X"].attrs["shape"]); ip=h["X/indptr"][:]; rows=[]
    for r0 in range(0,ci.size,3000):
        idx=ci[r0:r0+3000]; blk=np.zeros((idx.size,G),dtype=np.float32)
        for k,r in enumerate(idx):
            s,e=int(ip[r]),int(ip[r+1]); blk[k,h["X/indices"][s:e]]=h["X/data"][s:e]
        tot=blk.sum(1,keepdims=True); tot[tot==0]=1; rows.append(np.log1p(blk*(1e6/tot)))
    X=np.concatenate(rows); del rows
X=X[:, gi]                       # restrict to the DE universe (10,779 genes) -- what matters, and smaller
X-=X.mean(0,keepdims=True); sd=X.std(0); sd[sd==0]=1; X/=sd
tpos={g:i for i,g in enumerate(S.genes)}; tidx=np.array([tpos[str(p)] for p in S.perts if str(p) in tpos])
pmask=np.array([str(p) in tpos for p in S.perts]); print(f"{pmask.sum()}/85 targets are in the DE universe")
pool=(S.real_sig&S.adj)[pmask]; sT=np.sign(S.real_lfc)[pmask]; Hm=H[pmask]; consm=cons[pmask]
print(f"data ready [{time.time()-t0:.0f}s]; computing PCs on a 6000-cell subsample")
rng=np.random.default_rng(0); sub=rng.choice(X.shape[0],6000,replace=False)
U,s_,Vt=np.linalg.svd(X[sub],full_matrices=False); print(f"svd done [{time.time()-t0:.0f}s]; top-10 var share {np.sum(s_[:10]**2)/np.sum(s_**2):.3f}")
def acc(key,m): return np.mean(key[m]==sT[m])
ac=np.abs(consm); med=np.median(ac[pool&(ac>0)]); weak=pool&(ac<med)
print(f"\n{'PCs removed':>12} {'-sign(r) acc':>13} {'on weak-donor half':>19} {'H acc|agree':>12} {'H acc|disagree':>15} {'|r| median':>11}")
for k in (0,2,5,10,20,40):
    Xr = X - (X@Vt[:k].T)@Vt[:k] if k>0 else X
    C=(Xr[:,tidx].T@Xr)/Xr.shape[0]
    key=-np.sign(C); agree=(key==np.sign(Hm))
    print(f"{k:>12d} {acc(key,pool):>13.4f} {acc(key,weak):>19.4f} {acc(np.sign(Hm),pool&agree):>12.4f} {acc(np.sign(Hm),pool&~agree):>15.4f} {np.median(np.abs(C)):>11.4f}")
    if k==10: C10=C
# with the conditioned correlation: combine with cons on the weak half
key10=-np.sign(C10)
print(f"\nwith 10 PCs removed, on the WEAK-donor half:  cons {acc(np.sign(consm),weak):.4f}   H {acc(np.sign(Hm),weak):.4f}   -sign(r) {acc(key10,weak):.4f}")
for lam in (0.5,1,2,4):
    comb=np.sign(consm/ (np.abs(consm).mean()+1e-12) - lam*C10/(np.abs(C10).mean()+1e-12))
    print(f"   sign(cons_std - {lam}*r_std):  pool {acc(comb,pool):.4f}   weak half {acc(comb,weak):.4f}")
