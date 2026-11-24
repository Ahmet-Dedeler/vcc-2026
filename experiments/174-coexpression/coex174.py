"""174 -- does target<->gene co-expression across the RECIPIENT's control cells predict the sign
of the perturbation response?  A (target, gene) signal from the recipient that no donor carries.
h1-test, 18,400 controls; truth = reference DE table (164's de164.a15.npz)."""
import numpy as np, sys, time, h5py, scipy.sparse as sp
sys.path.insert(0, "experiments/153-de-bench"); sys.path.insert(0, "experiments/168-sharpening")
from h1_tiebreak import labels
from debench import Surface
t0=time.time()
S = Surface("experiments/164-callset-artefact/de164.a15.npz")
E = np.load("experiments/169-sharpening-six-member/work/h1test-effect.npz", allow_pickle=True)
F = np.load("experiments/172-sign-confidence-key/work/h1test-feats.npz", allow_pickle=True)
eg = {str(g): i for i, g in enumerate(E["genes"])}; et = {str(t): i for i, t in enumerate(E["targets"])}
gi = np.array([eg[g] for g in S.genes]); pi = np.array([et[p] for p in S.perts])
H = E["H"][pi][:, gi].astype(float); cons = F["cons"].astype(float)

# ---- control cells, log-CPM, dense (18400 x 18080 float32 = 1.3 GB) ----
with h5py.File("data/harness/h1-test.h5ad","r") as h:
    genes=labels(h,"var/_index"); lab=labels(h,"obs/target_gene")
    ci=np.flatnonzero(lab=="non-targeting"); n,G=(int(x) for x in h["X"].attrs["shape"])
    ip=h["X/indptr"][:]; rows=[]
    for r0 in range(0,ci.size,3000):
        idx=ci[r0:r0+3000]; blk=np.zeros((idx.size,G),dtype=np.float32)
        for k,r in enumerate(idx):
            s,e=int(ip[r]),int(ip[r+1]); blk[k,h["X/indices"][s:e]]=h["X/data"][s:e]
        tot=blk.sum(1,keepdims=True); tot[tot==0]=1
        rows.append(np.log1p(blk*(1e6/tot)))
    X=np.concatenate(rows); del rows
assert np.array_equal(genes, np.asarray([str(g) for g in E["genes"]]))
print(f"controls {X.shape} log-CPM [{time.time()-t0:.0f}s]")
X-=X.mean(0,keepdims=True); sd=X.std(0); sd[sd==0]=1; X/=sd
tpos={g:i for i,g in enumerate(genes)}
tidx=np.array([tpos[str(p)] for p in S.perts])
C = (X[:,tidx].T @ X)/X.shape[0]          # (85 targets x 18080 genes) Pearson r across cells
C = C[:, gi]                              # onto the DE-universe axis
print(f"co-expression matrix {C.shape}; |r| median {np.median(np.abs(C)):.4f}, q99 {np.quantile(np.abs(C),.99):.4f} [{time.time()-t0:.0f}s]")

pool = S.real_sig & S.adj
sT = np.sign(S.real_lfc)
for name,key in [("-sign(r)  [t activates g -> knockdown lowers g]", -np.sign(C)), ("+sign(r)", np.sign(C))]:
    print(f"\nsign accuracy on the reference pool, {name}: {np.mean(key[pool]==sT[pool]):.4f}")
absr=np.abs(C)
print("\n-sign(r) accuracy on the pool by |r| quantile (does STRONGER co-expression predict better?):")
q=np.quantile(absr[pool],[0,.5,.75,.9,.95,.99,1])
for a,b in zip(q[:-1],q[1:]):
    m=pool&(absr>=a)&(absr<b); print(f"   |r| in [{a:.4f},{b:.4f})  n={m.sum():6d}  acc {np.mean(-np.sign(C)[m]==sT[m]):.4f}")
# where the donors are UNINFORMATIVE (bottom half of |cons|), does co-expression add anything?
ac=np.abs(cons); med=np.median(ac[pool&(ac>0)])
weak=pool&(ac<med); strong=pool&(ac>=med)
print(f"\nwhere donor consensus is WEAK (|cons| < median):  cons acc {np.mean(np.sign(cons)[weak]==sT[weak]):.4f}   H acc {np.mean(np.sign(H)[weak]==sT[weak]):.4f}   -sign(r) acc {np.mean(-np.sign(C)[weak]==sT[weak]):.4f}")
print(f"where donor consensus is STRONG:                   cons acc {np.mean(np.sign(cons)[strong]==sT[strong]):.4f}   H acc {np.mean(np.sign(H)[strong]==sT[strong]):.4f}   -sign(r) acc {np.mean(-np.sign(C)[strong]==sT[strong]):.4f}")
# top-|r| genes per target: accuracy of -sign(r) there
print("\n-sign(r) accuracy on the pool restricted to each target's top-K genes by |r| (the co-expression 'module'):")
for K in (50,100,250,500,1000):
    m=np.zeros_like(pool)
    for p in range(85):
        top=np.argsort(-absr[p])[:K]; m[p,top]=True
    mm=m&pool; print(f"   top {K:4d} by |r|:  pool genes covered {mm.sum():5d}   acc {np.mean(-np.sign(C)[mm]==sT[mm]):.4f}   (H acc there {np.mean(np.sign(H)[mm]==sT[mm]):.4f})")
# does agreement between co-expression sign and H sign mark reliable genes?
agree=(-np.sign(C)==np.sign(H))
print(f"\nH acc where -sign(r) AGREES with H: {np.mean(np.sign(H)[pool&agree]==sT[pool&agree]):.4f} (n={int((pool&agree).sum())})   where it DISAGREES: {np.mean(np.sign(H)[pool&~agree]==sT[pool&~agree]):.4f} (n={int((pool&~agree).sum())})")
# reach counterfactual with |r| as key and with |H|*agree
absH=np.abs(H)
print(f"\nreach key-counterfactual:  |H| {S.reach(key=absH):.5f}   |r| {S.reach(key=absr):.5f}   |H|*(1+agree) {S.reach(key=absH*(1+agree)):.5f}   |H|*|r| {S.reach(key=absH*absr):.5f}")
np.savez_compressed("experiments/174-coexpression/coex174.npz", C=C.astype(np.float32), genes=S.genes, perts=S.perts)
