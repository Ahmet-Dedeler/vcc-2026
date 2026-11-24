import numpy as np, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from feats172 import FEATURES
from fit172 import PROD, Model, load_folds
from harness172 import H1Test
from run172 import train, prod_key
from scipy.stats import spearmanr
T=H1Test(); F=load_folds()
Xall=np.concatenate([F[h]["X"] for h in PROD]); yall=np.concatenate([F[h]["y"] for h in PROD])
tall=np.concatenate([F[h]["t"]+100000*i for i,h in enumerate(PROD)])
hf=np.load("h1test-feats.npz") if os.path.exists("h1test-feats.npz") else np.load(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"work","h1test-feats.npz"))
from fit172 import h1_features
Xh,_=h1_features(T); shape=T.H.shape
def show(tag, m, cols=None):
    XX = Xh if cols is None else Xh[:,cols]
    c=m.predict(XX).reshape(shape)
    co=m.clf.coef_.ravel()
    names=[FEATURES[i] for i in (cols if cols is not None else range(len(FEATURES)))]
    top=sorted(zip(names,co),key=lambda z:-abs(z[1]))[:6]
    sp=np.median([spearmanr(c[p,i],np.abs(hf["cons"])[p,i]).statistic for p in range(85) if (i:=np.flatnonzero(T.pool[p])).size>5])
    print(f"{tag:<26} reach {T.reach(c):.5f} d1 {T.deciles(c)[0]:.3f} spread {c.std():.3e} "
          f"rho_to|cons| {sp:+.3f}  top coefs " + " ".join(f"{n}={v:+.3f}" for n,v in top))
cfg=dict(kind="logistic",C=0.1); w="top25"
show("fitted (real labels)", train(cfg,w,Xall,yall,tall))
for seed in range(3):
    rng=np.random.default_rng(seed); yp=yall.copy()
    for p in np.unique(tall):
        m=tall==p; yp[m]=rng.permutation(yall[m])
    print(f"  seed {seed}: shuffled y mean {yp.mean():.4f} (orig {yall.mean():.4f}) "
          f"corr(y,rank_cons) real {np.corrcoef(yall,Xall[:,FEATURES.index('rank_cons')])[0,1]:+.5f} "
          f"shuf {np.corrcoef(yp,Xall[:,FEATURES.index('rank_cons')])[0,1]:+.5f}")
    show(f"null seed {seed}", train(cfg,w,Xall,yp,tall,seed=seed))
# zero-coefficient control: what does a literally constant key read?
print(f"constant key: reach {T.reach(np.zeros(shape)):.5f} d1 {T.deciles(np.zeros(shape))[0]:.3f}")
