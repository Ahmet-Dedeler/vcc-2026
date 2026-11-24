import numpy as np, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from feats172 import FEATURES
from fit172 import PROD, load_folds
from run172 import train
F=load_folds()
X=np.concatenate([F[h]["X"] for h in PROD]); y=np.concatenate([F[h]["y"] for h in PROD])
t=np.concatenate([F[h]["t"]+100000*i for i,h in enumerate(PROD)])
m=train(dict(kind="logistic",C=0.1),"top25",X,y,t)
co=m.clf.coef_.ravel()
print("M1 logistic @ C=0.1 w=top25 — standardised coefficients (log-odds per 1 sd), all 3 folds")
print(f"{'feature':<16}{'coef':>9}{'|coef|rank':>12}   univariate corr(y, feature) pooled over folds")
uc={f: np.corrcoef(y, X[:,i])[0,1] for i,f in enumerate(FEATURES)}
for r,(f,c) in enumerate(sorted(zip(FEATURES,co), key=lambda z:-abs(z[1])),1):
    print(f"{f:<16}{c:>+9.3f}{r:>12}{uc[f]:>+34.4f}")
