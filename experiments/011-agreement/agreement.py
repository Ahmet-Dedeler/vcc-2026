"""Experiment 011 — is cross-donor agreement a per-perturbation confidence signal?

Experiment 003 established that some perturbations transfer well across lines and others do
not (Nadig 2025: 44% are cell-type-preferential), and that donor effect **magnitude** does not
predict which (Spearman 0.063). That killed the obvious confidence weighting.

But we now have three donor lines. If K562, HCT116 and HEK293T *agree* on a perturbation's
response, that response is plausibly conserved machinery and should transfer to a fourth line.
If they disagree, it is line-specific and asserting it is a bet we should not take.

The key property: agreement among predictors is computable **without the target**, so whatever
it tells us here is usable on A/B/C where there is no ground truth at all.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, _norm_scale

DONORS = {"k562_gwps":"data/donors/k562_gwps","orion_hct116":"data/donors/orion_hct116",
          "orion_hek293t":"data/donors/orion_hek293t"}

banks={}
for n,p in DONORS.items():
    d=Donor.load(p); m,s=d.decompose()
    s=(s/_norm_scale(s)).astype(np.float32)
    gl=[str(g) for g in d.genes]
    banks[n]=dict(s=s,gpos={g:i for i,g in enumerate(gl)},gset=set(gl),
                  idx={q:i for i,q in enumerate(d.perts)})

vcc=[g for g in pd.read_csv("data/challenge/pert_counts.csv")["target_gene"]]
common=[g for g in vcc if all(g in b["idx"] for b in banks.values())]
print(f"{len(common)} panel perturbations present in all three lines\n")

u=lambda M: M/(np.linalg.norm(M,axis=1,keepdims=True)+1e-12)
rows=[]
for held in banks:
    others=[n for n in banks if n!=held]
    shared=sorted(set.intersection(*[banks[n]["gset"] for n in [held]+others]))
    panel=set(common)
    feat=[g for g in shared if g not in panel]
    col={n:np.fromiter((banks[n]["gpos"][g] for g in feat),np.int64,len(feat))
         for n in [held]+others}
    T=np.vstack([banks[held]["s"][banks[held]["idx"][p]][col[held]] for p in common])
    A=np.vstack([banks[others[0]]["s"][banks[others[0]]["idx"][p]][col[others[0]]] for p in common])
    B=np.vstack([banks[others[1]]["s"][banks[others[1]]["idx"][p]][col[others[1]]] for p in common])
    agree=np.sum(u(A)*u(B),axis=1)
    truth=np.sum(u(u(A)+u(B))*u(T),axis=1)
    rows.append(pd.DataFrame({"held":held,"pert":common,"agreement":agree,"transfer":truth}))

df=pd.concat(rows)
print(f"{'held-out line':18} {'n':>5} {'Spearman(agree,transfer)':>26} {'p':>10}")
for held,g in df.groupby("held"):
    rho,p=spearmanr(g.agreement,g.transfer)
    print(f"  {held:16} {len(g):5d} {rho:26.3f} {p:10.2e}")
rho,p=spearmanr(df.agreement,df.transfer)
print(f"  {'POOLED':16} {len(df):5d} {rho:26.3f} {p:10.2e}")

print("\nTransfer by agreement quintile (pooled):")
df["q"]=pd.qcut(df.agreement,5,labels=[1,2,3,4,5])
for q,g in df.groupby("q",observed=True):
    print(f"  Q{q}  agreement [{g.agreement.min():+.3f},{g.agreement.max():+.3f}]  "
          f"n={len(g):4d}  mean transfer {g.transfer.mean():+.4f}")
lo=df[df.q==1].transfer.mean(); hi=df[df.q==5].transfer.mean()
print(f"\n  Q5 - Q1 = {hi-lo:+.4f}")
df.to_csv("experiments/011-agreement/agreement.csv",index=False)
