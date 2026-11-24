"""Experiment 012 — turn cross-donor agreement into per-perturbation magnitude.

Experiment 011 showed agreement between donor lines predicts how well a perturbation transfers
(Spearman 0.324, Q5/Q1 = 3.2x). This asks the question that actually matters: **how much of the
achievable magnitude-calibration gain can agreement capture?**

The theory. For perturbation p with predicted direction d_p and true effect delta_p, the
squared-error-optimal magnitude is

    alpha*_p = <d_p, delta_p> / <d_p, d_p>

and at that magnitude the residual is ||delta_p||^2 (1 - cos^2_p). Using one global alpha for
every perturbation is strictly worse whenever cos_p varies -- and it varies a lot. But alpha*_p
needs the truth, which we do not have on A/B/C.

Agreement is a proxy for cos_p that needs no truth. So: fit alpha_hat_p from agreement on two
donor lines, apply it to a third, and measure what fraction of the oracle gap it closes.

Three arms:
  uniform   one global alpha (what we submit today)
  oracle    per-perturbation alpha*_p -- upper bound, uses the held-out truth, not usable
  agreement per-perturbation alpha_hat_p fitted from donor agreement -- usable on A/B/C

Also reports pds for each, which should be *unchanged*: cosine distance is invariant to scaling
one perturbation's own effect, so per-perturbation magnitude is a free move on our best metric.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, _norm_scale

DONORS={"k562_gwps":"data/donors/k562_gwps","orion_hct116":"data/donors/orion_hct116",
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
print(f"{len(common)} panel perturbations in all three lines\n")

def pds(pred,true):
    P=pred/(np.linalg.norm(pred,axis=1,keepdims=True)+1e-12)
    T=true/(np.linalg.norm(true,axis=1,keepdims=True)+1e-12)
    D=1.0-P@T.T; n=D.shape[1]; sd=np.diag(D)[:,None]
    k=(D<sd).sum(1)+0.5*((D==sd).sum(1)-1)
    return float((1.0-k/(n-1)).mean())

u=lambda M: M/(np.linalg.norm(M,axis=1,keepdims=True)+1e-12)
folds=[]
for held in banks:
    others=[n for n in banks if n!=held]
    shared=sorted(set.intersection(*[banks[n]["gset"] for n in [held]+others]))
    feat=[g for g in shared if g not in set(common)]
    col={n:np.fromiter((banks[n]["gpos"][g] for g in feat),np.int64,len(feat)) for n in [held]+others}
    T=np.vstack([banks[held]["s"][banks[held]["idx"][p]][col[held]] for p in common])
    A=np.vstack([banks[others[0]]["s"][banks[others[0]]["idx"][p]][col[others[0]]] for p in common])
    B=np.vstack([banks[others[1]]["s"][banks[others[1]]["idx"][p]][col[others[1]]] for p in common])
    agree=np.sum(u(A)*u(B),axis=1)
    # 4:1 blend, matching the LODO optimum
    D_=4.0*u(A)+1.0*u(B); D_=D_/(np.linalg.norm(D_,axis=1,keepdims=True)+1e-12)
    folds.append(dict(held=held,D=D_,T=T,agree=agree))

print(f"{'held-out':16} {'uniform':>9} {'agreement':>10} {'oracle':>9} {'captured':>9}  "
      f"{'pds unif':>9} {'pds agree':>10}")
rows=[]
for i,f in enumerate(folds):
    D_,T,agree=f["D"],f["T"],f["agree"]
    dd=np.sum(D_*D_,axis=1); dt=np.sum(D_*T,axis=1)
    a_star=dt/np.maximum(dd,1e-12)                    # oracle per-pert alpha
    a_uni=float(dt.sum()/dd.sum())                    # best single global alpha

    # Fit agreement -> alpha on the OTHER folds only; never on this one.
    tr=np.concatenate([folds[j]["agree"] for j in range(len(folds)) if j!=i])
    ta=np.concatenate([(np.sum(folds[j]["D"]*folds[j]["T"],axis=1)/
                        np.maximum(np.sum(folds[j]["D"]*folds[j]["D"],axis=1),1e-12))
                       for j in range(len(folds)) if j!=i])
    b1,b0=np.polyfit(tr,ta,1)
    a_hat=np.clip(b0+b1*agree, 0.0, None)
    a_hat=a_hat*(a_uni/max(a_hat.mean(),1e-12))       # preserve mean magnitude

    err=lambda a: float((np.linalg.norm(a[:,None]*D_-T,axis=1)**2).sum())
    e_u,e_a,e_o=err(np.full(len(T),a_uni)),err(a_hat),err(a_star)
    cap=(e_u-e_a)/(e_u-e_o) if e_u>e_o else 0.0
    rows.append(cap)
    print(f"  {f['held']:14} {e_u:9.1f} {e_a:10.1f} {e_o:9.1f} {100*cap:8.1f}%  "
          f"{pds(a_uni*D_,T):9.4f} {pds(a_hat[:,None]*D_,T):10.4f}")

print(f"\n  mean fraction of the oracle magnitude gain captured by agreement: {100*np.mean(rows):.1f}%")
print("  pds columns should match to ~4dp -- per-perturbation scaling is cosine-invariant.")
allag=np.concatenate([f["agree"] for f in folds])
allst=np.concatenate([np.sum(f["D"]*f["T"],axis=1)/np.maximum(np.sum(f["D"]*f["D"],axis=1),1e-12)
                      for f in folds])
r,p=spearmanr(allag,allst)
print(f"\n  Spearman(agreement, oracle alpha*) = {r:.3f}  (p={p:.1e})")
print(f"  oracle alpha* spread: p10 {np.percentile(allst,10):.3f}  median {np.median(allst):.3f}  "
      f"p90 {np.percentile(allst,90):.3f}")
