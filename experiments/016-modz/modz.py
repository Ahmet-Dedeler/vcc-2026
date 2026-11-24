"""Experiment 016 — MODZ per-perturbation donor weighting, on a 4-line LODO.

MODZ (CMap L1000, `cmapPy/math/agg_wt_avg.py`) weights each replicate by its mean correlation
with the *other* replicates, per signature. Substituting donors for replicates gives exactly
"weight each donor by how much it agrees with the others, for this perturbation" -- the
published prior art for experiment 011's finding, and the top recommendation from the research
sweep.

It was untestable until now: MODZ needs >=3 sources to differentiate (with 2, corr(A,B) is a
single number and the weights are always equal), and 3-line LODO leaves only 2 predictors.
Adding H1-2025 as a fourth line gives 3 predictors per fold on 226 shared perturbations.

    w_d = mean_{j != d} max(0, spearman(v_d, v_j)),  floored at 0.01, normalised
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, _norm_scale

LINES={"k562_gwps":"data/donors/k562_gwps","orion_hct116":"data/donors/orion_hct116",
       "orion_hek293t":"data/donors/orion_hek293t","h1_2025":"data/donors/h1_2025"}
banks={}
for n,p in LINES.items():
    d=Donor.load(p); m,s=d.decompose()
    s=(s/_norm_scale(s)).astype(np.float32)
    gl=[str(g) for g in d.genes]
    banks[n]=dict(s=s,gpos={g:i for i,g in enumerate(gl)},gset=set(gl),
                  idx={q:i for i,q in enumerate(d.perts)})
common=sorted(set.intersection(*[set(b["idx"]) for b in banks.values()]))
vcc=set(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"])
print(f"{len(common)} perturbations shared across all 4 lines "
      f"({len(set(common)&vcc)} of them on the VCC panel)\n")

def pds(pred,true):
    P=pred/(np.linalg.norm(pred,axis=1,keepdims=True)+1e-12)
    T=true/(np.linalg.norm(true,axis=1,keepdims=True)+1e-12)
    D=1.0-P@T.T; n=D.shape[1]; sd=np.diag(D)[:,None]
    k=(D<sd).sum(1)+0.5*((D==sd).sum(1)-1)
    return float((1.0-k/(n-1)).mean())

def modz_weights(V):
    """V: (n_donors, n_genes) for ONE perturbation. Returns MODZ weights."""
    k=len(V)
    C=np.zeros((k,k))
    for i in range(k):
        for j in range(i+1,k):
            C[i,j]=C[j,i]=spearmanr(V[i],V[j])[0]
    np.fill_diagonal(C,np.nan)
    raw=np.nanmean(np.clip(C,0,None),axis=1)
    raw=np.maximum(raw,0.01)
    return raw/raw.sum()

schemes=["equal","primary 4:1:1","MODZ","MODZ x primary"]
acc={s:[] for s in schemes}
print(f"{'held out':16} " + "  ".join(f"{s:>14}" for s in schemes))
for held in banks:
    others=[n for n in banks if n!=held]
    shared=sorted(set.intersection(*[banks[n]["gset"] for n in [held]+others]))
    feat=[g for g in shared if g not in set(common)]
    col={n:np.fromiter((banks[n]["gpos"][g] for g in feat),np.int64,len(feat)) for n in [held]+others}
    T=np.vstack([banks[held]["s"][banks[held]["idx"][p]][col[held]] for p in common])
    S={n:np.vstack([banks[n]["s"][banks[n]["idx"][p]][col[n]] for p in common]) for n in others}
    for n in others: S[n]=S[n]/(np.median(np.linalg.norm(S[n],axis=1)) or 1.0)
    # fixed weights: primary = whichever of the remaining lines is k562 if present, else first
    prim = "k562_gwps" if "k562_gwps" in others else others[0]
    fw=np.array([4.0 if n==prim else 1.0 for n in others]); fw/=fw.sum()
    row={}
    row["equal"]=pds(np.mean([S[n] for n in others],axis=0),T)
    row["primary 4:1:1"]=pds(sum(fw[i]*S[others[i]] for i in range(len(others))),T)
    Pm=np.zeros_like(T); Pmx=np.zeros_like(T)
    for r,p in enumerate(common):
        V=np.vstack([S[n][r] for n in others])
        w=modz_weights(V)
        Pm[r]=w@V
        wx=w*fw; wx/=wx.sum(); Pmx[r]=wx@V
    row["MODZ"]=pds(Pm,T); row["MODZ x primary"]=pds(Pmx,T)
    for s in schemes: acc[s].append(row[s])
    print(f"  {held:14} " + "  ".join(f"{row[s]:14.4f}" for s in schemes))

print(f"\n{'MACRO':16} " + "  ".join(f"{np.mean(acc[s]):14.4f}" for s in schemes))
best=max(schemes,key=lambda s:np.mean(acc[s]))
print(f"\nbest: {best} -> {np.mean(acc[best]):.4f}   "
      f"vs fixed primary {np.mean(acc['primary 4:1:1']):.4f} "
      f"({np.mean(acc[best])-np.mean(acc['primary 4:1:1']):+.4f})")
