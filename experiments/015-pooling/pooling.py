"""Experiment 015 — per-gene inverse-variance pooling across donors.

Everything so far blends donors with ONE weight per donor, applied to every gene. But donors
differ in how well they measure each *gene*: a gene with 200 counts in K562 and 8 in HCT116 is
well-measured in one and noise in the other, and a fixed 8:1:1 weight ignores that entirely.

Inverse-variance pooling is the minimum-variance estimator when combining unbiased estimates:

    blended_g = sum_d (delta_dg / var_dg) / sum_d (1 / var_dg)

with Poisson variance var_dg ~ 1/(n_cells_d * depth * p_dg) from each donor's own control profile.

This is NOT the EB shrinkage that failed in exp 009/012. That scaled magnitudes down (and cost
pds); this reallocates weight *between donors* per gene while leaving the overall scale alone.
Sidechain (public, and holder of the board's best `fid`) does exactly this.

Also tested: cross-donor magnitude dispersion as a confidence feature. Nadig 2025 found
cell-type-preferential perturbations have ~2x larger cross-line differences in effect magnitude
(p=4e-8) -- a signal orthogonal to the directional agreement of exp 011.
"""
import sys, json
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
    e=np.expm1(d.b_ctrl); pg=e/max(e.sum(),1e-12)
    ncell=json.loads(Path(p+"/meta.json").read_text())["cells_per_pert"]["median"]
    banks[n]=dict(s=s,gpos={g:i for i,g in enumerate(gl)},gset=set(gl),
                  idx={q:i for i,q in enumerate(d.perts)},
                  var=1.0/np.maximum(ncell*20000.0*pg,1e-3))
vcc=[g for g in pd.read_csv("data/challenge/pert_counts.csv")["target_gene"]]
common=[g for g in vcc if all(g in b["idx"] for b in banks.values())]

def pds(pred,true):
    P=pred/(np.linalg.norm(pred,axis=1,keepdims=True)+1e-12)
    T=true/(np.linalg.norm(true,axis=1,keepdims=True)+1e-12)
    D=1.0-P@T.T; n=D.shape[1]; sd=np.diag(D)[:,None]
    k=(D<sd).sum(1)+0.5*((D==sd).sum(1)-1)
    return float((1.0-k/(n-1)).mean())

print(f"{'scheme':34} {'macro PDS':>10}")
res={}
for scheme in ("fixed 8:1:1","fixed 4:1:1","inverse-variance","inv-var x fixed 8:1:1"):
    folds=[]
    for held in banks:
        others=[n for n in banks if n!=held]
        shared=sorted(set.intersection(*[banks[n]["gset"] for n in [held]+others]))
        feat=[g for g in shared if g not in set(common)]
        col={n:np.fromiter((banks[n]["gpos"][g] for g in feat),np.int64,len(feat)) for n in [held]+others}
        T=np.vstack([banks[held]["s"][banks[held]["idx"][p]][col[held]] for p in common])
        S={n:np.vstack([banks[n]["s"][banks[n]["idx"][p]][col[n]] for p in common]) for n in others}
        for n in others: S[n]=S[n]/(np.median(np.linalg.norm(S[n],axis=1)) or 1.0)
        V={n:banks[n]["var"][col[n]] for n in others}
        if scheme.startswith("fixed"):
            w=[float(x) for x in scheme.split()[1].split(":")][:2]
            P=(w[0]*S[others[0]]+w[1]*S[others[1]])/sum(w)
        else:
            iv={n:1.0/np.maximum(V[n],1e-12) for n in others}
            if "x fixed" in scheme:
                fw={others[0]:8.0,others[1]:1.0}
                iv={n:iv[n]*fw[n] for n in others}
            tot=sum(iv.values())
            P=sum(S[n]*(iv[n]/tot)[None,:] for n in others)
        folds.append(pds(P,T))
    res[scheme]=float(np.mean(folds))
    print(f"  {scheme:32} {res[scheme]:10.4f}")
print(f"\nbest: {max(res,key=res.get)} -> {max(res.values()):.4f}")

# Nadig's magnitude-dispersion signal
print("\nNadig magnitude-dispersion vs directional agreement, as confidence features:")
ag,md,tr=[],[],[]
u=lambda M: M/(np.linalg.norm(M,axis=1,keepdims=True)+1e-12)
for held in banks:
    others=[n for n in banks if n!=held]
    shared=sorted(set.intersection(*[banks[n]["gset"] for n in [held]+others]))
    feat=[g for g in shared if g not in set(common)]
    col={n:np.fromiter((banks[n]["gpos"][g] for g in feat),np.int64,len(feat)) for n in [held]+others}
    T=np.vstack([banks[held]["s"][banks[held]["idx"][p]][col[held]] for p in common])
    A=np.vstack([banks[others[0]]["s"][banks[others[0]]["idx"][p]][col[others[0]]] for p in common])
    B=np.vstack([banks[others[1]]["s"][banks[others[1]]["idx"][p]][col[others[1]]] for p in common])
    ag.append(np.sum(u(A)*u(B),axis=1))
    na,nb=np.linalg.norm(A,axis=1),np.linalg.norm(B,axis=1)
    md.append(-np.abs(na-nb)/(na+nb+1e-12))          # negative dispersion = agreement
    tr.append(np.sum(u(u(A)+u(B))*u(T),axis=1))
ag,md,tr=np.concatenate(ag),np.concatenate(md),np.concatenate(tr)
print(f"  Spearman(directional agreement, transfer) = {spearmanr(ag,tr)[0]:+.3f}")
print(f"  Spearman(-magnitude dispersion, transfer) = {spearmanr(md,tr)[0]:+.3f}")
print(f"  Spearman(agreement, -magnitude dispersion) = {spearmanr(ag,md)[0]:+.3f}   (are they redundant?)")
