"""Experiment 014 — is `jac` escapable? Everyone on the board is at chance.

Board fact (84 entries, verified): raw `de_wilcoxon_sig_jaccard` is 0.029-0.033 for every team,
against a baseline of 0.030. Nobody has moved it.

The mechanism. J = |R n R_hat| / |R u R_hat|. With minimum-variance ("even") emission the
Wilcoxon calls essentially the whole 5-CPM gate, so R_hat ~ 10,000 genes while the reference's
R is a few hundred: J -> |R|/10,000 ~ 0.03. That is exactly the baseline, and it is a
*self-inflicted* ceiling, not a difficulty of the biology.

If instead we call exactly n_real genes, chosen as our most confident, then

    J = p / (2 - p)   where p = precision = |R n R_hat| / n_real

so p = 0.3 gives J = 0.176 (scaled 0.40, +0.066 overall) and p = 0.5 gives J = 0.333
(scaled 0.82, +0.137 overall). The whole gap from rank 5 to rank 1 is 0.0515.

**Everything therefore hinges on one number: what precision can donor transfer achieve at
k = n_real?** This measures it, held-out, with no target information.

Note this is a different quantity from sign accuracy. `jac` asks *which genes respond*, not
which direction they move.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, _norm_scale

B_JAC, R_JAC = 0.030, 0.399     # official anchors
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

def scaled(J): return (J-B_JAC)/(R_JAC-B_JAC)

print("Precision of donor-predicted DE sets at k = n_real, held out\n")
print(f"{'held-out':16} {'n_real':>7} {'precision':>10} {'J=p/(2-p)':>11} {'scaled':>8} {'overall':>9}")
rows=[]
for held in banks:
    others=[n for n in banks if n!=held]
    shared=sorted(set.intersection(*[banks[n]["gset"] for n in [held]+others]))
    feat=[g for g in shared if g not in set(common)]
    col={n:np.fromiter((banks[n]["gpos"][g] for g in feat),np.int64,len(feat)) for n in [held]+others}
    for n_real in (300,):
        precs=[]
        for p in common:
            T=banks[held]["s"][banks[held]["idx"][p]][col[held]]
            A=banks[others[0]]["s"][banks[others[0]]["idx"][p]][col[others[0]]]
            B=banks[others[1]]["s"][banks[others[1]]["idx"][p]][col[others[1]]]
            blend=4.0*A/(np.linalg.norm(A)+1e-12)+1.0*B/(np.linalg.norm(B)+1e-12)
            R=set(np.argsort(-np.abs(T))[:n_real].tolist())          # "reference significant"
            Rhat=set(np.argsort(-np.abs(blend))[:n_real].tolist())   # our top-k
            precs.append(len(R&Rhat)/n_real)
        pr=float(np.mean(precs)); J=pr/(2-pr)
        rows.append((held,n_real,pr,J))
        print(f"  {held:14} {n_real:7d} {pr:10.4f} {J:11.4f} {scaled(J):8.3f} {scaled(J)/6:+9.4f}")

pr=np.mean([r[2] for r in rows]); J=pr/(2-pr)
print(f"\n  MACRO precision {pr:.4f} -> J {J:.4f} -> scaled {scaled(J):.3f} -> overall {scaled(J)/6:+.4f}")
print(f"  (our current jac scaled is +0.0036; the whole gap to rank 1 is 0.0515)")

print("\nSensitivity — how many genes should we call?")
held=list(banks)[0]; others=[n for n in banks if n!=held]
shared=sorted(set.intersection(*[banks[n]["gset"] for n in [held]+others]))
feat=[g for g in shared if g not in set(common)]
col={n:np.fromiter((banks[n]["gpos"][g] for g in feat),np.int64,len(feat)) for n in [held]+others}
print(f"  {'k called':>9} {'precision':>10} {'recall':>8} {'J':>8} {'scaled':>8}")
for k in (50,100,200,300,500,1000,3000,10000):
    if k>len(feat): break
    js=[]
    for p in common[:120]:
        T=banks[held]["s"][banks[held]["idx"][p]][col[held]]
        A=banks[others[0]]["s"][banks[others[0]]["idx"][p]][col[others[0]]]
        B=banks[others[1]]["s"][banks[others[1]]["idx"][p]][col[others[1]]]
        blend=4.0*A/(np.linalg.norm(A)+1e-12)+1.0*B/(np.linalg.norm(B)+1e-12)
        R=set(np.argsort(-np.abs(T))[:300].tolist())
        Rhat=set(np.argsort(-np.abs(blend))[:k].tolist())
        inter=len(R&Rhat); js.append((inter/k, inter/300, inter/len(R|Rhat)))
    a=np.array(js).mean(0)
    print(f"  {k:9d} {a[0]:10.4f} {a[1]:8.4f} {a[2]:8.4f} {scaled(a[2]):8.3f}")
