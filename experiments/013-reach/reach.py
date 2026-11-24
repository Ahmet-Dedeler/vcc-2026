"""Experiment 013 — can we rank genes better for `reach`?

`de_wilcoxon_direction_reach_raw` has the widest span of any metric (baseline 0.072, replicate
0.968) and we sit at 0.042 scaled. It asks: take the genes the reference called significant,
order them by OUR confidence, walk down, and find the deepest prefix whose predicted signs are
>=90% correct. Score = that depth / n_real.

So `reach` is purely a **ranking** question. Magnitude does not enter; only the order in which we
would stand behind our own sign calls. Today that order falls out of whatever p-values the
emitted counts produce, which is an accident of the emission model rather than a choice.

Candidate confidence signals, all computable without the target:
  |effect|          predicted effect magnitude (proxy for what we emit today)
  agreement         per-gene sign agreement between donors
  min|donor|        the weaker donor's magnitude -- a conservative floor
  |effect| x agree  magnitude gated by agreement
  geometric         sqrt(|A| * |B|) where both donors point the same way, else 0

Evaluated in LODO: hold out a line, take its top-N genes by |true effect| as a stand-in for
"reference called significant", rank them by each candidate, and compute reach exactly as the
spec defines it.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, _norm_scale

PURITY = 0.9
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

def reach_one(conf, pred_sign, true_sign):
    """Deepest prefix (by descending conf) whose sign purity stays >= PURITY."""
    o=np.argsort(-conf)
    match=(pred_sign[o]==true_sign[o]).astype(float)
    cum=np.cumsum(match)/np.arange(1,len(match)+1)
    ok=np.flatnonzero(cum>=PURITY)
    return (ok[-1]+1)/len(match) if ok.size else 0.0

N_DE=200   # stand-in for the reference's significant set
print(f"reach with top-{N_DE} genes by |true effect| as the reference set, purity>={PURITY}\n")
sig=["|effect|","agreement","min|donor|","|effect| x agree","geometric","random"]
acc={s:[] for s in sig}
rng=np.random.default_rng(0)

for held in banks:
    others=[n for n in banks if n!=held]
    shared=sorted(set.intersection(*[banks[n]["gset"] for n in [held]+others]))
    feat=[g for g in shared if g not in set(common)]
    col={n:np.fromiter((banks[n]["gpos"][g] for g in feat),np.int64,len(feat)) for n in [held]+others}
    per={s:[] for s in sig}
    for p in common:
        T=banks[held]["s"][banks[held]["idx"][p]][col[held]]
        A=banks[others[0]]["s"][banks[others[0]]["idx"][p]][col[others[0]]]
        B=banks[others[1]]["s"][banks[others[1]]["idx"][p]][col[others[1]]]
        blend=4.0*A/(np.linalg.norm(A)+1e-12)+1.0*B/(np.linalg.norm(B)+1e-12)
        top=np.argsort(-np.abs(T))[:N_DE]
        ts,ps=np.sign(T[top]),np.sign(blend[top])
        agree=(np.sign(A[top])==np.sign(B[top])).astype(float)
        cands={
            "|effect|":       np.abs(blend[top]),
            "agreement":      agree + 1e-6*np.abs(blend[top]),
            "min|donor|":     np.minimum(np.abs(A[top]),np.abs(B[top])),
            "|effect| x agree": np.abs(blend[top])*(0.1+agree),
            "geometric":      np.sqrt(np.abs(A[top])*np.abs(B[top]))*agree,
            "random":         rng.random(len(top)),
        }
        for s,c in cands.items(): per[s].append(reach_one(c,ps,ts))
    for s in sig: acc[s].append(np.mean(per[s]))
    print(f"  {held:16} " + "  ".join(f"{s}={np.mean(per[s]):.4f}" for s in sig))

print(f"\n{'signal':20} {'macro reach':>12}")
for s in sorted(sig,key=lambda x:-np.mean(acc[x])):
    print(f"  {s:18} {np.mean(acc[s]):12.4f}")
base=np.mean(acc["|effect|"]); best=max(sig,key=lambda x:np.mean(acc[x]))
print(f"\nbest: {best} -> {np.mean(acc[best]):.4f}  vs |effect| {base:.4f}  "
      f"({100*(np.mean(acc[best])-base)/max(base,1e-9):+.1f}%)")
