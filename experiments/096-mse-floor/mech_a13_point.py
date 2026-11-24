"""Exact operating point for the shrinkage arm on top of a12, on the CHALLENGE panel.

A'(gamma) and B'(gamma) are computed exactly from a10's shipped group sums plus the real control.
X'(gamma) needs the truth, so only the dimensionless ratio X'/X is transferred, measured on BOTH
harnesses and taken at the rule-1 (smaller) value. The optimum over s is found through the exact
identity with both rho branches -- never the sqrt(W) shortcut, which is what broke intervention 5.
"""
import sys
import numpy as np
ROOT="~/Code/virtual-cell-challenge-2026"
sys.path.insert(0,f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.metrics.delta import _across_pert_budget, _row_weights
from cell_eval2.prep import bulk_lognorm_means
from decompose import terms
from msekit import TS, group_stats
from mech_miller import donor_key

bs=np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve2.npz")
X0=float(bs['X']); sD=float(bs['sumD']); Kc=float(bs['Kc']); K=float(bs['K'])
BA,SP=0.987862879,0.950162566
GAM=[0.0,0.25,0.5,0.75,1.0,1.5,2.0]

# ---------- challenge side: A'(g), B'(g) exact ----------
z=np.load(f"{ROOT}/experiments/096-mse-floor/work/a10-stats.npz",allow_pickle=True)
ctrl=np.load(f"{ROOT}/experiments/096-mse-floor/work/ctrl-stats.npz",allow_pickle=True)
perts=z['perts'].astype(str); genes=z['genes'].astype(str); G=genes.size
gpos={g:i for i,g in enumerate(genes)}
Ag={g:[] for g in GAM}; Bg={g:[] for g in GAM}; surv={g:[] for g in GAM}
blank_tot=0
for c in "ABC":
    sel=np.flatnonzero(np.char.startswith(perts,f"{c}|"))
    tgt=np.asarray([p.split("|",1)[1] for p in perts[sel]])
    bctrl=bulk_lognorm_means(ctrl[f'sum_{c}'][None,:].astype(np.float64),TS)[0]
    cols={k:gpos[t] for k,t in enumerate(tgt) if t in gpos}
    w=_row_weights(sel.size,G,cols)
    pm=bulk_lognorm_means(z['sums'][sel].astype(np.float64),TS)
    dh=pm-bctrl[None,:]
    for k,col in cols.items(): dh[k,col]=0.0
    dk=donor_key(genes,list(tgt)); blank=(dk.sum(1)==0)
    blank_tot+=int(blank.sum())
    rank=np.argsort(np.argsort(-dk,axis=1),axis=1)
    for g in GAM:
        a=np.exp(-g*rank/2000.0); a[blank]=1.0
        D=a*dh
        aa=np.einsum("ij,ij->i",D,D)
        Ag[g].append(float((aa*w).sum()))
        Bg[g].append(_across_pert_budget(bctrl[None,:]+D,cols,w))
        surv[g].append(float(np.mean(np.linalg.norm(D,axis=1)/np.maximum(np.linalg.norm(dh,axis=1),1e-12))))
A={g:float(np.mean(Ag[g])) for g in GAM}; B={g:float(np.mean(Bg[g])) for g in GAM}
SV={g:float(np.mean(surv[g])) for g in GAM}
print(f"challenge panel: donor key blank on {blank_tot}/900 (context,target) rows")

# ---------- harness side: X'(g)/X multiplier ----------
mult={g:[] for g in GAM}
for real,pred,nm in [(f"{ROOT}/data/harness/h1-test.h5ad",
  f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad","h1"),
 (f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
  f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad","orion")]:
    rp,rs,rc,rj,rg=group_stats(real); pp,ps,pc,pj,pg=group_stats(pred)
    t=terms(pp,ps,pc,pj,rp,rs,rc,rj,rg); dh,d,w=t["dh"],t["d"],t["w"]
    for k,c in t["cols"].items(): dh[k,c]=0.0; d[k,c]=0.0
    tg=[str(x) for x in pp[t["keep"]]]
    dk=donor_key(rg,tg); blank=(dk.sum(1)==0)
    rank=np.argsort(np.argsort(-dk,axis=1),axis=1)
    X_un=float((np.einsum("ij,ij->i",dh,d)*w).sum())
    for g in GAM:
        a=np.exp(-g*rank/2000.0); a[blank]=1.0
        mult[g].append(float((np.einsum("ij,ij->i",a*dh,d)*w).sum())/X_un)
print(f"{'gamma':>6}{'surv amp':>10}{'A''':>10}{'B''':>10}{'X''/X h1':>10}{'orion':>8}"
      f"{'rule-1':>8}{'s*':>8}{'eff_scale':>11}{'raw mse':>10}{'score':>8}{'overall':>9}")
best=None
for g in GAM:
    m1=min(mult[g])
    Xp=X0*m1
    cl=K
    ss=np.linspace(0.05,6.0,4000)
    vals=1+(ss**2*A[g]-2*ss*Xp+Kc-np.minimum(cl,ss**2*B[g]))/sD
    i=int(np.argmin(vals)); s=ss[i]; v=vals[i]
    sc=min(max((BA-v)/SP,0),1)
    print(f"{g:>6.2f}{SV[g]:>10.4f}{A[g]:>10.5f}{B[g]:>10.5f}{mult[g][0]:>10.3f}{mult[g][1]:>8.3f}"
          f"{m1:>8.3f}{s:>8.4f}{s:>11.4f}{v:>10.4f}{sc:>8.4f}{sc/6:>+9.5f}")
    if best is None or v<best[1]: best=(g,v,s,sc)
print(f"\nBEST: gamma={best[0]}, effect_scale={best[2]:.4f}, raw mse {best[1]:.4f}, "
      f"score {best[3]:.4f}, overall {best[3]/6:+.5f}")
print(f"  (a12 unshrunk reference: raw 0.8882, score {min(max((BA-0.8882)/SP,0),1):.4f}, "
      f"overall {min(max((BA-0.8882)/SP,0),1)/6:+.5f})")
