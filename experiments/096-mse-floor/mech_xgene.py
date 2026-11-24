"""Is a per-gene shrinkage BUILDABLE? Fit s_g on one cell line, test it on another.

Rule 4 (RESEARCH-METHOD): a coefficient fitted against the reference's own lfc is unbuildable
at submission time. The only leak-free fit is against a DIFFERENT dataset. orion-hct116 (HCT116,
132 perturbations) and h1-test (H1 hESC, 85 perturbations) are two such datasets, scored against
the same production model, so `s_g` fitted on one and applied to the other is exactly the
construction a submission could ship.
"""
import sys
import numpy as np
ROOT="~/Code/virtual-cell-challenge-2026"
sys.path.insert(0,f"{ROOT}/experiments/096-mse-floor")
from decompose import terms
from msekit import group_stats

def load(real,pred):
    rp,rs,rc,rj,rg=group_stats(real); pp,ps,pc,pj,pg=group_stats(pred)
    t=terms(pp,ps,pc,pj,rp,rs,rc,rj,rg); dh,d,w=t["dh"],t["d"],t["w"]
    for k,c in t["cols"].items(): dh[k,c]=0.0; d[k,c]=0.0
    return dh,d,w,rg

def W_of(dh,d,w):
    n=float((np.einsum("ij,ij->i",dh,d)*w).sum()); q=float((np.einsum("ij,ij->i",dh,dh)*w).sum())
    return n*n/q if q>0 else 0.0

H=(f"{ROOT}/data/harness/h1-test.h5ad",
   f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad")
O=(f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
   f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad")
dhH,dH,wH,gH = load(*H); dhO,dO,wO,gO = load(*O)
pos={g:i for i,g in enumerate(gO)}
col=np.array([pos.get(g,-1) for g in gH]); keep=col>=0
print(f"gene axes: h1-test {gH.size}, orion {gO.size}, shared {keep.sum()}")
dhH,dHk,wHk = dhH[:,keep], dH[:,keep], wH
dhO2, dO2 = dhO[:,col[keep]], dO[:,col[keep]]
for src,dst,sn,dn,dhs,ds,ws,dhd,dd,wd in [("orion","h1-test",0,0,dhO2,dO2,wO,dhH,dHk,wHk),
                                          ("h1-test","orion",0,0,dhH,dHk,wHk,dhO2,dO2,wO)]:
    U=np.einsum("ij,ij->j",dhs*ws[:,None],dhs); V=np.einsum("ij,ij->j",dhs*ws[:,None],ds)
    A=float((np.einsum("ij,ij->i",dhs,dhs)*ws).sum()); X=float((np.einsum("ij,ij->i",dhs,ds)*ws).sum())
    g=X/A; ok=U>0; sg=np.full(U.size,g); sg[ok]=V[ok]/U[ok]
    base=W_of(dhd,dd,wd)
    print(f"\n  fit s_g on {src}  ->  evaluate on {dst}   (W_global on {dst} = {base:.6g}; "
          f"global coeff from {src} = {g:.4f})")
    print(f"    {'lambda':>10}{'W':>14}{'sqrt(W/W_glob)':>17}")
    for l in (0.0,1,3,10,30,100,1000,1e6):
        s=(sg+l*g)/(1.0+l)
        v=W_of(s[None,:]*dhd,dd,wd)
        print(f"    {l:>10.4g}{v:>14.6g}{np.sqrt(v/base):>17.3f}")
