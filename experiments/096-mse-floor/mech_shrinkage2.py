"""Which shrinkage families are BUILDABLE, and what each is worth (h1-test, truth known).

Everything is reported as `sqrt(W/W_global)`, the multiplier on the mse gain, because
`MSE_min = 1 + (Kc - 2 sqrt(claim * W))/sumD`.

Families:
  global            one scalar                                        buildable (rho cliff)
  per-gene ridge    s_g = (V_g/U_g + lam*s_glob)/(1+lam), lam swept    needs truth -> CV'd here
  per-pert oracle   s_p = <dh_p,d_p>/||dh_p||^2                        NOT buildable
  norm-equalise     s_p proportional to 1/||dh_p||                     BUILDABLE (own file only)
  norm-power        s_p proportional to ||dh_p||^(-t), t swept         BUILDABLE
"""
import sys
import numpy as np
ROOT="~/Code/virtual-cell-challenge-2026"
sys.path.insert(0,f"{ROOT}/experiments/096-mse-floor")
from decompose import terms
from msekit import group_stats

def load(real, pred):
    rp,rs,rc,rj,rg = group_stats(real); pp,ps,pc,pj,pg = group_stats(pred)
    t = terms(pp,ps,pc,pj,rp,rs,rc,rj,rg)
    dh,d,w = t["dh"], t["d"], t["w"]
    for k,c in t["cols"].items(): dh[k,c]=0.0; d[k,c]=0.0
    return t,dh,d,w

def W_of(dh,d,w):
    num=float((np.einsum("ij,ij->i",dh,d)*w).sum()); den=float((np.einsum("ij,ij->i",dh,dh)*w).sum())
    return num*num/den if den>0 else 0.0

for real,pred,name in [(f"{ROOT}/data/harness/h1-test.h5ad",
   f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad","h1-test"),
  (f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
   f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad","orion-thinned")]:
    t,dh,d,w = load(real,pred)
    P,G = dh.shape
    base = W_of(dh,d,w)
    print(f"\n===== {name}: {P} perturbations, {G} genes.  W_global = {base:.6g} =====")
    # ---- per-gene ridge, honestly cross-validated over perturbations ----
    rng=np.random.default_rng(0); lams=[0.0,0.3,1,3,10,30,100,300,1000,1e4,1e6]
    acc={l:[] for l in lams}
    for rep in range(30):
        perm=rng.permutation(P); fit,ev = perm[:P//2], perm[P//2:]
        U=np.einsum("ij,ij->j",dh[fit]*w[fit,None],dh[fit]); V=np.einsum("ij,ij->j",dh[fit]*w[fit,None],d[fit])
        A=float((np.einsum("ij,ij->i",dh[fit],dh[fit])*w[fit]).sum()); X=float((np.einsum("ij,ij->i",dh[fit],d[fit])*w[fit]).sum())
        g=X/A; ok=U>0; sg=np.full(G,g); sg[ok]=V[ok]/U[ok]
        for l in lams:
            s=(sg+l*g)/(1.0+l)
            acc[l].append(W_of(s[None,:]*dh[ev], d[ev], w[ev]))
    print("  per-gene ridge toward the global coefficient, 30 random half-splits:")
    print(f"    {'lambda':>10}{'W(out-of-sample)':>20}{'sqrt(W/W_glob)':>17}")
    for l in lams:
        v=float(np.mean(acc[l])); print(f"    {l:>10.4g}{v:>20.6g}{np.sqrt(v/base):>17.3f}")
    U=np.einsum("ij,ij->j",dh*w[:,None],dh); V=np.einsum("ij,ij->j",dh*w[:,None],d); ok=U>0
    print(f"    {'ORACLE':>10}{float((V[ok]**2/U[ok]).sum()):>20.6g}"
          f"{np.sqrt(float((V[ok]**2/U[ok]).sum())/base):>17.3f}")
    # ---- per-perturbation families ----
    nr=np.sqrt(np.einsum("ij,ij->i",dh,dh)); xp=np.einsum("ij,ij->i",dh,d)
    print("  per-perturbation families:")
    orc=float((w[nr>0]*(xp[nr>0]/nr[nr>0])**2).sum())
    print(f"    {'oracle s_p = <dh,d>/||dh||^2':<34}{orc:>14.6g}{np.sqrt(orc/base):>10.3f}x   NOT buildable")
    for tpow in (0.0,0.25,0.5,0.75,1.0,1.25):
        s=np.where(nr>0, nr**(-tpow), 0.0)
        v=W_of(s[:,None]*dh,d,w)
        print(f"    {'||dh_p||^(-%.2f)'%tpow:<34}{v:>14.6g}{np.sqrt(v/base):>10.3f}x   BUILDABLE")
