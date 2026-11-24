"""181 — is the board-jac -> recall inversion in CURRENT_STATE valid, and does the board's
own `reach` level constrain the challenge n_real independently?

(1) jac is the MEAN of per-perturbation Jaccards.  Inverting the MEAN with mean n_real /
    mean n_pred is a Jensen approximation.  Test it on h1-test, where every quantity is known.
(2) reach = mean_p k*/n_real_p.  If the challenge's n_real were ~240 (let alone 2,670) for
    every perturbation, reach would collapse.  Measure the reach the h1-test surface gives at
    matched power and compare with the board's 0.1742.
"""
import numpy as np, sys
sys.path.insert(0,'experiments/153-de-bench')
from debench import Surface
DE='experiments/177-resample-emitter/de177.%s.npz'
S=Surface(DE%'KA_a15_raw')

nr=np.array([int(S.real_sig[p].sum()) for p in range(85)])
npd=np.array([int((S.pred_sig[p]&S.adj[p]).sum()) for p in range(85)])
rec=[]; jacs=[]
for p in range(85):
    a=S.real_sig[p]; b=S.pred_sig[p]; u=int((a|b).sum()); i=int((a&b).sum())
    jacs.append(1.0 if u==0 else i/u)
    if a.sum(): rec.append(i/int(a.sum()))
jacs=np.array(jacs); rec=np.array(rec)
J=float(jacs.mean()); R=float(rec.mean())
print(f'h1-test truth: mean jac {J:.5f}  mean recall {R:.4f}  median n_real {int(np.median(nr))}  mean n_real {nr.mean():.0f}  median n_pred {int(np.median(npd))}')

def invert(jac, n_real, n_pred):
    return jac*(n_real+n_pred)/(n_real*(1+jac))
print()
print('(1) CURRENT_STATE inversion  r = jac(n_real+n_pred) / (n_real(1+jac))')
for lab, nrv in [('median n_real', np.median(nr)), ('mean n_real', nr.mean())]:
    r = invert(J, nrv, np.median(npd))
    print(f'    using {lab:14} = {nrv:7.0f} -> implied recall {r:.4f}   (true mean recall {R:.4f}, error {r-R:+.4f})')
print(f'    the inversion is off by {invert(J,np.median(nr),np.median(npd))/R:.2f}x using the median and '
      f'{invert(J,nr.mean(),np.median(npd))/R:.2f}x using the mean.')
print('    -> the map jac -> recall is NOT invertible from panel means when n_real is heavy-tailed.')

print()
print('    the same inversion applied to the board (jac 0.0301, n_pred 4424):')
for n in (191,240,400,800,2670,4716):
    print(f'      n_real {n:5d} -> implied recall {invert(0.0301, n, 4424):.4f}   '
          f'(CURRENT_STATE quotes .57-.71 at 191-240 and .078 at 2670)')

print()
print('(2) reach vs power — the board reads 0.1742 raw for a15/a23')
print('    debench proportional power matching (157 S3b / 172 S7.3), production order:')
for f in (0.5,1.0,1.5,2.0,3.0,3.5,5.0,10.0):
    tgt=int(round(f*np.median(nr)))
    Sf=Surface(DE%'KA_a15_raw'); 
    # proportional: take each pert's top round(f*n_real_p) by the reference's own p
    rp=np.where(np.isfinite(np.load(DE%'KA_a15_raw')['real_p']), np.load(DE%'KA_a15_raw')['real_p'], np.inf)
    sig=np.zeros(rp.shape,bool)
    for p in range(85):
        k=int(round(f*nr[p]))
        if k<=0: continue
        k=min(k,rp.shape[1])
        sig[p,np.argpartition(rp[p],k-1)[:k]]=True
    Sf.real_sig=sig
    for p,j in enumerate(Sf.self_idx):
        if j>=0: Sf.real_sig[p,j]=False
    print(f'      f={f:4.1f}  median n_real {int(np.median(sig.sum(1))):5d}   reach {Sf.reach():.5f}')
print('    board raw reach = 0.17420 (a15) / 0.17416 (a23).  A uniform challenge n_real of 240,')
print('    let alone 2,670, is not compatible with that level on this surface.')
