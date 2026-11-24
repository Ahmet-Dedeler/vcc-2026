"""How much sampling correction can the emitter legitimately claim, per (spread, phi)?

`claim = sum_p w_p min(C_hat_p, C_p)` is the ONLY term in the numerator that emission moves --
the quasi emitter repairs each gene's panel total exactly, so the profile, and therefore `A`,
`B` and `X`, are invariant. #247 caps each row at the reference's own `C_p`, measured here at
14.54 / 14.77 / 16.37 (400 cells, one construct, per context).
"""
import sys, time
import numpy as np
from scipy.sparse import csr_matrix
ROOT="~/Code/virtual-cell-challenge-2026"
sys.path.insert(0,f"{ROOT}/src"); sys.path.insert(0,f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.moments import jackknife_correction
import mech_alpha_chal as M

D = M.load(f"{ROOT}/work" if False else f"{ROOT}/experiments/096-mse-floor/work/v11-stats.npz")
cpn = dict(np.load(f"{ROOT}/experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
GRID = [(0.25,0.0),(0.50,0.0),(1.00,0.0),(0.25,4.0),(1.00,4.0),(0.15,16.0),(0.25,16.0),(1.00,16.0),(1.00,40.0)]
NSUB = 20
print(f"{'spread':>7}{'phi':>6}{'mean C_hat':>12}{'C_p':>8}{'C_hat/C_p':>11}{'claim':>10}{'claim/sumD':>12}")
sD = float(np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve.npz")["sumD"])
for spread, phi in GRID:
    rows=[]
    for c in M.CTX:
        d=D[c]; idx=np.linspace(0,len(d['tgt'])-1,NSUB).astype(int)
        vals=[]
        for j in range(idx.size):
            prof = d['dh'][idx[j]]*0 + np.expm1(d['bctrl']+d['dh'][idx[j]])
            np.clip(prof,0,None,out=prof); prof/=prof.sum()
            libs=np.full(M.NCELL, d['S'][idx[j]]/M.NCELL)
            v,col,ip=M.emit_counts(prof,libs,0.0,np.random.default_rng(2000+j),mode="quasi",
                                   spread=spread, overdispersion=phi)
            X=csr_matrix((v.astype(np.float64),col.astype(np.int64),ip.astype(np.int64)),shape=(M.NCELL,d['G']))
            vals.append(jackknife_correction(X,np.zeros(M.NCELL,np.intp),1,M.TS,chunk=64)[0])
        cp=float(cpn[c].mean())
        rows.append((float(np.mean(vals)), cp, float(np.mean(np.minimum(vals,cp)))))
    ch=float(np.mean([r[0] for r in rows])); cp=float(np.mean([r[1] for r in rows]))
    claim=float(np.mean([300*r[2]/18532 for r in rows]))
    print(f"{spread:>7.2f}{phi:>6.1f}{ch:>12.4f}{cp:>8.4f}{ch/cp:>11.4f}{claim:>10.6f}{claim/sD:>12.4f}", flush=True)
