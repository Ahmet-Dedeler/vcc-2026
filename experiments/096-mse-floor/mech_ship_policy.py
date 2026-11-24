"""Can `claim` be raised to the #247 cap WITHOUT cutting the called set?

`C_hat_p` sums over all 18,533 genes; the call set is ~4,750. So over-dispersing only the genes
that are NOT called raises the claimable correction while leaving the called genes' p-values
where they are. That is exactly what `emission_spread_vector`'s `keep_top` / `keep_spread` /
`tail_spread` already do, and `emit_counts` takes a per-gene `overdispersion` too.

Policy: top-K genes by |predicted log-effect| at (keep_spread, keep_phi); the rest at
(tail_spread, tail_phi).
"""
import sys
import numpy as np
from scipy.sparse import csr_matrix
ROOT="~/Code/virtual-cell-challenge-2026"
sys.path.insert(0,f"{ROOT}/src"); sys.path.insert(0,f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.moments import jackknife_correction
import mech_alpha_chal as M
D=M.load(f"{ROOT}/experiments/096-mse-floor/work/v11-stats.npz")
cpn=dict(np.load(f"{ROOT}/experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
bs=np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve.npz")
sD=float(bs["sumD"]); A=float(bs["v11_A"]); B=float(bs["v11_B"]); Kc=float(bs["v11_Kc"]); X=float(bs["X"])
POL=[(8000,0.15,4.0,1.00,4.0),(8000,0.15,16.0,1.00,16.0),(5000,0.15,4.0,1.00,4.0),(8000,0.15,4.0,0.60,4.0),(10000,0.15,4.0,1.00,4.0),(8000,0.25,4.0,1.00,4.0),(12000,0.15,4.0,1.00,4.0),
     (8000,0.15,16.0,1.00,4.0)]
NSUB=16
print(f"{'keepK':>7}{'keep(sp,phi)':>15}{'tail(sp,phi)':>15}{'C_hat':>9}{'/C_p':>7}{'claim':>9}"
      f"{'s*':>8}{'ALPHA*':>9}{'mse@s*':>9}{'score':>8}{'overall':>9}")
for K,ks,kp,ts,tp in POL:
    rows=[]
    for c in M.CTX:
        d=D[c]; idx=np.linspace(0,len(d['tgt'])-1,NSUB).astype(int); vals=[]
        for j in range(idx.size):
            e=d['dh'][idx[j]]
            prof=np.expm1(d['bctrl']+e); np.clip(prof,0,None,out=prof); prof/=prof.sum()
            top=np.argpartition(-np.abs(e),K)[:K]
            spread=np.full(e.size,ts); spread[top]=ks
            phi=np.full(e.size,tp); phi[top]=kp
            libs=np.full(M.NCELL,d['S'][idx[j]]/M.NCELL)
            v,col,ip=M.emit_counts(prof,libs,0.0,np.random.default_rng(4000+j),mode="quasi",
                                   spread=spread,overdispersion=phi)
            Xm=csr_matrix((v.astype(np.float64),col.astype(np.int64),ip.astype(np.int64)),shape=(M.NCELL,d['G']))
            vals.append(jackknife_correction(Xm,np.zeros(M.NCELL,np.intp),1,M.TS,chunk=64)[0])
        cp=float(cpn[c].mean()); rows.append((float(np.mean(vals)),cp,float(np.mean(np.minimum(vals,cp)))))
    ch=float(np.mean([r[0] for r in rows])); cp=float(np.mean([r[1] for r in rows]))
    claim=float(np.mean([300*r[2]/18532 for r in rows]))
    s=np.sqrt(claim/B); mse=1+(s*s*(A-B)+Kc-2*s*X)/sD
    sc=min(max((0.987862879-mse)/0.950162566,0),1)
    print(f"{K:>7}{f'({ks},{kp})':>15}{f'({ts},{tp})':>15}{ch:>9.3f}{ch/cp:>7.3f}{claim:>9.5f}"
          f"{s:>8.4f}{1.0796*s:>9.4f}{mse:>9.4f}{sc:>8.4f}{sc/6:>9.5f}", flush=True)
