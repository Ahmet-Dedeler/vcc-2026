"""The free lever: genes BELOW the DE metrics' CPM gate still feed `C_hat_p`.

`filter_gene_min_cpm_cell = 5.0` is DE-only (`distances.py:98`, `de_compute.py:414`): the four
DE members never see a gene whose group mean CPM is under 5, and `pds`/`nmae` are invariant to
emission anyway because the quasi emitter repairs each gene's panel total exactly. But
`expr_mse_unbiased_capped`'s `C_hat_p` is a delete-1 jackknife over ALL 18,533 genes.

So scatter can be manufactured in the sub-gate genes and claimed against #247's cap while every
DE member and `pds` stay bit-identical. This measures how much is available there.
"""
import sys
import numpy as np
from scipy.sparse import csr_matrix
ROOT="~/Code/virtual-cell-challenge-2026"
sys.path.insert(0,f"{ROOT}/src"); sys.path.insert(0,f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.moments import jackknife_correction
import mech_alpha_chal as M
D=M.load(f"{ROOT}/experiments/096-mse-floor/work/v11-stats.npz")
ctrl=np.load(f"{ROOT}/experiments/096-mse-floor/work/ctrl-stats.npz",allow_pickle=True)
cpn=dict(np.load(f"{ROOT}/experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
bs=np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve.npz")
sD=float(bs["sumD"]);A=float(bs["v11_A"]);B=float(bs["v11_B"]);Kc=float(bs["v11_Kc"]);X=float(bs["X"])
NSUB=16
for c in M.CTX:
    s=ctrl[f"sum_{c}"].astype(np.float64); cpm=1e6*s/s.sum()
    print(f"ctx {c}: genes with control CPM <=5: {(cpm<=5).sum()} of {cpm.size} "
          f"({100*(cpm<=5).mean():.1f} %), holding {100*s[cpm<=5].sum()/s.sum():.2f} % of counts; "
          f"<=2: {(cpm<=2).sum()} ({100*s[cpm<=2].sum()/s.sum():.2f} % of counts)")
POL=[("sub-gate CPM<=5 at (1.0, 0)",5.0,1.0,0.0),("sub-gate CPM<=5 at (1.0, 4)",5.0,1.0,4.0),
     ("sub-gate CPM<=5 at (1.0, 16)",5.0,1.0,16.0),("sub-gate CPM<=2 at (1.0, 16)",2.0,1.0,16.0),
     ("sub-gate CPM<=5 at (1.0, 64)",5.0,1.0,64.0)]
print(f"\nbase policy on the gated genes stays a10's (spread 0.15, phi 16)")
print(f"{'policy':>32}{'C_hat':>9}{'/C_p':>7}{'claim':>9}{'s*':>8}{'ALPHA*':>9}{'mse@s*':>9}{'score':>8}{'overall':>9}")
for name,thr,ts,tp in POL:
    rows=[]
    for c in M.CTX:
        d=D[c]; sc=ctrl[f"sum_{c}"].astype(np.float64); cpm=1e6*sc/sc.sum(); sub=cpm<=thr
        idx=np.linspace(0,len(d['tgt'])-1,NSUB).astype(int); vals=[]
        for j in range(idx.size):
            e=d['dh'][idx[j]]
            prof=np.expm1(d['bctrl']+e); np.clip(prof,0,None,out=prof); prof/=prof.sum()
            spread=np.full(e.size,0.15); spread[sub]=ts
            phi=np.full(e.size,16.0);   phi[sub]=tp
            libs=np.full(M.NCELL,d['S'][idx[j]]/M.NCELL)
            v,col,ip=M.emit_counts(prof,libs,0.0,np.random.default_rng(5000+j),mode="quasi",
                                   spread=spread,overdispersion=phi)
            Xm=csr_matrix((v.astype(np.float64),col.astype(np.int64),ip.astype(np.int64)),shape=(M.NCELL,d['G']))
            vals.append(jackknife_correction(Xm,np.zeros(M.NCELL,np.intp),1,M.TS,chunk=64)[0])
        cp=float(cpn[c].mean()); rows.append((float(np.mean(vals)),cp,float(np.mean(np.minimum(vals,cp)))))
    ch=float(np.mean([r[0] for r in rows])); cp=float(np.mean([r[1] for r in rows]))
    claim=float(np.mean([300*r[2]/18532 for r in rows]))
    s=np.sqrt(claim/B); mse=1+(s*s*(A-B)+Kc-2*s*X)/sD
    sco=min(max((0.987862879-mse)/0.950162566,0),1)
    print(f"{name:>32}{ch:>9.3f}{ch/cp:>7.3f}{claim:>9.5f}{s:>8.4f}{1.0796*s:>9.4f}{mse:>9.4f}{sco:>8.4f}{sco/6:>9.5f}",flush=True)
