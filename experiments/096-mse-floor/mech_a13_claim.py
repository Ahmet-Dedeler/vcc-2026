"""Close the one omission every route shares: does `claim` move under shrinkage?

`claim = sum_p w_p min(C_hat_p, C_p)` is a functional of the EMITTED cells, and s* = sqrt(claim/B')
sits on it. All three routes held it at the unshrunk value. Re-emit the recommended arm and measure.
"""
import importlib.util as ilu, sys
import numpy as np
from scipy.sparse import csr_matrix
ROOT="~/Code/virtual-cell-challenge-2026"
sys.path.insert(0,f"{ROOT}/src"); sys.path.insert(0,f"{ROOT}/experiments/096-mse-floor")
_s=ilu.spec_from_file_location("mkpm",f"{ROOT}/scripts/make_submission_pm_cd4.py")
MK=ilu.module_from_spec(_s); sys.modules['mkpm']=MK; _s.loader.exec_module(MK)
from vcc2026.multi_donor import MultiDonorModel
from vcc2026.donor_transfer import Donor, emit_counts
from cell_eval2.moments import jackknife_correction
from cell_eval2.prep import bulk_lognorm_means
from msekit import TS
PROD=MK.PRODUCTION_DONORS
z=np.load(f"{ROOT}/experiments/096-mse-floor/work/a10-stats.npz",allow_pickle=True)
genes=z['genes'].astype(str); G=genes.size; gpos={g:i for i,g in enumerate(genes)}
perts=z['perts'].astype(str)
ctrl=np.load(f"{ROOT}/experiments/096-mse-floor/work/ctrl-stats.npz",allow_pickle=True)
cpn=dict(np.load(f"{ROOT}/experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
mdl=MultiDonorModel({n:Donor.load(f"{ROOT}/data/donors/{n}") for n in PROD},genes,gene_norm=False,eb=False)
print(f"{'arm':<34}{'C_hat':>9}{'/C_p':>7}{'claim':>9}{'vs cap':>9}")
for gg,gp,esc,lbl in ((0.0,0.0,0.9375,"a12 (no shrink)"),
                      (0.25,0.0,1.4068,"a13 gene-only  g=0.25"),
                      (0.25,1.0,2.4809,"a13 RECOMMENDED g=0.25,p=1.0"),
                      (1.0,1.0,3.1766,"a13 aggressive  g=1.0,p=1.0")):
    rows=[]
    for c in "ABC":
        sel=np.flatnonzero(np.char.startswith(perts,f"{c}|")); tgt=[p.split("|",1)[1] for p in perts[sel]]
        bctrl=bulk_lognorm_means(ctrl[f'sum_{c}'][None,:].astype(np.float64),TS)[0]
        dh=bulk_lognorm_means(z['sums'][sel].astype(np.float64),TS)-bctrl[None,:]
        KK=np.vstack([mdl.donor_agreement(str(t)) for t in tgt])
        zs,car=MK.donor_cosine_key(mdl,[str(t) for t in tgt])
        zp=MK.pert_shrink_vector(zs,car,gp) if gp>0 else np.ones(len(tgt))
        idx=np.linspace(0,len(tgt)-1,14).astype(int); vals=[]
        for j in idx:
            al=MK.agreement_shrink_vector(KK[j],gg,protect=gpos.get(str(tgt[j])))
            eff=esc*zp[j]*al*dh[j]
            prof=np.clip(np.expm1(bctrl+eff),0,None); prof/=prof.sum()
            top=np.argpartition(-np.abs(eff),8000)[:8000]
            sp=np.full(G,1.0); ph=np.full(G,4.0); sp[top]=0.25
            v,col,ip=emit_counts(prof,np.full(400,z['sums'][sel][j].sum()/400),0.0,
                                 np.random.default_rng(12000+j),mode="quasi",spread=sp,overdispersion=ph)
            M=csr_matrix((v.astype(np.float64),col.astype(np.int64),ip.astype(np.int64)),shape=(400,G))
            vals.append(jackknife_correction(M,np.zeros(400,np.intp),1,TS,chunk=64)[0])
        cp=float(cpn[c].mean()); rows.append((float(np.mean(vals)),cp,float(np.mean(np.minimum(vals,cp)))))
    ch=float(np.mean([r[0] for r in rows])); cp=float(np.mean([r[1] for r in rows]))
    claim=float(np.mean([300*r[2]/18532 for r in rows]))
    print(f"{lbl:<34}{ch:>9.2f}{ch/cp:>7.2f}{claim:>9.5f}{claim/0.246529:>9.4f}")
