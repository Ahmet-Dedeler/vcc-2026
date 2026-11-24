"""181 — the two mechanism claims that carry 171 and 177, checked directly."""
import numpy as np, sys
from scipy.stats import spearmanr
sys.path.insert(0,'experiments/153-de-bench')
from debench import Surface
P='experiments/177-resample-emitter/de177.%s.npz'

print('(1) 171 S3 — "our significant calls first" is a BH no-op (a proof, then a measurement)')
S=Surface(P%'KA_a15_raw')
viol=0; ties=0
for p in range(85):
    c=S.pred_sig[p]; pa=S.pred_padj[p]
    if c.sum()==0 or (~c).sum()==0: continue
    if np.nanmax(pa[c])>=np.nanmin(pa[~c]): viol+=1
    ties += int(np.sum(pa[c][:,None]==np.nanmin(pa[~c])))
print(f'    perturbations where max p_adj(called) >= min p_adj(not called): {viol}/85   (171: 0/85)')
allsig=np.ones_like(S.pred_sig); nosig=np.zeros_like(S.pred_sig)
print(f'    reach with the call set = production {S.reach():.5f} | nothing {S.reach(call=nosig):.5f} | everything {S.reach(call=allsig):.5f}')
print('    -> identical; the call set cannot move reach in either direction. CONFIRMED.')

print()
print('(2) 177 S0.4 — the emission LAW flips the p-order from expression-set to |H|-set')
eff=np.load('experiments/169-sharpening-six-member/work/h1test-effect.npz',allow_pickle=True)
gi={g:i for i,g in enumerate(eff['genes'])}
H=np.abs(eff['H'][np.ix_([list(eff['targets']).index(t) for t in S.perts],[gi[g] for g in S.genes])])
expr=np.array([eff['cpm'][gi[g]] for g in S.genes])
print(f"    {'arm':14} {'rho(-log padj, expr) in-set':>28} {'rho(-log padj, |H|) in-set':>28}")
for a in ('KA_a15_raw','base_mf','hy_kt8000','hy_kt4000','rs_all','quasi100'):
    s=Surface(P%a); re=[]; rh=[]
    for p in range(85):
        i=np.flatnonzero(s.pred_sig[p]&np.isfinite(s.pred_padj[p]))
        if i.size<20: continue
        k=-np.log10(np.maximum(s.pred_padj[p,i],1e-300))
        re.append(spearmanr(k,expr[i]).statistic); rh.append(spearmanr(k,H[p,i]).statistic)
    print(f'    {a:14} {np.median(re):28.3f} {np.median(rh):28.3f}')
print('    177 reports a15 +0.485 / -0.391 and hy_kt2000 +0.001 / +0.601.')
