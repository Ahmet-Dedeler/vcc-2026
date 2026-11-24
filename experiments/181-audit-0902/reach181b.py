"""181 — the test 177 did not run: is the hy_kt8000 reach gain carried by the small-n_real
stratum the challenge does not have?  171 S7a(i) killed --keep-top 8300 on exactly this
(gain 3x smaller on n_real >= 20).  172 S7.2: 94.4% of h1-test reach lives in 32 targets
with n_real <= 20.  The challenge's n_real is ~240 for every perturbation.
"""
import numpy as np, sys
sys.path.insert(0,'experiments/153-de-bench')
from debench import Surface
DE='experiments/177-resample-emitter/de177.%s.npz'

def reach_per_pert(s, key=None):
    out=np.full(len(s.perts), np.nan)
    for p in range(len(s.perts)):
        pool=np.flatnonzero(s.real_sig[p])
        if pool.size==0: continue
        order=s._order(p,pool,key,None); pool=pool[order]; pool=pool[s.adj[p,pool]]
        if pool.size==0: out[p]=0.0; continue
        ok=np.sign(s.pred_lfc[p,pool])==np.sign(s.real_lfc[p,pool])
        pur=np.cumsum(ok)/np.arange(1,pool.size+1)
        deep=np.flatnonzero(pur>=0.9); k=int(deep[-1]+1) if deep.size else 0
        out[p]=k/int(s.real_sig[p].sum())
    return out

arms={a:Surface(DE%a) for a in ['base_mf','base_mf_b','hy_kt8000','hy_kt4000','rs_all','KA_a15_raw']}
R={a:reach_per_pert(s) for a,s in arms.items()}
nr=np.array([int(arms['base_mf'].real_sig[p].sum()) for p in range(85)])
base=np.nanmean([R['base_mf'],R['base_mf_b']],axis=0)   # twin mean, as 177 does

print(f'n_real: min {nr.min()} med {int(np.median(nr))} max {nr.max()};  '
      f'perts with n_real<=13: {(nr<=13).sum()},  <=20: {(nr<=20).sum()},  >=20: {(nr>=20).sum()}')
print()
print('Delta reach vs the base TWIN MEAN, by n_real stratum (mean over the perturbations in it,')
print('renormalised to the full 85-pert panel mean the metric actually reports):')
print(f"{'arm':11} {'full':>9} {'n_real<=13':>11} {'14-800':>9} {'>800':>9} {'n_real>=20':>11} {'n_real>=100':>12}")
strat = {'<=13': nr<=13, '14-800': (nr>13)&(nr<=800), '>800': nr>800, '>=20': nr>=20, '>=100': nr>=100}
for a in ['hy_kt8000','hy_kt4000','rs_all','KA_a15_raw']:
    d = R[a]-base
    row=[np.nanmean(d)]
    for k in ['<=13','14-800','>800','>=20','>=100']:
        m=strat[k]; row.append(np.nanmean(d[m]))
    print(f'{a:11} '+' '.join(f'{v:+9.5f}' if i==0 else f'{v:+11.5f}' if i in(1,4) else f'{v:+9.5f}' for i,v in enumerate(row[:4]))+f' {row[4]:+11.5f} {row[5]:+12.5f}')

print()
print('Share of the FULL-panel delta contributed by each stratum (sum of d over stratum / sum over all):')
for a in ['hy_kt8000','hy_kt4000','rs_all']:
    d=R[a]-base; tot=np.nansum(d)
    parts={k: np.nansum(d[m])/tot for k,m in strat.items() if k in ('<=13','14-800','>800')}
    print(f"  {a:11} n_real<=13 {parts['<=13']:6.1%}   14-800 {parts['14-800']:6.1%}   >800 {parts['>800']:6.1%}   (n perts {int((nr<=13).sum())}/{int(((nr>13)&(nr<=800)).sum())}/{int((nr>800).sum())})")

print()
print('Top-5 |delta| perturbations, hy_kt8000 (171 S7a(i) found 69% of its arm in five perts):')
d=R['hy_kt8000']-base; idx=np.argsort(-np.abs(d))[:5]
tot=np.nansum(np.abs(d))
for i in idx:
    print(f'   {arms["base_mf"].perts[i]:10} n_real {nr[i]:5d}  d {d[i]:+.4f}  {abs(d[i])/tot:5.1%} of the |sum|')
print(f'   top-5 share of the signed panel sum: {np.nansum(d[idx])/np.nansum(d):.1%}')
print(f'   delta with the top-1 removed: {np.nanmean(np.delete(d, idx[0])):+.5f}  (full {np.nanmean(d):+.5f})')

print()
print('--- the same question for the |H| COUNTERFACTUAL, which is what 177 says it realises ---')
eff=np.load('experiments/169-sharpening-six-member/work/h1test-effect.npz',allow_pickle=True)
gi={g:i for i,g in enumerate(eff['genes'])}; ti={t:i for i,t in enumerate(eff['targets'])}
S=arms['KA_a15_raw']
H=np.abs(eff['H'][np.ix_([ti[t] for t in S.perts],[gi[g] for g in S.genes])])
rk=reach_per_pert(S, key=H); r0=reach_per_pert(S)
d=rk-r0
print(f'  |H| counterfactual full panel: {np.nanmean(d):+.5f}  (0.23533-0.18283 = +0.05250)')
for k,m in strat.items():
    print(f'    stratum {k:7} n={int(m.sum()):3d}  mean d {np.nanmean(d[m]):+.5f}  share of sum {np.nansum(d[m])/np.nansum(d):6.1%}')
