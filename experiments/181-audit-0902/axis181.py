"""181 — recompute 175's per-donor top PC (the 'platform axis') and save the loadings.
KA rows (175 addendum): HEK293T top-1 variance share 19.7%, K562 5.6%, HCT116 8.6%;
cos(HEK top1, HCT top1) = 0.872, both vs K562 0.23-0.31."""
import numpy as np, sys, itertools
sys.path.insert(0,'experiments/153-de-bench'); sys.path.insert(0,'experiments/172-sign-confidence-key')
from debench import Surface
from donors172 import DonorHandle

S = Surface('experiments/164-callset-artefact/de164.a15.npz')
genes = np.asarray([str(g) for g in S.genes]); targets=[str(p) for p in S.perts]
rng = np.random.default_rng(0)
PC1={}; VAR={}
for n in ('k562_gwps','orion_hct116','orion_hek293t'):
    h = DonorHandle(n)
    tcols,dcols = h.project(genes)
    sub = np.sort(rng.choice(h.n_perts, min(4000,h.n_perts), replace=False))
    M = np.zeros((sub.size,len(genes)))
    for r,i in enumerate(sub): M[r,tcols] = (np.asarray(h.delta[i],dtype=np.float64)-h.m)[dcols]
    U,s_,Vt = np.linalg.svd(M, full_matrices=False)
    PC1[n]=Vt[0].copy(); VAR[n]=(s_**2/np.sum(s_**2))[:3]
    print(f'{n:<15} top-3 variance share: {VAR[n][0]:.4f} {VAR[n][1]:.4f} {VAR[n][2]:.4f}', flush=True)
    del M,U,Vt,h
print()
for a,b in itertools.combinations(PC1,2):
    print(f'  cos(top1 {a}, top1 {b}) = {abs(float(PC1[a]@PC1[b])):.3f}')
np.savez('experiments/181-audit-0902/work/pc1.npz', genes=genes,
         **{n:PC1[n] for n in PC1}, **{f'var_{n}':VAR[n] for n in VAR})
v = PC1['orion_hek293t']
v = v * (-1.0 if v[np.argmax(np.abs(v))] > 0 else 1.0)
o = np.argsort(v)
print('\ntop-30 most NEGATIVE loadings (175: AUTS2 GPC6 NRXN3 PRKG1 WWOX RAD51B RORA ERC2 CACNA2D3 PARD3B ADGRL3 FAT3 MAGI1 PTPRG ...):')
print('  ', ' '.join(genes[o[:30]]))
print('top-30 most POSITIVE loadings:')
print('  ', ' '.join(genes[o[-30:]]))
