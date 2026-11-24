"""181 — two mask questions the audit was asked to test.

(a) 175: 'the Orion pair adds ~nothing over K562 alone' is measured on the ALL-THREE-CARRIED
    mask, which is a subset of K562's 8,248-gene axis.  What happens on the pool pairs K562
    cannot reach at all?  If Orion signs those above chance the claim is mask-limited.
(b) 176: the CRISPR-KO +0.019 landmark gain was capped by 'landmarks are ~13% of pool pairs'.
    Are the 867 landmark genes OVER-represented in the reference's DE pool (which would raise
    the pooled effect above 13%)?
"""
import numpy as np, sys
sys.path.insert(0,'experiments/153-de-bench'); sys.path.insert(0,'experiments/172-sign-confidence-key')
from debench import Surface
from donors172 import DonorHandle

S = Surface('experiments/164-callset-artefact/de164.a15.npz')
pool = S.real_sig & S.adj; sT = np.sign(S.real_lfc)
genes = np.asarray([str(g) for g in S.genes]); targets=[str(p) for p in S.perts]
print(f'reference DE pool: {int(pool.sum())} (pert,gene) pairs over {len(targets)} targets x {len(genes)} genes')

D={}
for n in ('k562_gwps','orion_hct116','orion_hek293t'):
    h=DonorHandle(n); Sd,c = h.rows(targets,genes); D[n]=(Sd.astype(np.float64), c, set(map(str,h.genes)))
    print(f'  {n:<15} carries {int(c.sum())}/85 targets, {len(h.genes)} genes on its own axis')

K,Hc,He = D['k562_gwps'][0], D['orion_hct116'][0], D['orion_hek293t'][0]
nzK, nzH, nzE = np.sign(K)!=0, np.sign(Hc)!=0, np.sign(He)!=0
allc = np.all([D[n][1] for n in D],axis=0)
m3 = pool & allc[:,None] & nzK & nzH & nzE
print(f'\nKA (175): all-three mask n={int(m3.sum())} (175 reports 39,945)')
def a(k,m): return float(np.mean(np.sign(k)[m]==sT[m]))
print(f'  K562 {a(K,m3):.4f}  HCT116 {a(Hc,m3):.4f}  HEK293T {a(He,m3):.4f}  '
      f'orion-pair {a(Hc+He,m3):.4f}  equal {a(K+Hc+He,m3):.4f}   (175: .5670 .5375 .5308 .5425 .5697)')

print('\n(a) the pool pairs K562 CANNOT reach — where the Orion pair is the whole consensus')
onk = np.isin(genes, list(D['k562_gwps'][2]))
print(f'    genes on K562 axis: {int(onk.sum())}/{len(genes)}; pool pairs off it: {int((pool & ~onk[None,:]).sum())}')
strata = {
 'K562 nonzero + both orion nonzero (175 mask)': m3,
 'K562 ZERO/absent, both orion nonzero'       : pool & ~nzK & nzH & nzE,
 'K562 ZERO/absent, >=1 orion nonzero'        : pool & ~nzK & (nzH|nzE),
 'gene off K562 axis entirely, >=1 orion nz'  : pool & ~onk[None,:] & (nzH|nzE),
}
print(f"    {'stratum':46} {'n':>8} {'K562':>8} {'orion pair':>11} {'equal(3)':>9}")
for k,m in strata.items():
    if m.sum()==0: continue
    print(f'    {k:46} {int(m.sum()):8d} {a(K,m):8.4f} {a(Hc+He,m):11.4f} {a(K+Hc+He,m):9.4f}')

print('\n    what production actually delivers: the FULL pool, consensus over whichever donors carry')
mfull = pool & (nzK|nzH|nzE)
print(f'      full pool with >=1 donor  n={int(mfull.sum())}  equal-mean {a(K+Hc+He,mfull):.4f}   K562-only key {a(K,mfull):.4f}')
print(f'      K562 alone can sign only {int((pool&nzK).sum())} of {int(mfull.sum())} = {(pool&nzK).sum()/mfull.sum():.1%} of the signable pool')
# what the shipped H actually reads
print(f'      shipped emitted pred_lfc on the full pool: {float(np.mean(np.sign(S.pred_lfc)[pool]==sT[pool])):.4f}')

print('\n(b) 176 landmark over-representation in the reference DE pool')
lm = np.load('experiments/176-cmap-donors/work/landmark_mask.npy')
print(f'    landmark_mask shape {lm.shape} dtype {lm.dtype} sum {int(lm.sum())}')
if lm.shape[0]==len(genes):
    lmg = lm.astype(bool)
    base = lmg.mean()
    inpool = pool[:, lmg].sum()/pool.sum()
    print(f'    landmark genes are {base:.3%} of the DE-universe axis')
    print(f'    landmark pairs are {inpool:.3%} of the reference DE pool  -> enrichment {inpool/base:.2f}x')
    print(f'    176 assumed ~13% of pool pairs are landmarks; measured {inpool:.1%}')
    cov = 82/300
    for share in (base, inpool):
        print(f'      coverage arithmetic at share {share:.3f}: 82/300 x {share:.3f} x 0.019 = {cov*share*0.019:+.6f}')
