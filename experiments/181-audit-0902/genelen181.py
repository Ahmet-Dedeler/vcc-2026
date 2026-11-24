"""181 — test the 'long-gene / nuclear-fraction technical axis' identification (175 addendum,
now in AGENTS.md's harness warning) as a CORRELATION rather than a top-30 gene-name eyeball.

The confound the eyeball cannot see: long intron-rich genes are also LOW-expression genes, and
short housekeeping genes are HIGH-expression.  A pure expression axis would produce exactly the
same top-30 list.  So the test is the PARTIAL correlation of the loading with log(genomic span)
controlling for log(control CPM), and vice versa.
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr

z = np.load('experiments/181-audit-0902/work/pc1.npz', allow_pickle=True)
genes = np.array([str(g) for g in z['genes']])
sp = pd.read_csv('experiments/181-audit-0902/work/gene_spans.csv')
span = dict(zip(sp.gene, sp.span)); bt = dict(zip(sp.gene, sp.biotype))

eff = np.load('experiments/169-sharpening-six-member/work/h1test-effect.npz', allow_pickle=True)
cpm = dict(zip([str(g) for g in eff['genes']], eff['cpm']))

for donor in ('orion_hek293t','orion_hct116','k562_gwps'):
    v = z[donor].astype(float)
    v = v * (-1.0 if v[np.argmax(np.abs(v))] > 0 else 1.0)
    m = np.array([g in span and g in cpm and cpm[g] > 0 for g in genes])
    L = np.log10(np.array([span[g] for g in genes[m]], float))
    E = np.log10(np.array([cpm[g] for g in genes[m]], float))
    V = v[m]
    def part(a, b, c):     # partial corr of a,b controlling c  (Spearman on ranks)
        from scipy.stats import rankdata
        A,B,C = (rankdata(x) for x in (a,b,c))
        A = A - np.polyval(np.polyfit(C,A,1),C); B = B - np.polyval(np.polyfit(C,B,1),C)
        return pearsonr(A,B)[0]
    print(f'--- {donor}  (n={m.sum()} genes with span and cpm) ---')
    print(f'   spearman(loading, log10 span) = {spearmanr(V,L).statistic:+.4f}')
    print(f'   spearman(loading, log10 cpm ) = {spearmanr(V,E).statistic:+.4f}')
    print(f'   spearman(log10 span, log10 cpm) = {spearmanr(L,E).statistic:+.4f}   <- the confound')
    print(f'   PARTIAL spearman(loading, span | cpm) = {part(V,L,E):+.4f}')
    print(f'   PARTIAL spearman(loading, cpm  | span) = {part(V,E,L):+.4f}')
    # decile tables
    q = np.quantile(L,[i/10 for i in range(11)])
    print('   mean loading by span decile: ', ' '.join(f'{V[(L>=q[i])&(L<=q[i+1])].mean():+.4f}' for i in range(10)))
    q2 = np.quantile(E,[i/10 for i in range(11)])
    print('   mean loading by cpm  decile: ', ' '.join(f'{V[(E>=q2[i])&(E<=q2[i+1])].mean():+.4f}' for i in range(10)))
    print()

# how long are the top-30 vs the axis as a whole
v = z['orion_hek293t'].astype(float); v = v*(-1.0 if v[np.argmax(np.abs(v))]>0 else 1.0)
m = np.array([g in span for g in genes]); o = np.argsort(v[m]); gm = genes[m]
S = np.array([span[g] for g in gm], float)
print(f'median genomic span: whole axis {np.median(S)/1e3:.0f} kb; '
      f'top-30 most-negative loadings {np.median(S[o[:30]])/1e3:.0f} kb; '
      f'top-30 most-positive {np.median(S[o[-30:]])/1e3:.0f} kb')
Cm = np.array([cpm.get(g,np.nan) for g in gm], float)
print(f'median control CPM : whole axis {np.nanmedian(Cm):.1f}; '
      f'top-30 negative {np.nanmedian(Cm[o[:30]]):.1f}; top-30 positive {np.nanmedian(Cm[o[-30:]]):.1f}')
