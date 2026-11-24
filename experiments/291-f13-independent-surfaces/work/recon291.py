"""291 recon — COUNTS AND COLUMN AVAILABILITY ONLY (no feature-outcome number).
Run before the preregistration is written, to fix n and the excluded-source list."""
import json
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path("~/Code/virtual-cell-challenge-2026")
z = np.load(ROOT/'experiments/290-remainder-response-features/work/persrc290.npz', allow_pickle=True)
genes = [str(g) for g in z['genes']]; keys=[str(x) for x in z['keys']]; lines=[str(x) for x in z['lines']]
held = set(genes)
k562_srcs = [k for k,l in zip(keys,lines) if l=='K562']
print('K562 sources', k562_srcs, '-> K562-free banks', len(keys)-len(k562_srcs))

inv = json.load(open(ROOT/'experiments/268-broad16-context-residuals/work/inventory.json'))
tg = inv['targets']; controls=set(inv['controls']); complete=inv['complete']
print('broad targets', len(tg), 'of which in held map', len(set(tg)&held))
print('  analysis (non-OR) targets', len([t for t in tg if t not in controls]),
      'in held map', len(set(t for t in tg if t not in controls)&held))
print('  complete', len(complete), 'in held map', len(set(complete)&held))
print('  complete non-OR in held', len(set(c for c in complete if c not in controls)&held))
pt = pd.read_csv(ROOT/'experiments/268-broad16-context-residuals/work/per_target_decomp.csv')
print('per_target_decomp rows', len(pt), 'in held', len(set(pt['target'])&held))

k5 = pd.read_csv(ROOT/'data/external/cd4_cell2026/suppl_tables_2026-08-24/K562_comparison.suppl_table.csv', low_memory=False)
print('\nK562 table rows', len(k5), 'conditions', sorted(k5['condition'].unique()))
print('unique genes', k5['target_contrast_gene_name'].nunique(), 'in held map',
      len(set(k5['target_contrast_gene_name'])&held))
for c in ['Rest','Stim8hr','Stim48hr']:
    s = k5[k5['condition']==c]
    n_r = s['logfc_pearson_r'].notna().sum()
    n_h = len(set(s.dropna(subset=['logfc_pearson_r'])['target_contrast_gene_name'])&held)
    print(f'  {c:9s} rows {len(s)} with r {n_r} genes-in-held-with-r {n_h}',
          'donor_corr non-null', s['donor_correlation_mean'].notna().sum(),
          'n_degs_K562 non-null', s['n_degs_MASH_K562'].notna().sum())
# 284 component / family labels coverage
p284 = pd.read_csv(ROOT/'experiments/284-conservation-family-split/work/pred284.csv', low_memory=False)
print('\npred284 rows', len(p284), 'components', p284['component'].nunique())
b = set(tg)&held; g = set(k5['target_contrast_gene_name'])&held
print('components spanned by Broad-16 unit', p284[p284['gene'].isin(b)]['component'].nunique())
print('components spanned by GWCD4i unit', p284[p284['gene'].isin(g)]['component'].nunique())
print('overlap Broad16 ∩ GWCD4i units', len(b&g))
