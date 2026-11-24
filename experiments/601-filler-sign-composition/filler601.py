"""601: do the reference signs of NON-DE high-CPM genes ("filler" calls) follow a compositional shift?

fid counts our calls against the reference lfc sign whether or not the reference calls the gene. Our ~2,200-gene CPM
block is mostly filler for small-n_real targets. Hypothesis (council 600): per-cell CPM sums to 1e6, so when a target's
DE genes lose (gain) net mass, every other gene's CPM rises (falls) a little; for high-CPM genes that shift can be
larger than the sampling noise of a 400-cell mean, so the filler's reference sign is predictable from the DE set's net
mass change. Analysis of h1-test's own truth (structure only; no parameter is fitted here).
"""
import numpy as np, pandas as pd, anndata as ad, pyarrow.parquet as pq, scipy.sparse as sp, json, sys

a = ad.read_h5ad('data/harness/h1-test.h5ad', backed='r')
genes = np.asarray(a.var_names)
ctrl_idx = np.where(a.obs['target_gene'].values == 'non-targeting')[0]
X = a.X[ctrl_idx]                      # sparse rows of control cells
X = sp.csr_matrix(X, dtype=np.float64)
lib = np.asarray(X.sum(1)).ravel()
cpm = sp.diags(1e6 / lib) @ X
m_ctrl = pd.Series(np.asarray(cpm.mean(0)).ravel(), index=genes)
mc = pd.Series(np.asarray(X.mean(0)).ravel(), index=genes)   # mean count per cell
print('controls', len(ctrl_idx), 'median lib', np.median(lib))

de = pq.read_table('experiments/323-current-family-validation/work/reference_de.parquet').to_pandas()
rows = []
for t, d in de.groupby('target'):
    d = d[d.feature != t].copy()
    d = d[np.isfinite(d.log2_fold_change) & (d.log2_fold_change != 0)]
    d['m'] = m_ctrl.reindex(d.feature).values
    d['cnt'] = mc.reindex(d.feature).values
    sig = d.p_adj < 0.05
    S = d[sig]; N = d[~sig]
    dM_sig = float((S.m * (2.0 ** S.log2_fold_change - 1)).sum())
    dM_all = float((d.m * (2.0 ** d.log2_fold_change - 1)).sum())
    top = d.sort_values('m', ascending=False).head(2200)
    F = top[top.p_adj >= 0.05]
    upF = float((F.log2_fold_change > 0).mean()) if len(F) else np.nan
    rule = -np.sign(dM_sig) if dM_sig != 0 else 0
    acc_rule = float(((np.sign(F.log2_fold_change) == rule)).mean()) if len(F) and rule != 0 else np.nan
    # high-count filler only (mean count/cell >= 5): where a small shift beats a 400-cell mean's noise
    Fh = F[F.cnt >= 5]
    upFh = float((Fh.log2_fold_change > 0).mean()) if len(Fh) else np.nan
    acc_rule_h = float(((np.sign(Fh.log2_fold_change) == rule)).mean()) if len(Fh) and rule != 0 else np.nan
    rows.append(dict(target=t, n_real=int(sig.sum()), dM_sig=dM_sig, dM_all=dM_all, nF=len(F), upF=upF,
                     acc_rule=acc_rule, nFh=len(Fh), upFh=upFh, acc_rule_h=acc_rule_h,
                     medlfcF=float(F.log2_fold_change.median()) if len(F) else np.nan))
r = pd.DataFrame(rows).sort_values('n_real')
pd.set_option('display.width', 220)
print(r.to_string(index=False, float_format=lambda x: f'{x:.4f}'))
print()
print('targets', len(r))
print('filler up-fraction: mean %.4f  median %.4f' % (r.upF.mean(), r.upF.median()))
print('all-down accuracy on filler (mean over targets): %.4f' % (1 - r.upF).mean())
print('all-up   accuracy on filler (mean over targets): %.4f' % r.upF.mean())
print('rule -sign(dM_sig) accuracy on filler: %.4f   high-count filler: %.4f' % (r.acc_rule.mean(), r.acc_rule_h.mean()))
print('spearman(upF, -dM_sig) = %.3f' % r.upF.corr(-r.dM_sig, method='spearman'))
print('spearman(medlfcF, -dM_sig) = %.3f' % r.medlfcF.corr(-r.dM_sig, method='spearman'))
print('dM_all (should be ~0 by the CPM constraint): median %.1f CPM, |max| %.1f' % (r.dM_all.median(), r.dM_all.abs().max()))
r.to_csv('experiments/601-filler-sign-composition/filler601.csv', index=False)
