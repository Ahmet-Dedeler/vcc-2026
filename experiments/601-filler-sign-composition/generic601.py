"""601b: how much of the filler's reference sign is a per-gene GENERIC direction shared across targets (H1 truth)?
Leave-one-target-out majority sign per gene (the ceiling of any per-gene, target-agnostic sign rule in this line),
compared with all-down and the compositional rule; restricted to the top-2,200-CPM genes outside each target's DE set."""
import numpy as np, pandas as pd, anndata as ad, pyarrow.parquet as pq, scipy.sparse as sp

a = ad.read_h5ad('data/harness/h1-test.h5ad', backed='r')
genes = np.asarray(a.var_names)
ci = np.where(a.obs['target_gene'].values == 'non-targeting')[0]
X = sp.csr_matrix(a.X[ci], dtype=np.float64); lib = np.asarray(X.sum(1)).ravel()
m_ctrl = pd.Series(np.asarray((sp.diags(1e6 / lib) @ X).mean(0)).ravel(), index=genes)
de = pq.read_table('experiments/323-current-family-validation/work/reference_de.parquet').to_pandas()
de = de[(de.feature != de.target) & np.isfinite(de.log2_fold_change) & (de.log2_fold_change != 0)]
top = m_ctrl.reindex(de.feature.unique()).sort_values(ascending=False).head(2200).index
L = de.pivot(index='target', columns='feature', values='log2_fold_change')[top]
P = de.pivot(index='target', columns='feature', values='p_adj')[top]
S = np.sign(L.values)                       # target x gene reference signs
filler = (P.values >= 0.05) & np.isfinite(L.values)
up = (S > 0).astype(float); n = np.isfinite(L.values).astype(float)
tot_up = np.nansum(np.where(np.isfinite(L.values), up, 0), 0); tot_n = n.sum(0)
acc_loto, acc_down, acc_dm = [], [], []
for i in range(S.shape[0]):
    f = filler[i]
    if f.sum() == 0: continue
    loo_up = (tot_up - up[i] * n[i]) / np.maximum(tot_n - n[i], 1)
    pred = np.where(loo_up > 0.5, 1, -1)
    acc_loto.append(float((pred[f] == S[i, f]).mean()))
    acc_down.append(float((S[i, f] < 0).mean()))
print('targets', len(acc_loto))
print('filler accuracy, mean over targets:  all-down %.4f   LOTO per-gene majority %.4f' % (np.mean(acc_down), np.mean(acc_loto)))
cons = np.abs(tot_up / np.maximum(tot_n, 1) - 0.5)
for q in (0.1, 0.2):
    sel = cons >= q
    acc = []
    for i in range(S.shape[0]):
        f = filler[i] & sel
        if f.sum() < 5: continue
        loo_up = (tot_up - up[i] * n[i]) / np.maximum(tot_n - n[i], 1)
        pred = np.where(loo_up > 0.5, 1, -1)
        acc.append(float((pred[f] == S[i, f]).mean()))
    print(f'  genes with |up-frac - 0.5| >= {q}: {int(sel.sum())} genes, LOTO majority accuracy on their filler {np.mean(acc):.4f}')
