"""601d: is H1's generic filler sign present for NULL targets too (a protocol signature) or only for strong ones?"""
import numpy as np, pandas as pd, anndata as ad, pyarrow.parquet as pq, scipy.sparse as sp
a = ad.read_h5ad('data/harness/h1-test.h5ad', backed='r'); genes = np.asarray(a.var_names)
ci = np.where(a.obs['target_gene'].values == 'non-targeting')[0]
X = sp.csr_matrix(a.X[ci], dtype=np.float64); lib = np.asarray(X.sum(1)).ravel()
m_ctrl = pd.Series(np.asarray((sp.diags(1e6 / lib) @ X).mean(0)).ravel(), index=genes)
de = pq.read_table('experiments/323-current-family-validation/work/reference_de.parquet').to_pandas()
de = de[(de.feature != de.target) & np.isfinite(de.log2_fold_change) & (de.log2_fold_change != 0)]
nreal = de.assign(s=de.p_adj < 0.05).groupby('target').s.sum()
block = m_ctrl.reindex(de.feature.unique()).sort_values(ascending=False).head(2200).index
L = de.pivot(index='target', columns='feature', values='log2_fold_change')[block]
P = de.pivot(index='target', columns='feature', values='p_adj')[block]
S = np.sign(L.values); fil = (P.values >= 0.05) & np.isfinite(L.values)
tg = L.index.values; nr = nreal.reindex(tg).values
strata = {'null (n_real<=5)': nr <= 5, 'weak (6-100)': (nr > 5) & (nr <= 100), 'strong (>100)': nr > 100}
up = (S > 0).astype(float); n = np.isfinite(L.values).astype(float)
for src_name, src in strata.items():
    # majority computed ONLY from the source stratum's targets (leave-one-out within it)
    tu = (np.where(np.isfinite(L.values), up, 0) * src[:, None]).sum(0); tn = (n * src[:, None]).sum(0)
    for ev_name, ev in strata.items():
        acc = []; dn = []
        for i in np.where(ev)[0]:
            f = fil[i]
            u, m_ = tu - (up[i] * n[i] if src[i] else 0), tn - (n[i] if src[i] else 0)
            pred = np.where(u / np.maximum(m_, 1) > 0.5, 1, -1)
            acc.append(float((pred[f] == S[i, f]).mean())); dn.append(float((S[i, f] < 0).mean()))
        print(f'majority from {src_name:16s} ({int(src.sum()):2d} tg) -> evaluated on {ev_name:16s}: acc {np.mean(acc):.4f}  (all-down {np.mean(dn):.4f})')
