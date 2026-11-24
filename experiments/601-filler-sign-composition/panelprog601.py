"""601e: does a PANEL-RESTRICTED donor programme (mean of the donor rows of THIS panel's targets) orient the filler
better than the donor's genome-wide mean? H1 truth evaluates; donors are other lines (no h1_2025)."""
import numpy as np, pandas as pd, anndata as ad, pyarrow.parquet as pq, scipy.sparse as sp
a = ad.read_h5ad('data/harness/h1-test.h5ad', backed='r'); genes = np.asarray(a.var_names)
ci = np.where(a.obs['target_gene'].values == 'non-targeting')[0]
X = sp.csr_matrix(a.X[ci], dtype=np.float64); lib = np.asarray(X.sum(1)).ravel()
m_ctrl = pd.Series(np.asarray((sp.diags(1e6 / lib) @ X).mean(0)).ravel(), index=genes)
de = pq.read_table('experiments/323-current-family-validation/work/reference_de.parquet').to_pandas()
de = de[(de.feature != de.target) & np.isfinite(de.log2_fold_change) & (de.log2_fold_change != 0)]
panel = sorted(de.target.unique())
block = m_ctrl.reindex(de.feature.unique()).sort_values(ascending=False).head(2200).index
L = de.pivot(index='target', columns='feature', values='log2_fold_change')[block]
P = de.pivot(index='target', columns='feature', values='p_adj')[block]
S = np.sign(L.values); fil = (P.values >= 0.05) & np.isfinite(L.values)
def acc(pred, rows=None):
    p = pred.reindex(block).values; out = []; dn = []
    for i in range(S.shape[0]):
        if rows is not None and not rows[i]: continue
        f = fil[i] & np.isfinite(p) & (p != 0)
        if f.sum() < 20: continue
        out.append(float((np.sign(p[f]) == S[i, f]).mean())); dn.append(float((S[i, f] < 0).mean()))
    return np.mean(out), np.mean(dn)
for dn in ['k562_gwps', 'nadig_jurkat', 'rpe1', 'orion_hct116', 'orion_hek293t', 'k562_essential']:
    try:
        D = np.load(f'data/donors/{dn}/delta.npy', mmap_mode='r')
        dg = [l.strip() for l in open(f'data/donors/{dn}/genes.txt')]
        pt = [l.strip() for l in open(f'data/donors/{dn}/perts.txt')]
    except Exception as e:
        print(dn, 'skip', e); continue
    pt = pd.Index(pt); idx = [i for i, t in enumerate(pt) if t in set(panel)]
    glob = pd.Series(np.asarray(D.mean(0), dtype=float), index=dg); glob = glob[~glob.index.duplicated()]
    if len(idx) < 10:
        print(f'{dn:14s} covers {len(idx)} panel targets -- skip panel mean'); continue
    pm = pd.Series(np.asarray(D[idx].mean(0), dtype=float), index=dg); pm = pm[~pm.index.duplicated()]
    # LOTO-free: the panel mean uses donor rows of the panel's targets; evaluating target i includes donor row of i,
    # which is legitimate (on the day we have donor rows for every covered panel target)
    a_g = acc(np.sign(glob)); a_p = acc(np.sign(pm))
    print(f'{dn:14s} covers {len(idx):3d}/{len(panel)} panel targets | genome-wide mean sign {a_g[0]:.4f} | panel-mean sign {a_p[0]:.4f} | all-down {a_p[1]:.4f}')
