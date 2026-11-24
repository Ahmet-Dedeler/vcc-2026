"""601c: can a line's OWN control cells orient the generic (target-agnostic) per-gene response sign of the filler?

160 found the line's mean perturbation response sits mostly in the top PCs of its own control cells (captured 0.71 /
0.89 at k 50; it was useless for the bulk because the mean response is pds's hub). The DE channel's filler signs are
exactly that hub (d4cd5f: fid = precision, pinned by target-independent block signs). Here: H1 truth evaluates;
predictors use only H1 controls + donor lines' mean responses (no h1_2025 donor: that would leak).
"""
import numpy as np, pandas as pd, anndata as ad, pyarrow.parquet as pq, scipy.sparse as sp, json
from sklearn.utils.extmath import randomized_svd
rng = np.random.default_rng(601)
a = ad.read_h5ad('data/harness/h1-test.h5ad', backed='r')
genes = np.asarray(a.var_names)
ci = np.where(a.obs['target_gene'].values == 'non-targeting')[0]
X = sp.csr_matrix(a.X[ci], dtype=np.float64); lib = np.asarray(X.sum(1)).ravel()
cpm = sp.diags(1e6 / lib) @ X
m_ctrl = pd.Series(np.asarray(cpm.mean(0)).ravel(), index=genes)
de = pq.read_table('experiments/323-current-family-validation/work/reference_de.parquet').to_pandas()
de = de[(de.feature != de.target) & np.isfinite(de.log2_fold_change) & (de.log2_fold_change != 0)]
block = m_ctrl.reindex(de.feature.unique()).sort_values(ascending=False).head(2200).index
L = de.pivot(index='target', columns='feature', values='log2_fold_change')[block]
P = de.pivot(index='target', columns='feature', values='p_adj')[block]
S = np.sign(L.values); filler = (P.values >= 0.05) & np.isfinite(L.values)

# control PCs in the scorer's per-cell space log1p(TS * p), cell-centred, on the 6,000 highest-CPM genes
hv = m_ctrl.sort_values(ascending=False).head(6000).index
gi = pd.Index(genes).get_indexer(hv)
Y = np.log1p(5e4 * (sp.diags(1.0 / lib) @ X)[:, gi].toarray()).astype(np.float32)
Y -= Y.mean(0)
U, s, Vt = randomized_svd(Y, n_components=50, random_state=601)
V = pd.DataFrame(Vt.T, index=hv)            # gene loadings
print('control PCs: top-5 var share', np.round(s[:5] ** 2 / (Y ** 2).sum(), 4))

def acc(pred):                               # pred: Series gene -> +1/-1 (NaN = undefined)
    p = pred.reindex(block).values
    out, dn = [], []
    for i in range(S.shape[0]):
        f = filler[i] & np.isfinite(p)
        if f.sum() < 20: continue
        out.append(float((np.sign(p[f]) == S[i, f]).mean())); dn.append(float((S[i, f] < 0).mean()))
    return np.mean(out), np.mean(dn), int(np.isfinite(p).sum())

res = {}
# LOTO ceiling and all-down on the full block (reference)
up = (S > 0).astype(float); n = np.isfinite(L.values).astype(float)
tot_up = np.where(np.isfinite(L.values), up, 0).sum(0); tot_n = n.sum(0)
lo = []
for i in range(S.shape[0]):
    f = filler[i]
    if f.sum() == 0: continue
    loo = (tot_up - up[i] * n[i]) / np.maximum(tot_n - n[i], 1)
    lo.append(float((np.where(loo > 0.5, 1, -1)[f] == S[i, f]).mean()))
print('ceiling: LOTO same-line majority %.4f' % np.mean(lo))

for dn in ['k562_gwps', 'nadig_jurkat', 'rpe1', 'orion_hct116', 'orion_hek293t']:
    try:
        D = np.load(f'data/donors/{dn}/delta.npy', mmap_mode='r')
        dg = [l.strip() for l in open(f'data/donors/{dn}/genes.txt')]
    except Exception as e:
        print(dn, 'skip', e); continue
    m = pd.Series(np.asarray(D.mean(0), dtype=np.float64), index=dg)   # donor mean response (log-pseudobulk delta)
    m = m[~m.index.duplicated()]
    a_raw = acc(np.sign(m))
    # recipient-adapted: project the donor mean response onto H1's control PC span (top k), genes in hv
    mv = m.reindex(hv).fillna(0.0).values
    line = {}
    for k in (1, 3, 5, 10, 20, 50):
        Vk = V.values[:, :k]
        proj = pd.Series(Vk @ (Vk.T @ mv), index=hv)
        line[k] = acc(np.sign(proj))
    # PC1 oriented by the donor mean response
    pc1 = V[0]; orient = np.sign(float(pc1.values @ mv)) or 1.0
    a_pc1 = acc(np.sign(pc1 * orient))
    print(f'{dn:14s} raw donor-mean sign {a_raw[0]:.4f} (all-down same genes {a_raw[1]:.4f}, n {a_raw[2]})  '
          + '  '.join(f'proj k{k} {v[0]:.4f}' for k, v in line.items()) + f'  | PC1 oriented {a_pc1[0]:.4f}')
    res[dn] = dict(raw=a_raw, proj={k: v for k, v in line.items()}, pc1=a_pc1)
json.dump(res, open('experiments/601-filler-sign-composition/ctrlprog601.json', 'w'), indent=1, default=float)
