import numpy as np, anndata as ad, scanpy as sc, scipy.sparse as sp
from scipy import stats
import sys
sys.path.insert(0, ".")
from de_kernel import ctrl_tiesums, wilcoxon_z

rng = np.random.default_rng(0)
m, n, G = 400, 60, 500
Xc = sp.random(m, G, density=0.25, format="csr", random_state=1,
               data_rvs=lambda k: rng.integers(1, 30, k).astype(np.float64))
Xp = sp.random(n, G, density=0.30, format="csr", random_state=2,
               data_rvs=lambda k: rng.integers(1, 30, k).astype(np.float64))
X = sp.vstack([Xp, Xc]).tocsr()
labels = np.array(["pert"] * n + ["ctrl"] * m)
A = ad.AnnData(X=X.copy(), obs={"g": labels})
A.X = A.X.astype(float)
sc.pp.normalize_total(A, target_sum=1e6)
lin = A.X.copy()
sc.pp.log1p(A)
sc.tl.rank_genes_groups(A, groupby="g", reference="ctrl", method="wilcoxon",
                        tie_correct=True, n_genes=A.n_vars)
res = A.uns["rank_genes_groups"]
names = np.asarray(res["names"]["pert"], dtype=str)
pv_sc = np.asarray(res["pvals"]["pert"], dtype=float)
order = np.argsort(names.astype(int) if names[0].isdigit() else names)
# map by gene name
gene_ids = np.asarray(A.var_names, dtype=str)
pos = {g: i for i, g in enumerate(gene_ids)}
pv_ref = np.empty(G); pv_ref[[pos[x] for x in names]] = pv_sc

D = np.asarray(A.X.todense(), dtype=np.float64)   # log1p CPM
ctrl = np.ascontiguousarray(D[n:].T)              # (G, m)
pert = np.ascontiguousarray(D[:n].T)              # (G, n)
cs = np.sort(ctrl, axis=1)
ct = ctrl_tiesums(cs)
z = wilcoxon_z(cs, ct, pert)
pv_mine = 2 * stats.distributions.norm.sf(np.abs(z))

print("max abs diff p:", np.abs(pv_mine - pv_ref).max())
print("max rel diff p:", (np.abs(pv_mine - pv_ref) / np.maximum(pv_ref, 1e-300)).max())
print("exact equal:", np.array_equal(pv_mine, pv_ref))
# also z-scores
zs = np.asarray(res["scores"]["pert"], dtype=float)
z_ref = np.empty(G); z_ref[[pos[x] for x in names]] = zs
print("max abs diff z:", np.abs(z - z_ref).max())
