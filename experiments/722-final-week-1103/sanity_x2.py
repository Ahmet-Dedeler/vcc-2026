"""Independent semantic sanity check of a final-entry prediction h5ad (11-03).

Checks what postbuild does not: that each context's cells look like that context's own controls,
that each target label carries its own knockdown, and basic count hygiene. Streams the CSR by row
blocks (backed mode), so memory stays small.
"""
import sys, json
import numpy as np, pandas as pd, anndata as ad, h5py
from scipy import sparse

pred_path, ctl_dir = sys.argv[1], sys.argv[2]
genes = pd.read_csv(f"{ctl_dir}/gene_names.csv")["gene_name"].tolist()
perts = pd.read_csv(f"{ctl_dir}/pert_counts.csv")["target_gene"].tolist()
gidx = {g: i for i, g in enumerate(genes)}

A = ad.read_h5ad(pred_path, backed="r")
obs = A.obs
assert list(A.var_names) == genes, "var order differs from gene_names.csv"
ctxs = sorted(obs["context"].astype(str).unique())
out = {"n_obs": int(A.n_obs), "contexts": ctxs,
       "per_ctx_target_counts_ok": bool(all(
           (obs[obs["context"].astype(str) == c]["target_gene"].astype(str).value_counts() == 400).all()
           and set(obs[obs["context"].astype(str) == c]["target_gene"].astype(str)) == set(perts) for c in ctxs))}

# Stream the prediction: per (context, target) summed counts, per-cell totals, max entry, integer check.
f = h5py.File(pred_path, "r")
X = f["X"]
indptr = X["indptr"][:]
ctx_arr = obs["context"].astype(str).values
tg_arr = obs["target_gene"].astype(str).values
keys = pd.Series(list(zip(ctx_arr, tg_arr)))
key_codes, key_uni = pd.factorize(keys)
G = len(genes)
sums = np.zeros((len(key_uni), G), dtype=np.float64)
totals = np.zeros(A.n_obs)
nonint = 0; neg = 0; maxv = 0.0
B = 8000
for s in range(0, A.n_obs, B):
    e = min(s + B, A.n_obs)
    a, b = indptr[s], indptr[e]
    data = X["data"][a:b]; ind = X["indices"][a:b]
    M = sparse.csr_matrix((data, ind, indptr[s:e + 1] - a), shape=(e - s, G))
    totals[s:e] = np.asarray(M.sum(1)).ravel()
    nonint += int(np.count_nonzero(data != np.round(data))); neg += int((data < 0).sum())
    maxv = max(maxv, float(data.max()) if data.size else 0.0)
    codes = key_codes[s:e]
    agg = sparse.csr_matrix((np.ones(e - s), (codes, np.arange(e - s))), shape=(len(key_uni), e - s))
    sums += (agg @ M).toarray()
out.update(nonint=nonint, negative=neg, max_entry=maxv, max_cell_total=float(totals.max()))

# Controls: per-context pseudobulk + per-gene control mean CPM.
ctl = {}
for c in ctxs:
    C = ad.read_h5ad(f"{ctl_dir}/context_{c}.h5ad", backed="r")
    tot = np.zeros(G); lib = []
    for s in range(0, C.n_obs, B):
        M = C.X[s:min(s + B, C.n_obs)]
        M = M if sparse.issparse(M) else sparse.csr_matrix(M)
        tot += np.asarray(M.sum(0)).ravel(); lib.append(np.asarray(M.sum(1)).ravel())
    ctl[c] = (tot, np.concatenate(lib))
def logcpm(v): return np.log1p(v / v.sum() * 1e6)

# Context identity: prediction pseudobulk of context X vs controls of each context.
ident = {}
for c in ctxs:
    pb = sums[[i for i, k in enumerate(key_uni) if k[0] == c]].sum(0)
    ident[c] = {d: float(np.corrcoef(logcpm(pb), logcpm(ctl[d][0]))[0, 1]) for d in ctxs}
    ident[c]["lib_median_pred"] = float(np.median(totals[ctx_arr == c]))
    ident[c]["lib_median_ctl"] = float(np.median(ctl[c][1]))
out["identity"] = ident

# Own-target knockdown: CPM of the target gene in its own rows vs the context's controls.
kd = {}
for c in ctxs:
    ctl_cpm = ctl[c][0] / ctl[c][0].sum() * 1e6
    ratios = []
    for i, k in enumerate(key_uni):
        if k[0] != c or k[1] not in gidx: continue
        j = gidx[k[1]]
        r = (sums[i, j] / sums[i].sum() * 1e6) / max(ctl_cpm[j], 1e-9)
        ratios.append(r)
    ratios = np.array(ratios)
    kd[c] = {"n": int(len(ratios)), "median_ratio": float(np.median(ratios)),
             "frac_below_0.5": float((ratios < 0.5).mean()), "frac_above_1": float((ratios > 1).mean())}
out["own_target_knockdown"] = kd

# Label-shuffle check: own-target ratio vs the same gene's ratio in a random other target's rows.
print(json.dumps(out, indent=1))
