"""Find the cell->target assignment rule that reproduces the on-disk pseudobulk_sum_by_ko column (known-answer for the parse)."""
import sys, time, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
line = sys.argv[1]; gefile = sys.argv[2]; gfile = sys.argv[3]
t0 = time.time()
d = np.load("pb_cache.npz", allow_pickle=True); X, genes, lines, targets = d["X"], d["genes"], d["lines"], d["targets"]
cols = np.where(lines == line)[0]; pb = pd.DataFrame(X[:, cols], index=genes, columns=targets[cols])
g = pd.read_csv(gfile, index_col=0); g.columns = [c.split("-")[0] for c in g.columns]
a = sc.read_10x_h5(gefile, gex_only=False); a.obs_names = [b.split("-")[0] for b in a.obs_names]
print(line, "GE cells", a.n_obs, "guide barcodes", g.shape[1], f"{time.time()-t0:.0f}s")
common = [b for b in a.obs_names if b in set(g.columns)]
a = a[common].copy(); G = g[common].to_numpy()  # 200 x cells
gnames = np.array(g.index); gtarget = np.array([n.split("_")[0] for n in gnames])
order = np.argsort(-G, axis=0); top = order[0]; second = order[1]
topc = G[top, np.arange(G.shape[1])]; secc = G[second, np.arange(G.shape[1])]; tot = G.sum(0)
# gene mapping: unique-name genes present in the pseudobulk axis
vn = a.var_names.to_numpy(); uniq = pd.Series(vn).duplicated(keep=False).to_numpy() == False
pbset = set(pb.index[~pb.iloc[:, 0].isna()]); keep = np.array([u and (v in pbset) for u, v in zip(uniq, vn)])
Xc = sp.csr_matrix(a.X)[:, keep]; gk = vn[keep]; pbk = pb.loc[gk]
print("genes compared", len(gk), "of h5", len(vn), "pb non-NaN", len(pbset))
def evaluate(mask, label):
    tg = gtarget[top]; res = []; ncell = []
    for t in pbk.columns:
        m = mask & (tg == t)
        s = np.asarray(Xc[m].sum(0)).ravel(); ref = pbk[t].to_numpy()
        res.append((np.array_equal(s, ref), float(np.abs(s - ref).sum() / max(ref.sum(), 1)), float(np.corrcoef(s, ref)[0, 1]) if s.sum() > 0 else 0)); ncell.append(int(m.sum()))
    ex = sum(r[0] for r in res); dev = np.median([r[1] for r in res]); cor = np.median([r[2] for r in res])
    print(f"{label:38s} cells/KO mean {np.mean(ncell):7.2f}  exact {ex:3d}/{len(res)}  median rel.dev {dev:.4f}  median corr {cor:.5f}", flush=True)
    return ex
evaluate((tot > 0) & (secc == 0), "single nonzero guide")
for k in (1, 2, 3, 5, 10, 20):
    evaluate((topc >= k) & (secc == 0), f"top>={k} & second==0")
for k in (1, 3, 5, 10):
    for f in (0.5, 0.8, 0.9):
        evaluate((topc >= k) & (topc / np.maximum(tot, 1) >= f), f"top>={k} & frac>={f}")
cq = pd.read_csv("~/Code/virtual-cell-challenge-2026/data/external/depmap26q1/extracted/processed/cell_quality_covariates.csv").set_index("cell_line")
print("reference n_cells_per_ko", cq.loc[line, "n_cells_per_ko"], "avg_umis", cq.loc[line, "avg_umis_per_ko"])
lib = np.asarray(sp.csr_matrix(a.X).sum(1)).ravel()
print("h5 library size quantiles", np.quantile(lib, [0, .1, .5, .9, 1]))
print(f"done {time.time()-t0:.0f}s")
