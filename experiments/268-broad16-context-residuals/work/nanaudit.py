import json, numpy as np
d = np.load("pb_cache.npz", allow_pickle=True); X, genes, lines, targets = d["X"], d["genes"], d["lines"], d["targets"]
res = json.load(open("reserve.json")); work = res["working"]
ul = sorted(set(lines))
nanl = {}
for l in ul:
    cols = np.where(lines == l)[0]; sub = X[:, cols]
    rows_nan = np.isnan(sub).any(1); rows_allnan = np.isnan(sub).all(1)
    nanl[l] = (int(rows_nan.sum()), int(rows_allnan.sum()), int((~rows_nan).sum()))
print("per line: genes with any NaN / all NaN / complete:", nanl)
wmask = np.ones(len(genes), bool); amask = np.ones(len(genes), bool)
for l in ul:
    cols = np.where(lines == l)[0]; m = ~np.isnan(X[:, cols]).any(1)
    amask &= m
    if l in work: wmask &= m
print("genes complete across 16 lines:", int(amask.sum()), " across 12 working:", int(wmask.sum()))
# within-line: is NaN gene-wise (a gene missing for all columns of a line) or scattered?
l = ul[0]; cols = np.where(lines == l)[0]; sub = X[:, cols]; nn = np.isnan(sub).sum(1)
print(l, "NaN count per gene distribution (0, partial, all):", int((nn == 0).sum()), int(((nn > 0) & (nn < len(cols))).sum()), int((nn == len(cols)).sum()))
depth = np.nansum(X[wmask], 0)
print("depth on 12-line-complete genes: min %.3g median %.3g max %.3g" % (depth.min(), np.median(depth), depth.max()))
for l in ul:
    cols = np.where(lines == l)[0]
    print(f"{l:9s} depth median {np.median(depth[cols]):.3g}  min {depth[cols].min():.3g}  OR-control depth median {np.median(depth[[c for c in cols if targets[c].startswith('OR')]]):.3g}")
