"""813: quantiles of named numeric obs columns of a remote h5ad via HTTP range reads (reuses h5peek813.RangeFile).
Usage: uv run python obsq813.py <url> <out.json> <cap_mb> col1 [col2 ...]"""
import json, sys
sys.argv, rest = sys.argv[:1] + sys.argv[1:4], sys.argv[4:]
import importlib.util, numpy as np, h5py
spec = importlib.util.spec_from_file_location("p", __file__.replace("obsq813.py", "h5peek813.py"))
src = open(spec.origin).read().split("rf = RangeFile(URL)")[0]
ns = {"__name__": "x"}; exec(src, ns)
rf = ns["RangeFile"](sys.argv[1]); out = {}
with h5py.File(rf, "r") as f:
    for c in rest:
        a = f["obs"][c][:]
        q = np.quantile(a, [0, .1, .25, .5, .75, .9, 1]).tolist()
        out[c] = {"n": int(a.size), "mean": float(a.mean()), "q0_10_25_50_75_90_100": q}
        if a.dtype.kind in "iu" and a.max() < 50:
            out[c]["hist"] = np.bincount(a.astype(int)).tolist()
out["bytes_fetched"] = rf.fetched
json.dump(out, open(sys.argv[2], "w"), indent=1); print(json.dumps(out, indent=1))
