"""Extract the 17,495 needed pseudobulk rows from the downloaded GWCD4i.pseudobulk_merged.h5ad into a
compact dense float32 bank, then the big file can be deleted.  Streaming, block-wise; peak RSS < 1.5 GiB."""
import numpy as np, h5py, time, json, resource, hashlib
ROOT = "~/Code/virtual-cell-challenge-2026"
PARTS = [f"parts/p{i}" for i in range(8)]
z = np.load("obs280.npz", allow_pickle=True)
indptr = z["indptr"]; pert = z["pert"].astype(str); gtype = z["gtype"].astype(str)
panel = sorted(set(l.split(",")[0] for l in open(f"{ROOT}/data/challenge/pert_counts.csv").read().strip().split("\n")[1:]))
sel = np.where((np.isin(pert, panel) & (gtype == "targeting")) | (gtype != "targeting"))[0]
n_g = len(z["genes"])
print(f"selecting {len(sel):,} rows x {n_g} genes -> {len(sel)*n_g*4/1e9:.2f} GB float32")
from parts_file import PartsFile
pf = PartsFile(PARTS)
assert pf.size == 44566657140, pf.size
f = h5py.File(pf, "r")
D, I = f["X/data"], f["X/indices"]
out = np.lib.format.open_memmap("bank280_X.npy", mode="w+", dtype=np.float32, shape=(len(sel), n_g))
t0 = time.time()
runs = np.split(sel, np.where(np.diff(sel) != 1)[0] + 1)
w = 0
maxval = 0.0
for k, run in enumerate(runs):
    a, b = int(indptr[run[0]]), int(indptr[run[-1] + 1])
    dat = np.asarray(D[a:b]); idx = np.asarray(I[a:b])
    off = indptr[run[0]]
    for r in run:
        s, e = int(indptr[r] - off), int(indptr[r + 1] - off)
        out[w, idx[s:e]] = dat[s:e]
        w += 1
    maxval = max(maxval, float(dat.max()) if len(dat) else 0.0)
    if k % 300 == 0:
        print(f"  run {k}/{len(runs)} rows {w}/{len(sel)} {time.time()-t0:.0f}s "
              f"rss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**3:.2f} GiB", flush=True)
out.flush(); del out
np.save("bank280_rows.npy", sel)
print(f"done {time.time()-t0:.0f}s, max count value {maxval:.0f} (float32 exact to 16,777,216); "
      f"rss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**3:.2f} GiB")
