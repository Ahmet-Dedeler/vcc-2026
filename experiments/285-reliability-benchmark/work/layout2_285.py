"""285 -- the CD4 object's X/data and X/indices base offsets, WITHOUT walking the chunk index.

Only `get_chunk_info(0)` (one call per dataset) is used; linearity is then validated the way 280's KA-0
validated its extraction: a row's raw bytes, decoded under the linear assumption, must have exactly
`n_genes_by_counts[r]` entries, strictly increasing indices below n_genes, and a data sum equal to
`obs/total_counts[r]`.  Fifty scattered rows across the whole file.
"""
import json, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import h5py, numpy as np
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/280-gwcd4i-donor-replicates/work")
from remote_h5 import RemoteFile
WORK = Path("~/Code/virtual-cell-challenge-2026/experiments/285-reliability-benchmark/work")
B = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
NAME = sys.argv[1] if len(sys.argv) > 1 else "D2_Stim48hr.assigned_guide.h5ad"
STEM = NAME.split(".")[0]; url = f"{B}/{NAME}"
rf = RemoteFile(url, block=1 << 22, max_blocks=32); f = h5py.File(rf, "r")
ip = np.load(WORK / f"indptr_{STEM}.npy")
z = np.load(WORK / f"probe2_{STEM}.npz", allow_pickle=True)
tot = np.asarray(z["total"], dtype=np.float64); ngc = np.asarray(z["ngc"], dtype=np.int64)
t0 = time.time(); base = {}
for nm in ("X/data", "X/indices"):
    d = f[nm]; inf = d.id.get_chunk_info(0)
    base[nm] = int(inf.byte_offset)
    print(f"{nm}: {d.shape[0]:,} {d.dtype}, chunk0 offset {inf.byte_offset:,} size {inf.size:,}", flush=True)
ng = int(f["X"].attrs["shape"][1])
print(f"[{time.time()-t0:.0f}s, {rf.n_requests} req, {rf.n_bytes/1e6:.0f} MB]", flush=True)
def get(lo, hi):
    req = urllib.request.Request(url, headers={"Range": f"bytes={lo}-{hi-1}"})
    for a in range(4):
        try:
            with urllib.request.urlopen(req, timeout=180) as r: return r.read()
        except Exception:
            if a == 3: raise
rng = np.random.default_rng(20260906)
rows = np.sort(rng.choice(ngc.size, 50, replace=False))
bad = 0; t0 = time.time()
def check(r):
    a, b = int(ip[r]), int(ip[r+1])
    dat = np.frombuffer(get(base["X/data"] + a*4, base["X/data"] + b*4), dtype="<f4")
    idx = np.frombuffer(get(base["X/indices"] + a*8, base["X/indices"] + b*8), dtype="<i8")
    ok = (len(dat) == ngc[r] == len(idx) and (np.diff(idx) > 0).all() and idx.min() >= 0 and idx.max() < ng
          and float(dat.sum()) == tot[r] and (dat > 0).all())
    return r, ok, len(dat), float(dat.sum()), float(tot[r]), int(idx.max())
with ThreadPoolExecutor(max_workers=12) as ex:
    res = list(ex.map(check, rows))
for r, ok, n, s, t, mx in res[:6]:
    print(f"  row {r:>9}: nnz {n} (obs {ngc[r]}), sum {s:.0f} vs total_counts {t:.0f}, max index {mx} -> {'OK' if ok else 'BAD'}")
nbad = sum(1 for _, ok, *_ in res if not ok)
print(f"KA-layout: {len(res)-nbad}/{len(res)} scattered rows decode exactly under the linear assumption "
      f"({time.time()-t0:.0f}s) -> {'PASS' if nbad == 0 else 'FAIL'}")
json.dump(dict(base=base, n_genes=ng, checked=len(res), bad=nbad,
               rows=[[int(r), bool(o), int(n)] for r, o, n, *_ in res]),
          open(WORK / f"layout2_{STEM}.json", "w"), indent=1)
