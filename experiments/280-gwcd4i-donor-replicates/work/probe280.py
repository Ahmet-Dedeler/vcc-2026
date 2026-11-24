"""Inspect the GWCD4i S3 objects' structure over HTTP range requests -- no download."""
import sys, time, json
import numpy as np, h5py
from remote_h5 import RemoteFile

B = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"

def walk(g, prefix="", depth=0, maxdepth=3, out=None):
    for k, v in g.items():
        p = f"{prefix}/{k}"
        if isinstance(v, h5py.Dataset):
            out.append((p, "D", str(v.shape), str(v.dtype), str(v.chunks), str(v.compression)))
        else:
            out.append((p, "G", dict(v.attrs).get("encoding-type", b"").__str__(), "", "", ""))
            if depth < maxdepth:
                walk(v, p, depth + 1, maxdepth, out)

name = sys.argv[1]
t0 = time.time()
rf = RemoteFile(f"{B}/{name}")
print(f"{name}: {rf.size:,} bytes")
f = h5py.File(rf, "r")
out = []
walk(f, maxdepth=int(sys.argv[2]) if len(sys.argv) > 2 else 3, out=out)
for p, kind, a, b, c, d in out:
    print(f"  {kind} {p:<70} {a:<24} {b:<10} chunks={c:<18} {d}")
print(f"[{time.time()-t0:.1f}s, {rf.n_requests} range requests, {rf.n_bytes/1e6:.1f} MB fetched]")
