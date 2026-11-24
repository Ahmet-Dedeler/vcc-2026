import numpy as np, h5py, time
from remote_h5 import RemoteFile
B="https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
t0=time.time()
rf=RemoteFile(f"{B}/GWCD4i.DE_stats.by_donors.h5mu", block=1<<22, max_blocks=64)
f=h5py.File(rf,"r")
print("root keys:", list(f.keys()), f"size {rf.size:,}")
mods=list(f["mod"].keys()); print("modalities:", mods)
m=f[f"mod/{mods[0]}"]
print(" keys:", list(m.keys()))
for k in ("layers","obs","var"):
    if k in m: print(f"  {k}:", list(m[k].keys())[:14])
d=m["layers/log_fc"]
print("  log_fc:", d.shape, d.dtype, "chunks", d.chunks, "compression", d.compression)
print(f"[{time.time()-t0:.0f}s {rf.n_requests} req {rf.n_bytes/1e6:.0f} MB]")
