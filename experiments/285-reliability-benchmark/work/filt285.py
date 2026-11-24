import sys, time, signal
from pathlib import Path
import h5py, numpy as np
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/280-gwcd4i-donor-replicates/work")
from remote_h5 import RemoteFile
B = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
rf = RemoteFile(f"{B}/D2_Stim48hr.assigned_guide.h5ad", block=1<<22, max_blocks=64)
f = h5py.File(rf, "r")
for nm in ("X/data","X/indices","X/indptr"):
    d = f[nm]; p = d.id.get_create_plist(); nf = p.get_nfilters()
    print(nm, d.dtype, d.shape, "chunks", d.chunks, "nfilters", nf, "layout", p.get_layout(),
          "alloc", p.get_alloc_time(), flush=True)
    for i in range(nf):
        print("   filter", p.get_filter(i), flush=True)
t0=time.time()
print("data[0:8]", f["X/data"][0:8], f"{time.time()-t0:.1f}s", flush=True)
t0=time.time()
print("indices[0:8]", f["X/indices"][0:8], f"{time.time()-t0:.1f}s", flush=True)
print("requests", rf.n_requests, "MB", rf.n_bytes/1e6)
