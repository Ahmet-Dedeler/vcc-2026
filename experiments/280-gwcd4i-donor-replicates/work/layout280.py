"""Is X/data and X/indices laid out linearly in chunk order?  If so we can compute exact byte ranges."""
import numpy as np, h5py, time
from remote_h5 import RemoteFile
B = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
rf = RemoteFile(f"{B}/GWCD4i.pseudobulk_merged.h5ad", block=1 << 20, max_blocks=512)
f = h5py.File(rf, "r")
t0 = time.time()
for nm in ("X/data", "X/indices"):
    d = f[nm]
    dsid = d.id
    n = dsid.get_num_chunks()
    print(f"{nm}: shape {d.shape} chunks {d.chunks} n_chunks {n} itemsize {d.dtype.itemsize} compression {d.compression}")
    infos = [dsid.get_chunk_info(i) for i in (0, 1, 2, 100, n // 2, n - 2, n - 1)]
    base = infos[0].byte_offset
    cb = d.chunks[0] * d.dtype.itemsize
    for i, inf in zip((0, 1, 2, 100, n // 2, n - 2, n - 1), infos):
        pred = base + i * cb
        print(f"   chunk {i:>6}: offset {inf.byte_offset:>15,}  size {inf.size:>9,}  predicted {pred:>15,}  delta {inf.byte_offset-pred:>+12,}")
print(f"[{time.time()-t0:.1f}s, {rf.n_requests} req, {rf.n_bytes/1e6:.1f} MB]")
