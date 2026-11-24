"""Benchmark: how much network traffic and time to pull N selected rows via h5py slicing?"""
import numpy as np, h5py, time
from remote_h5 import RemoteFile
B = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
z = np.load("obs280.npz", allow_pickle=True)
indptr = z["indptr"]; pert = z["pert"]; gtype = z["gtype"]
ROOT = "~/Code/virtual-cell-challenge-2026"
panel = set(l.split(",")[0] for l in open(f"{ROOT}/data/challenge/pert_counts.csv").read().strip().split("\n")[1:])
sel = np.where(np.isin(pert, list(panel)) & (gtype == "targeting") | (gtype != "targeting"))[0]
runs = np.split(sel, np.where(np.diff(sel) != 1)[0] + 1)
for blk in (1 << 19, 1 << 21):
    rf = RemoteFile(f"{B}/GWCD4i.pseudobulk_merged.h5ad", block=blk, max_blocks=64)
    f = h5py.File(rf, "r")
    D, I = f["X/data"], f["X/indices"]
    r0, b0 = rf.n_requests, rf.n_bytes
    t0 = time.time(); nnz = 0
    for run in runs[:40]:
        a, b = int(indptr[run[0]]), int(indptr[run[-1] + 1])
        _ = D[a:b]; _ = I[a:b]; nnz += b - a
    dt = time.time() - t0
    got = (rf.n_bytes - b0)
    print(f"block {blk>>10:>5} KiB: 40 runs, {nnz:,} nnz (payload {nnz*16/1e6:.0f} MB) in {dt:.1f}s, "
          f"{rf.n_requests-r0} req, {got/1e6:.0f} MB fetched, over-fetch x{got/(nnz*16):.2f}, {got/1e6/dt:.1f} MB/s")
