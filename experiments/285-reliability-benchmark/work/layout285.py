"""285 -- is X/data and X/indices of the CD4 cell-level object laid out LINEARLY in chunk order?

If it is, element i of an uncompressed chunked 1-D dataset sits at `base + i * itemsize`, so a CSR row r is the
byte range [base + indptr[r]*itemsize, base + indptr[r+1]*itemsize) and no chunk-index walk is needed.
280's `layout280.py` made the same check on the pseudobulk object.  Also benchmarks a parallel row fetch.
"""
from __future__ import annotations

import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/280-gwcd4i-donor-replicates/work")
from remote_h5 import RemoteFile  # noqa: E402

WORK = Path("~/Code/virtual-cell-challenge-2026/experiments/285-reliability-benchmark/work")
B = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
NAME = sys.argv[1] if len(sys.argv) > 1 else "D2_Stim48hr.assigned_guide.h5ad"
STEM = NAME.split(".")[0]
url = f"{B}/{NAME}"

rf = RemoteFile(url, block=1 << 22, max_blocks=32)
f = h5py.File(rf, "r")
ip = np.load(WORK / f"indptr_{STEM}.npy")
t0 = time.time()
bases = {}
for nm in ("X/data", "X/indices"):
    d = f[nm]
    dsid = d.id
    n = dsid.get_num_chunks()
    cb = d.chunks[0] * d.dtype.itemsize
    probe = (0, 1, 2, 1000, n // 2, n - 2, n - 1)
    infos = [dsid.get_chunk_info(i) for i in probe]
    base = infos[0].byte_offset
    ok = all(inf.byte_offset == base + i * cb for i, inf in zip(probe, infos))
    print(f"{nm}: {d.shape[0]:,} {d.dtype} in {n:,} chunks of {d.chunks[0]}, base {base:,}, chunk bytes {cb}, "
          f"LINEAR {ok}")
    for i, inf in zip(probe, infos):
        print(f"   chunk {i:>9}: offset {inf.byte_offset:>16,}  size {inf.size:>8,}  predicted "
              f"{base + i*cb:>16,}  delta {inf.byte_offset - (base + i*cb):>+14,}")
    bases[nm] = (base, ok, d.dtype.itemsize)
print(f"[{time.time()-t0:.0f}s, {rf.n_requests} req, {rf.n_bytes/1e6:.0f} MB]")

# verify against h5py's own read of a few rows, and benchmark a parallel fetch
plan = np.load(WORK / f"plan_{STEM}.npy")
bd, okd, isd = bases["X/data"]
bi, oki, isi = bases["X/indices"]


def get(lo, hi):
    req = urllib.request.Request(url, headers={"Range": f"bytes={lo}-{hi-1}"})
    for a in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception:
            if a == 3:
                raise


if okd and oki:
    rows = plan[np.linspace(0, plan.size - 1, 5).astype(int)]
    for r in rows:
        a, b = int(ip[r]), int(ip[r + 1])
        raw = get(bd + a * isd, bd + b * isd)
        mine = np.frombuffer(raw, dtype="<f4")
        theirs = f["X/data"][a:b]
        rawi = get(bi + a * isi, bi + b * isi)
        minei = np.frombuffer(rawi, dtype="<i8")
        theirsi = f["X/indices"][a:b]
        print(f"  row {r}: nnz {b-a}, data max|d| {np.abs(mine-theirs).max():.3g}, "
              f"indices max|d| {np.abs(minei-theirsi).max()}")
    t0 = time.time()
    sub = plan[:400]
    with ThreadPoolExecutor(max_workers=16) as ex:
        res = list(ex.map(lambda r: (get(bd + int(ip[r]) * isd, bd + int(ip[r + 1]) * isd),
                                     get(bi + int(ip[r]) * isi, bi + int(ip[r + 1]) * isi)), sub))
    dt = time.time() - t0
    got = sum(len(a) + len(b) for a, b in res)
    print(f"  parallel fetch: {len(sub)} rows, {got/1e6:.1f} MB in {dt:.1f}s with 16 threads -> "
          f"{len(sub)/dt:.0f} rows/s; {plan.size} rows would take {plan.size/(len(sub)/dt)/60:.1f} min")
