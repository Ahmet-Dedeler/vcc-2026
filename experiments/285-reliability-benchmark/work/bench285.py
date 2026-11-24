"""285 -- marginal cost of a SCATTERED CSR row read through h5py over HTTP range requests (the only route
left after the chunk layout turned out to be interleaved, not linear)."""
import sys, time
from pathlib import Path
import h5py, numpy as np
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/280-gwcd4i-donor-replicates/work")
from remote_h5 import RemoteFile
WORK = Path("~/Code/virtual-cell-challenge-2026/experiments/285-reliability-benchmark/work")
B = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
rf = RemoteFile(f"{B}/D2_Stim48hr.assigned_guide.h5ad", block=1 << 20, max_blocks=512)
f = h5py.File(rf, "r"); D, I = f["X/data"], f["X/indices"]
ip = np.load(WORK / "indptr_D2_Stim48hr.npy")
plan = np.load(WORK / "plan_D2_Stim48hr.npy")
z = np.load(WORK / "probe2_D2_Stim48hr.npz", allow_pickle=True)
tot = np.asarray(z["total"], dtype=np.float64); ngc = np.asarray(z["ngc"], dtype=np.int64)
for tag, rows in (("warm-up", plan[:3]), ("30 in file order", plan[1000:1030]), ("30 scattered", plan[::len(plan)//30][:30])):
    t0 = time.time(); r0, b0 = rf.n_requests, rf.n_bytes; ok = 0
    for r in rows:
        a, b = int(ip[r]), int(ip[r + 1])
        d = D[a:b]; i = I[a:b]
        ok += int(d.size == ngc[r] and float(d.sum()) == tot[r] and (np.diff(i) > 0).all())
    dt = time.time() - t0
    print(f"{tag}: {len(rows)} rows in {dt:.1f}s = {dt/len(rows)*1000:.0f} ms/row, {rf.n_requests-r0} req, "
          f"{(rf.n_bytes-b0)/1e6:.0f} MB, KA {ok}/{len(rows)}; projected for 62,842 rows: "
          f"{62842*dt/len(rows)/60:.0f} min single-process", flush=True)
