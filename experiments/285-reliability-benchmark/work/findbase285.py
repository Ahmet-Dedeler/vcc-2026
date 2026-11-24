"""285 -- locate X/data and X/indices in the CD4 object by BYTE PATTERN, then pin the exact base with a
known-answer row.  h5py's chunk index is unusable here (1.57 M chunks: `get_chunk_info(0)` alone did not
return in 10 min over HTTP), but the datasets are uncompressed, so their bytes are directly recognisable:
`indices` are int64 in [0, 18130) (six zero high bytes), `data` are float32 positive integers.
The exact base is then fixed by requiring a LANDMARK ROW to decode with exactly n_genes_by_counts[r]
entries summing to obs/total_counts[r] -- 280's KA-0 identity."""
import json, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
WORK = Path("~/Code/virtual-cell-challenge-2026/experiments/285-reliability-benchmark/work")
B = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
NAME = sys.argv[1] if len(sys.argv) > 1 else "D2_Stim48hr.assigned_guide.h5ad"
STEM = NAME.split(".")[0]; url = f"{B}/{NAME}"
SIZE = 137554355000; NNZ = 11417773379; NG = 18130
ip = np.load(WORK / f"indptr_{STEM}.npy")
z = np.load(WORK / f"probe2_{STEM}.npz", allow_pickle=True)
tot = np.asarray(z["total"], dtype=np.float64); ngc = np.asarray(z["ngc"], dtype=np.int64)
def get(lo, hi):
    req = urllib.request.Request(url, headers={"Range": f"bytes={lo}-{hi-1}"})
    for a in range(5):
        try:
            with urllib.request.urlopen(req, timeout=180) as r: return r.read()
        except Exception:
            if a == 4: raise
            time.sleep(1)
def classify(off):
    raw = get(off, off + 8192)
    i8 = np.frombuffer(raw, dtype="<i8"); f4 = np.frombuffer(raw, dtype="<f4")
    idxish = float(((i8 >= 0) & (i8 < NG)).mean())
    datish = float(((f4 > 0) & (f4 < 1e5) & (f4 == np.floor(f4))).mean())
    return off, idxish, datish
t0 = time.time()
grid = [int(x) for x in np.linspace(0, SIZE - 8192, 60)]
with ThreadPoolExecutor(max_workers=12) as ex:
    prof = list(ex.map(classify, grid))
for off, a, b in prof:
    print(f"  {off/1e9:8.2f} GB: index-like {a:.2f}  data-like {b:.2f}", flush=True)
print(f"[grid {time.time()-t0:.0f}s]", flush=True)
# refine the data/indices boundary by bisection on the classifier
def is_idx(off): return classify(off)[1] > 0.9
def is_dat(off): return classify(off)[2] > 0.9
regions = {}
for nm, pred in (("data", is_dat), ("indices", is_idx)):
    hits = [off for off, a, b in prof if (b > 0.9 if nm == "data" else a > 0.9)]
    if not hits: print(f"{nm}: no grid hit"); continue
    lo_hit = min(hits)
    lo, hi = max(0, lo_hit - (grid[1] - grid[0])), lo_hit
    while hi - lo > 4096:
        mid = (lo + hi) // 2
        if pred(mid): hi = mid
        else: lo = mid
    regions[nm] = hi
    print(f"{nm}: starts within [{lo:,}, {hi:,}] ({hi/1e9:.3f} GB)", flush=True)
# pin the exact base with a landmark row
out = {"regions": regions}
for nm, itemsize, dt in (("data", 4, "<f4"), ("indices", 8, "<i8")):
    if nm not in regions: continue
    approx = regions[nm]
    win_lo = max(0, approx - (1 << 21)); win = get(win_lo, approx + (1 << 21))
    r = 5000
    a, b = int(ip[r]), int(ip[r + 1]); n = b - a
    found = None
    for o in range(0, len(win) - itemsize * n, itemsize):
        cand = win_lo + o - a * itemsize
        if cand < 0: continue
        pos = o
        if nm == "data":
            v = np.frombuffer(win[pos:pos + 4 * ngc[0]], dtype=dt)
            if v.size == ngc[0] and (v > 0).all() and float(v.sum()) == tot[0]:
                found = win_lo + o; break
        else:
            v = np.frombuffer(win[pos:pos + 8 * ngc[0]], dtype=dt)
            if v.size == ngc[0] and (np.diff(v) > 0).all() and v[0] >= 0 and v[-1] < NG:
                found = win_lo + o; break
    print(f"{nm}: base candidate {found:,}" if found else f"{nm}: base NOT pinned in the window", flush=True)
    out[f"base_{nm}"] = found
json.dump(out, open(WORK / f"findbase_{STEM}.json", "w"), indent=1, default=int)
