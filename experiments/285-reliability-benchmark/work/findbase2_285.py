"""285 -- locate X/data and X/indices by byte pattern, with 8-BYTE-ALIGNED probes, and pin the exact base
against the first 8 values of row 0 read through h5py (data [1,2,1,1,1,3,1,1]; indices [15375,1393,13172,
15891,200,1375,10264,9695]) and the obs `total_counts` / `n_genes_by_counts` identity on scattered rows."""
import json, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
WORK = Path("~/Code/virtual-cell-challenge-2026/experiments/285-reliability-benchmark/work")
B = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
NAME = "D2_Stim48hr.assigned_guide.h5ad"; STEM = NAME.split(".")[0]; url = f"{B}/{NAME}"
SIZE = 137554355000; NG = 18130
ip = np.load(WORK / f"indptr_{STEM}.npy")
z = np.load(WORK / f"probe2_{STEM}.npz", allow_pickle=True)
tot = np.asarray(z["total"], dtype=np.float64); ngc = np.asarray(z["ngc"], dtype=np.int64)
D0 = np.array([1,2,1,1,1,3,1,1], dtype=np.float32); I0 = np.array([15375,1393,13172,15891,200,1375,10264,9695], dtype=np.int64)
def get(lo, hi):
    req = urllib.request.Request(url, headers={"Range": f"bytes={lo}-{hi-1}"})
    for a in range(5):
        try:
            with urllib.request.urlopen(req, timeout=180) as r: return r.read()
        except Exception:
            if a == 4: raise
            time.sleep(1)
def classify(off):
    off = (off // 8) * 8
    raw = get(off, off + 8192)
    i8 = np.frombuffer(raw, dtype="<i8"); f4 = np.frombuffer(raw, dtype="<f4")
    return off, float(((i8 >= 0) & (i8 < NG)).mean()), float(((f4 > 0) & (f4 < 1e5) & (f4 == np.floor(f4))).mean())
t0 = time.time()
grid = [(int(x)//8)*8 for x in np.linspace(0, SIZE - 8192, 80)]
with ThreadPoolExecutor(max_workers=16) as ex: prof = list(ex.map(classify, grid))
for off, a, b in prof: print(f"  {off/1e9:8.2f} GB: index-like {a:.2f}  data-like {b:.2f}", flush=True)
print(f"[grid {time.time()-t0:.0f}s]", flush=True)
out = {}
for nm, key, itemsize in (("data", 2, 4), ("indices", 1, 8)):
    hits = [o for o, a, b in prof if (b if nm == "data" else a) > 0.9]
    if not hits: print(f"{nm}: NO grid hit"); continue
    step = grid[1] - grid[0]
    lo, hi = max(0, min(hits) - step), min(hits)
    pred = (lambda o: classify(o)[2] > 0.9) if nm == "data" else (lambda o: classify(o)[1] > 0.9)
    while hi - lo > 65536:
        mid = ((lo + hi) // 2 // 8) * 8
        if pred(mid): hi = mid
        else: lo = mid
    print(f"{nm}: region starts in [{lo:,}, {hi:,}]", flush=True)
    win_lo = max(0, lo - (1 << 20)); win = get(win_lo, hi + (1 << 20))
    ref = D0.tobytes() if nm == "data" else I0.tobytes()
    p = win.find(ref)
    base = win_lo + p if p >= 0 else None
    print(f"{nm}: row-0 prefix found at {base:,}" if base is not None else f"{nm}: row-0 prefix NOT found", flush=True)
    out[f"base_{nm}"] = base
# validate on 60 scattered rows with the total_counts identity
if out.get("base_data") is not None and out.get("base_indices") is not None:
    bd, bi = out["base_data"], out["base_indices"]
    rng = np.random.default_rng(20260906); rows = np.sort(rng.choice(ngc.size, 60, replace=False))
    def check(r):
        a, b = int(ip[r]), int(ip[r+1])
        d = np.frombuffer(get(bd + a*4, bd + b*4), dtype="<f4")
        i = np.frombuffer(get(bi + a*8, bi + b*8), dtype="<i8")
        return bool(d.size == ngc[r] == i.size and (d > 0).all() and float(d.sum()) == tot[r]
                    and (np.diff(i) > 0).all() and i[0] >= 0 and i[-1] < NG)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=16) as ex: res = list(ex.map(check, rows))
    nb = sum(1 for x in res if not x)
    print(f"KA-layout: {len(res)-nb}/{len(res)} scattered rows decode exactly ({time.time()-t0:.0f}s) -> "
          f"{'PASS' if nb == 0 else 'FAIL'}", flush=True)
    out["ka_layout"] = dict(checked=len(res), bad=nb, ok=bool(nb == 0))
json.dump(out, open(WORK / f"findbase2_{STEM}.json", "w"), indent=1, default=int)
