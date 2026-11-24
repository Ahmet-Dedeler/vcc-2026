import sys, numpy as np
from pathlib import Path
ROOT = Path("~/Code/virtual-cell-challenge-2026"); HERE = ROOT/"experiments/299-reliable-training-rows"
sys.path.insert(0, str(ROOT/"experiments/272-noise-corrected-operator"))
import op272 as O
for line in ("k562_gwps","orion_hct116","orion_hek293t"):
    z = np.load(HERE/f"work/rel/{line}.npz", allow_pickle=True)
    names=[str(x) for x in z["names"]]; rnp=z["r_nopanel"]
    z271 = np.load(ROOT/f"experiments/271-donor-reliability-estimation/work/half/rows-{line}.npz", allow_pickle=True)
    t271=[str(t) for t in z271["targets"]]; ok=np.asarray(z271["avail_A"],bool)&np.asarray(z271["avail_B"],bool)
    r271=O.row_cos(z271["A"].astype(np.float64), z271["B"].astype(np.float64))
    pos={n:i for i,n in enumerate(names)}
    sh=[(pos[n],t271.index(n)) for n in names if n in t271 and ok[t271.index(n)]]
    m=np.array([rnp[i] for i,_ in sh]); t=np.array([r271[k] for _,k in sh]); d=np.abs(m-t)
    print(f"{line}: n {len(sh)} max|d| {d.max():.3e} median|d| {np.median(d):.3e}")
