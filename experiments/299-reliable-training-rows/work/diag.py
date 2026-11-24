import sys, json, numpy as np
from pathlib import Path
ROOT = Path("~/Code/virtual-cell-challenge-2026"); HERE = ROOT/"experiments/299-reliable-training-rows"
sys.path.insert(0, str(ROOT/"experiments/162-retrodicting-validator")); sys.path.insert(0, str(ROOT/"experiments/272-noise-corrected-operator"))
import plolo, op272 as O
plolo.WORK = HERE/"work/plolo"
E062, OFFICIAL, PANEL = plolo.E062, plolo.OFFICIAL, list(plolo.PANEL)
line = sys.argv[1]
z = np.load(HERE/f"work/rel/{line}.npz", allow_pickle=True)
names = [str(x) for x in z["names"]]; r = z["r"]; cols = z["cols"]
z271 = np.load(ROOT/f"experiments/271-donor-reliability-estimation/work/half/rows-{line}.npz", allow_pickle=True)
t271 = [str(t) for t in z271["targets"]]
ok271 = np.asarray(z271["avail_A"],bool)&np.asarray(z271["avail_B"],bool)
r271 = O.row_cos(z271["A"].astype(np.float64), z271["B"].astype(np.float64))
pos = {n:i for i,n in enumerate(names)}
sh = [(pos[n], t271.index(n)) for n in names if n in t271 and ok271[t271.index(n)]]
mine = np.array([r[i] for i,_ in sh]); th = np.array([r271[k] for _,k in sh])
d = np.abs(mine-th)
print(f"{line}: n {len(sh)} max|d| {d.max():.3e} median|d| {np.median(d):.3e} p90 {np.percentile(d,90):.3e} corr {np.corrcoef(mine,th)[0,1]:.6f}")
k = int(np.argmax(d)); i,j = sh[k]
print("  worst row", names[i], "mine", mine[k], "271", th[k])
nA = np.linalg.norm(z271["A"].astype(np.float64), axis=1)
print("  its 271 normA", nA[j], "median normA", np.median(nA[ok271]))
# elementwise: rebuild the A row for this target from the half bank and compare
print("  n cols mine", cols.size, "n nonzero in 271 A row", int((z271['A'][j]!=0).sum()))
nz = np.flatnonzero(z271['A'][j]!=0)
print("  271 A nonzero positions inside cols?", np.isin(nz, cols).mean())
# ratio of cosines using only cols in 271's own arrays
r271c = O.row_cos(z271["A"].astype(np.float64)[:,cols], z271["B"].astype(np.float64)[:,cols])
print("  271 restricted-to-cols vs full:", float(np.nanmax(np.abs(r271c[ok271]-r271[ok271]))))
