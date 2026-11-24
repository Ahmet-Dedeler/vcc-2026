"""450j: does a NEAR line (same cell type, different line) carry the programme? H1 hESC truth; iPSC banks (kolf2, Tian 2019
iPSC, Feng iPSC VCC-panel subset) vs the same-line bank (Pan H1) and CLEAN4. P_down programme (396 rule), h1-test targets
excluded, common coverage with CLEAN4 and every bank."""
import runpy
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.argv = [sys.argv[0]]
g = runpy.run_path(str(HERE / "sameline450.py"), run_name="not_main")
P, lfc, gate, block, adj, truth, T, on_axis = (g[k] for k in ("P", "lfc", "gate", "block", "adj", "truth", "T", "on_axis"))


def pdown_bank(d, delta_name="delta.npy"):
    D = np.load(d / delta_name, mmap_mode="r")
    perts = [l.strip() for l in open(d / "perts.txt")] if (d / "perts.txt").exists() else None
    genes = [l.strip() for l in open(d / "genes.txt")]
    keep = np.ones(D.shape[0], bool) if perts is None else np.array([p not in set(T) for p in perts])
    X = np.asarray(D[np.where(keep)[0]], dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        pd = np.where((X != 0).sum(0) > 0, (X < 0).sum(0) / np.maximum((X != 0).sum(0), 1), np.nan)
    return on_axis(genes, pd), int(keep.sum())


banks = {
    "kolf2 iPSC (near line)": ROOT / "data/external/kolf2",
    "Tian 2019 iPSC (near line)": ROOT / "experiments/225-source-sweep/work/tables/tian2019_ipsc",
}
S = {}
for name, d in banks.items():
    v, n = pdown_bank(d)
    S[name] = (np.where(np.isnan(v), np.nan, np.where(v > 0.5, -1.0, 1.0)), n)
S["Pan H1 (same line)"] = (P["pan_pdown"], 2324)
S["CLEAN4 (4 cancer lines)"] = (P["clean4"], 0)
both = np.all([np.isfinite(s) for s, _ in S.values()], axis=0)
print(f"genes covered by all predictors: {both.sum()}")
for name, (s, n) in S.items():
    out = []
    for nm, sub in (("gate-block", gate & block), ("fid-like block", adj & block)):
        mm = sub & both[None, :]
        out.append(f"{nm} {((np.broadcast_to(s, lfc.shape) == truth) & mm).sum() / mm.sum():.4f}")
    print(f"  {name:30s} rows {n:5d}  " + "   ".join(out))
