"""425c: where does the analytic MDE gain live? add sign*c/sqrt(count) (c 0.05, 0.1) restricted to one control-count bin
at a time; and the same with a count floor (genes >= x counts/cell only)."""
from pathlib import Path
import numpy as np, polars as pl
ROOT = Path(__file__).resolve().parents[2]
W323 = ROOT / "experiments/323-current-family-validation/work"
r = np.load(W323 / "reference.npz", allow_pickle=True)
s = np.load(ROOT / "experiments/410-pertscale-revisit/work/h1/stats-base395.npz", allow_pickle=True)
genes = r["genes"].astype(str)[r["keep"].astype(bool)]; targets = s["targets"].astype(str)
L = np.nan_to_num(np.asarray(s["lfc"], float))
count = np.asarray(r["ref_lin"], float) * np.median(np.asarray(r["libraries"], float)) / 1e6
ref = pl.read_parquet(W323 / "reference_de.parquet"); gi = {g: i for i, g in enumerate(genes)}
gates = []
for k, t in enumerate(targets):
    m = ref.filter(pl.col("target") == t)
    g = m["feature"].to_numpy().astype(str); q = m["p_adj"].to_numpy().astype(float); l = m["log2_fold_change"].to_numpy().astype(float)
    ok = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if ok.sum() >= 10:
        gates.append((k, np.array([gi[x] for x in g[ok]]), l[ok]))
def nmae(Lh):
    return float(np.mean([np.abs(Lh[k, idx] - l).sum() / np.abs(l).sum() for k, idx, l in gates]))
base = nmae(L); sg = np.sign(L); add = 1 / np.sqrt(np.maximum(count, 1e-3))
edges = [0, 0.3, 0.6, 1.0, 1.5, 2.25, 3.5, 6, 1e9]
print(f"base {base:.5f}")
for c in (0.05, 0.1):
    print(f"c {c}: per-bin Δnmae (add only in that bin)")
    for a, b in zip(edges[:-1], edges[1:]):
        m = (count >= a) & (count < b)
        print(f"   {a:5.2f}-{b:8.2f} ({m.sum():5d} genes): {nmae(L + sg * c * add * m[None, :]) - base:+.5f}")
    for fl in (0.3, 0.6, 1.0, 1.5):
        m = count >= fl
        print(f"   floor count >= {fl}: {nmae(L + sg * c * add * m[None, :]) - base:+.5f}")
    for cap in (0.05, 0.1, 0.15):
        print(f"   add capped at {cap} log2: {nmae(L + sg * np.minimum(c * add, cap)[None, :]) - base:+.5f}")
