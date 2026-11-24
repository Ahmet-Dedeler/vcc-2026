"""427f: reach under alternative oracle-free ORDERS of our calls (signs = our delivered lfc, unchanged).
Orders (called genes first, as the scorer does): delivered (q, p); |our lfc| desc; |donor avg| desc; donor agreement
count then |donor avg|; our p but demoting genes whose sign disagrees with the donor average."""
from pathlib import Path
import numpy as np, polars as pl
ROOT = Path(__file__).resolve().parents[2]
W323 = ROOT / "experiments/323-current-family-validation/work"
r = np.load(W323 / "reference.npz", allow_pickle=True)
genes = r["genes"].astype(str)[r["keep"].astype(bool)]; gi = {g: i for i, g in enumerate(genes)}
s = np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz"); targets = s["targets"].astype(str)
L = np.nan_to_num(s["lfc"]); Q = s["q"]; P = s["p"]
DON = ["k562_gwps", "k562_essential", "rpe1", "nadig_jurkat", "orion_hct116", "orion_hek293t"]
D = np.full((len(DON), len(targets), len(genes)), np.nan, dtype=np.float32)
for di, d in enumerate(DON):
    bd = ROOT / "data/donors" / d
    Pn = [l.strip() for l in open(bd / "perts.txt")]; G = [l.strip() for l in open(bd / "genes.txt")]
    pi = {p: i for i, p in enumerate(Pn)}; gcol = {g: i for i, g in enumerate(G)}
    M = np.load(bd / "delta.npy", mmap_mode="r")
    cols = np.array([gcol.get(g, -1) for g in genes]); okc = cols >= 0
    for k, t in enumerate(targets):
        if t in pi:
            D[di, k, okc] = np.asarray(M[pi[t]], dtype=np.float32)[cols[okc]]
with np.errstate(invalid="ignore"):
    avg = np.nan_to_num(np.nanmean(D, 0))
agree = np.nansum(np.sign(D) == np.sign(L)[None], 0)
ref = pl.read_parquet(W323 / "reference_de.parquet")
pools = []
for k, t in enumerate(targets):
    m = ref.filter(pl.col("target") == t)
    g = m["feature"].to_numpy().astype(str); q = m["p_adj"].to_numpy().astype(float); l = m["log2_fold_change"].to_numpy().astype(float)
    ok = (q < .05) & (g != t)
    if ok.sum():
        pools.append((k, np.array([gi[x] for x in g[ok]]), l[ok]))
def reach(score):
    """score: higher = earlier, applied within the called block and within the uncalled block separately."""
    vals = []
    for k, idx, lr in pools:
        called = Q[k, idx] < .05
        order = np.lexsort((-score[k, idx], ~called))
        adj = np.isfinite(lr) & (lr != 0)
        c = n = best = 0
        for j in order:
            if not adj[j]:
                continue
            n += 1; c += int(np.sign(L[k, idx[j]]) == np.sign(lr[j]))
            if c / n >= 0.9:
                best = n
        vals.append(best / idx.size)
    return float(np.mean(vals))
base_score = -np.log(np.maximum(P, 1e-300))
print(f"delivered order (by p): {reach(base_score):.5f}")
print(f"by |our lfc|:            {reach(np.abs(L)):.5f}")
print(f"by |donor avg|:          {reach(np.abs(avg)):.5f}")
print(f"by donor agreement, then |avg|: {reach(agree * 10 + np.abs(avg) / (np.abs(avg).max() + 1e-9)):.5f}")
dis = np.sign(avg) != np.sign(L)
print(f"by p, demote sign-disagreeing with donor avg: {reach(base_score - 1e6 * dis):.5f}")
print(f"by p, demote only if |avg| > 0.05 and disagree:  {reach(base_score - 1e6 * (dis & (np.abs(avg) > 0.05))):.5f}")
