"""427e: what nmae does a plain DONOR-AVERAGE log2FC get on H1? l_hat_{p,g} = w * mean_d(delta_{d,p,g}) / ln 2 over the six
non-H1 donors covering (p, g) (0 where none), for a global w. Compared with our delivered lfc (0.952). If simple donor
magnitudes beat ours, our DE magnitudes (not signs) are the defect."""
from pathlib import Path
import numpy as np, polars as pl
ROOT = Path(__file__).resolve().parents[2]
W323 = ROOT / "experiments/323-current-family-validation/work"
r = np.load(W323 / "reference.npz", allow_pickle=True)
genes = r["genes"].astype(str)[r["keep"].astype(bool)]; gi = {g: i for i, g in enumerate(genes)}
s = np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz"); targets = s["targets"].astype(str)
L = np.nan_to_num(s["lfc"])
DON = ["k562_gwps", "k562_essential", "rpe1", "nadig_jurkat", "orion_hct116", "orion_hek293t"]
D = np.full((len(DON), len(targets), len(genes)), np.nan, dtype=np.float32)
for di, d in enumerate(DON):
    bd = ROOT / "data/donors" / d
    P = [l.strip() for l in open(bd / "perts.txt")]; G = [l.strip() for l in open(bd / "genes.txt")]
    pi = {p: i for i, p in enumerate(P)}; gcol = {g: i for i, g in enumerate(G)}
    M = np.load(bd / "delta.npy", mmap_mode="r")
    cols = np.array([gcol.get(g, -1) for g in genes]); okc = cols >= 0
    for k, t in enumerate(targets):
        if t in pi:
            D[di, k, okc] = np.asarray(M[pi[t]], dtype=np.float32)[cols[okc]]
ref = pl.read_parquet(W323 / "reference_de.parquet")
gates = []
for k, t in enumerate(targets):
    m = ref.filter(pl.col("target") == t)
    g = m["feature"].to_numpy().astype(str); q = m["p_adj"].to_numpy().astype(float); l = m["log2_fold_change"].to_numpy().astype(float)
    ok = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if ok.sum() >= 10:
        gates.append((k, np.array([gi[x] for x in g[ok]]), l[ok]))
def nmae(Lh):
    return float(np.mean([np.abs(Lh[k, idx] - l).sum() / np.abs(l).sum() for k, idx, l in gates]))
with np.errstate(invalid="ignore"):
    avg = np.nanmean(D, 0) / np.log(2)
avg = np.nan_to_num(avg)
K = np.concatenate([np.full(i.size, k) for k, i, _ in gates]); I = np.concatenate([i for _, i, _ in gates]); LL = np.concatenate([l for _, _, l in gates])
print(f"ours: nmae {nmae(L):.5f}; sign acc on gate {np.mean(np.sign(L[K, I]) == np.sign(LL)):.4f}; median |ours|/|l| {np.median(np.abs(L[K, I]) / np.abs(LL)):.3f}")
cov = avg[K, I] != 0
print(f"donor avg: sign acc on gate {np.mean(np.sign(avg[K, I][cov]) == np.sign(LL[cov])):.4f} (covered {cov.mean():.3f}); median |avg|/|l| {np.median(np.abs(avg[K, I][cov]) / np.abs(LL[cov])):.3f}; corr(|avg|,|l|) {np.corrcoef(np.abs(avg[K, I][cov]), np.abs(LL[cov]))[0, 1]:+.3f}; corr(|ours|,|l|) {np.corrcoef(np.abs(L[K, I]), np.abs(LL))[0, 1]:+.3f}")
for w in (0.25, 0.5, 1.0, 1.5, 2.0, 3.0):
    print(f"  donor avg x{w}: nmae {nmae(w * avg):.5f}   | blend 0.5*ours + 0.5*donor x{w}: {nmae(0.5 * L + 0.5 * w * avg):.5f}")

count = r["ref_lin"] * np.median(r["libraries"]) / 1e6
add = np.where(count >= 0.3, np.minimum(0.1 / np.sqrt(np.maximum(count, 1e-3)), 0.1), 0.0)[None, :]
base = nmae(L)
print(f"\nMDE (c0.1 cap0.1 min0.3) analytic, base {base:.5f}:")
print(f"  on our signs:                    {nmae(L + np.sign(L) * add) - base:+.5f}")
for wd in (0.5, 1.0, 1.5, 2.0):
    B = 0.5 * L + 0.5 * wd * avg
    print(f"  on blend(ours, donor x{wd}) signs+mag: {nmae(B + np.sign(B) * add) - base:+.5f}   (blend alone {nmae(B) - base:+.5f}); sign acc {np.mean(np.sign(B[K, I]) == np.sign(LL)):.4f}")
for wd in (1.0, 1.5):
    print(f"  on donor x{wd} alone + add:         {nmae(wd * avg + np.sign(avg) * add) - base:+.5f}")
