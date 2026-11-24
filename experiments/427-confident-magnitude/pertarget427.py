"""427b: per-target oracle magnitude scale on H1, and what predicts it. For each gated target p: k_p* = argmin_k nmae_p(k *
ours) (grid), the target's sign accuracy on its gate, n_real, and oracle-free features: our row norm, the donors'
effect norms / agreement for p. If a feature predicts k_p*, a per-target scale is an oracle-free lever (388's per-target
oracle scalar was worth -0.023)."""
from pathlib import Path
import numpy as np, polars as pl
from scipy.stats import spearmanr
ROOT = Path(__file__).resolve().parents[2]
W323 = ROOT / "experiments/323-current-family-validation/work"
r = np.load(W323 / "reference.npz", allow_pickle=True)
genes = r["genes"].astype(str)[r["keep"].astype(bool)]; gi = {g: i for i, g in enumerate(genes)}
s = np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz"); targets = s["targets"].astype(str)
L = np.nan_to_num(s["lfc"])
ref = pl.read_parquet(W323 / "reference_de.parquet")
DON = ["k562_gwps", "orion_hct116", "orion_hek293t", "k562_essential", "rpe1", "nadig_jurkat"]
dnorm = {d: np.full(len(targets), np.nan) for d in DON}
for d in DON:
    bd = ROOT / "data/donors" / d
    P = [l.strip() for l in open(bd / "perts.txt")]; pi = {p: i for i, p in enumerate(P)}
    M = np.load(bd / "delta.npy", mmap_mode="r")
    for k, t in enumerate(targets):
        if t in pi:
            dnorm[d][k] = float(np.linalg.norm(np.asarray(M[pi[t]], dtype=np.float64)))
grid = np.array([0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8])
rows = []
for k, t in enumerate(targets):
    m = ref.filter(pl.col("target") == t)
    g = m["feature"].to_numpy().astype(str); q = m["p_adj"].to_numpy().astype(float); l = m["log2_fold_change"].to_numpy().astype(float)
    ok = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if ok.sum() < 10:
        continue
    idx = np.array([gi[x] for x in g[ok]]); ll = l[ok]; p_ = L[k, idx]
    errs = [np.abs(kk * p_ - ll).sum() / np.abs(ll).sum() for kk in grid]
    j = int(np.argmin(errs))
    acc = np.mean(np.sign(p_) == np.sign(ll))
    rows.append(dict(t=t, n_real=int(ok.sum()), kstar=grid[j], gain=errs[j] - errs[2], acc=acc, ournorm=float(np.linalg.norm(L[k])),
                     **{f"n_{d}": dnorm[d][k] for d in DON}))
import pandas as pd
df = pd.DataFrame(rows)
print(f"gated targets {len(df)}; oracle per-target k: mean gain {df.gain.mean():+.5f}; k* distribution {df.kstar.value_counts().sort_index().to_dict()}")
for c in ["n_real", "acc", "ournorm"] + [f"n_{d}" for d in DON]:
    v = df[c].to_numpy(); okv = np.isfinite(v)
    rho = spearmanr(v[okv], df.kstar.to_numpy()[okv]).correlation
    print(f"  spearman(k*, {c:18s}) = {rho:+.3f}  (n {okv.sum()})")
print(df.sort_values("kstar")[["t", "n_real", "kstar", "gain", "acc", "ournorm", "n_orion_hct116", "n_k562_gwps"]].to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# 427c: equalise delivered magnitude across targets: lfc_p *= (median norm / norm_p)^beta (oracle-free)
gates = []
for k, t in enumerate(targets):
    m = ref.filter(pl.col("target") == t)
    g = m["feature"].to_numpy().astype(str); q = m["p_adj"].to_numpy().astype(float); l = m["log2_fold_change"].to_numpy().astype(float)
    ok = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if ok.sum() >= 10:
        gates.append((k, np.array([gi[x] for x in g[ok]]), l[ok]))
def nmae(Lh):
    return float(np.mean([np.abs(Lh[k, idx] - l).sum() / np.abs(l).sum() for k, idx, l in gates]))
base = nmae(L)
nrm = np.linalg.norm(L, axis=1); med = np.median(nrm)
print(f"\nbase nmae {base:.5f}")
for beta in (0.25, 0.5, 0.75, 1.0):
    for glob in (1.0, 1.25):
        sc = (med / np.maximum(nrm, 1e-9)) ** beta * glob
        print(f"  equalise beta {beta:4.2f} x{glob}: {nmae(L * sc[:, None]) - base:+.5f}   (scale range {sc.min():.2f}..{sc.max():.2f})")

# 427d: does donor CONSISTENCY for target p (mean pairwise corr of donor rows on shared genes) predict its sign accuracy / k*?
from itertools import combinations
rowsD = {}
for d in DON:
    bd = ROOT / "data/donors" / d
    P = [l.strip() for l in open(bd / "perts.txt")]; G = [l.strip() for l in open(bd / "genes.txt")]
    pi = {p: i for i, p in enumerate(P)}; M = np.load(bd / "delta.npy", mmap_mode="r")
    rowsD[d] = (pi, {g: i for i, g in enumerate(G)}, M)
cons = []
for k, t in enumerate(targets):
    vecs = []
    for d in DON:
        pi, gcol, M = rowsD[d]
        if t in pi:
            vecs.append((d, gcol, np.asarray(M[pi[t]], dtype=np.float64)))
    cs = []
    for (d1, g1, v1), (d2, g2, v2) in combinations(vecs, 2):
        common = [g for g in g1 if g in g2]
        a1 = v1[[g1[g] for g in common]]; a2 = v2[[g2[g] for g in common]]
        if a1.std() > 0 and a2.std() > 0:
            cs.append(np.corrcoef(a1, a2)[0, 1])
    cons.append(np.mean(cs) if cs else np.nan)
cons = np.array(cons)
df["cons"] = [cons[list(targets).index(t)] for t in df.t]
okc = np.isfinite(df.cons.to_numpy())
print(f"\ndonor consistency: spearman with acc {spearmanr(df.cons[okc], df.acc[okc]).correlation:+.3f}, with k* {spearmanr(df.cons[okc], df.kstar[okc]).correlation:+.3f}, with gain {spearmanr(df.cons[okc], df.gain[okc]).correlation:+.3f} (n {okc.sum()})")
# rule: scale by consistency tercile
for lo_s, hi_s in ((1.0, 1.5), (1.0, 2.0), (0.75, 1.5), (0.75, 2.0)):
    q1, q2 = np.nanquantile(cons, [1/3, 2/3])
    sc = np.where(np.isfinite(cons), np.where(cons >= q2, hi_s, np.where(cons >= q1, (lo_s + hi_s) / 2, lo_s)), 1.0)
    print(f"  scale by consistency tercile (lo {lo_s}, hi {hi_s}): {nmae(L * sc[:, None]) - base:+.5f}")
