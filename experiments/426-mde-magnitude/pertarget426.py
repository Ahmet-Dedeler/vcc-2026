"""426c analytic: scale the MDE add per target by a predicted effect strength s_p = (row norm of our delivered LFC /
median)^beta, i.e. strong perturbations get a larger add (their gate |lfc| sit further above the MDE). H1, 426 base stats."""
from pathlib import Path
import numpy as np, polars as pl
ROOT = Path(__file__).resolve().parents[2]
W323 = ROOT / "experiments/323-current-family-validation/work"
r = np.load(W323 / "reference.npz", allow_pickle=True)
genes = r["genes"].astype(str)[r["keep"].astype(bool)]
s = np.load(Path(__file__).parent / "work/h1/stats-base_s0_bps.npz"); targets = s["targets"].astype(str)
L = np.nan_to_num(s["lfc"]); Q = s["q"]
count = r["ref_lin"] * np.median(r["libraries"]) / 1e6
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
base = nmae(L); sg = np.sign(L)
norm = np.linalg.norm(L, axis=1); strength = norm / np.median(norm)
ncall = (Q < .05).sum(1); sc = ncall / np.median(ncall)
print(f"base {base:.5f}; corr(our row norm, real gate mean|l|) = "
      f"{np.corrcoef([norm[k] for k,_,_ in gates], [np.abs(l).mean() for _,_,l in gates])[0,1]:+.3f}; "
      f"corr(n calls, n_real) = {np.corrcoef([ncall[k] for k,_,_ in gates], [l.size for _,_,l in gates])[0,1]:+.3f}")
for c in (0.1,):
    for cap in (0.1, 0.2):
        for beta in (0.0, 0.5, 1.0):
            for lab, st in (("norm", strength), ("ncall", sc)):
                if beta == 0 and lab == "ncall":
                    continue
                add = np.where(count >= 0.3, np.minimum(c / np.sqrt(np.maximum(count, 1e-3))[None, :] * st[:, None] ** beta, cap), 0.0)
                print(f"  c {c} cap {cap} beta {beta} ({lab:5s}): {nmae(L + sg * add) - base:+.5f}")
