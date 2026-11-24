"""412 — oracle-free fold-change magnitude rules for nmae on H1 (the DE channel is free under ballast).
Base = the shipped 395 construction's delivered H1 statistics (410/work/h1/stats-base395.npz). nmae per the scorer:
gate = reference q < 0.05, finite nonzero lfc, own target removed, >= 10 genes; mean over gated targets.
Every rule reads only OUR delivered lfc / q / p (no reference information)."""
import numpy as np, polars as pl, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
W323 = ROOT / "experiments/323-current-family-validation/work"
r = np.load(W323 / "reference.npz", allow_pickle=True)
s = np.load(ROOT / "experiments/410-pertscale-revisit/work/h1/stats-base395.npz", allow_pickle=True)
genes = r["genes"].astype(str)[r["keep"].astype(bool)]
targets = s["targets"].astype(str)
ref = pl.read_parquet(W323 / "reference_de.parquet")
gi = {g: i for i, g in enumerate(genes)}
L = np.asarray(s["lfc"], float); Q = np.asarray(s["q"], float); P = np.asarray(s["p"], float)
gates = []
for k, t in enumerate(targets):
    m = ref.filter(pl.col("target") == t)
    g = m["feature"].to_numpy().astype(str); q = m["p_adj"].to_numpy().astype(float); l = m["log2_fold_change"].to_numpy().astype(float)
    ok = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if ok.sum() >= 10:
        gates.append((k, np.array([gi[x] for x in g[ok]]), l[ok]))
def nmae(Lhat):
    v = []
    for k, idx, l in gates:
        p = np.nan_to_num(Lhat[k, idx], nan=0.0)
        v.append(np.abs(p - l).sum() / np.abs(l).sum())
    return float(np.mean(v))
base = nmae(L)
print(f"gated targets {len(gates)}; base nmae {base:.5f} (410 H1 base395 row: 0.95210)")
out = {"base": base, "rules": {}}
def rec(name, Lh):
    v = nmae(Lh); out["rules"][name] = v; print(f"  {name:36s} nmae {v:.5f}  Δ {v - base:+.5f}")
for k in (1.25, 1.5, 2.0, 3.0):
    rec(f"global x{k}", L * k)
rank = np.argsort(np.argsort(np.nan_to_num(P, nan=1.0), axis=1), axis=1)   # 0 = our most confident gene per target
for N in (50, 100, 200, 400):
    for k in (2.0, 3.0, 5.0):
        rec(f"top{N} by p x{k}", np.where(rank < N, L * k, L))
for f in (0.1, 0.2, 0.3, 0.5):
    sig = Q < 0.05
    rec(f"floor |lfc|>={f} on our calls", np.where(sig & (np.abs(L) < f) & (L != 0), np.sign(L) * f, L))
(ROOT / "experiments/412-nmae-magnitude-rules/work/rules412.json").write_text(json.dumps(out, indent=1))
