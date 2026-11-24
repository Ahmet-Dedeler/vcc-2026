"""425b robustness: dose scan of the MDE-structured magnitude rules, exponent scan, and a split-half over targets
(fit c on one half of the gated targets, read the other half) to check the H1 optimum is not a per-target artefact."""
import json
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
W323 = ROOT / "experiments/323-current-family-validation/work"
r = np.load(W323 / "reference.npz", allow_pickle=True)
s = np.load(ROOT / "experiments/410-pertscale-revisit/work/h1/stats-base395.npz", allow_pickle=True)
genes = r["genes"].astype(str)[r["keep"].astype(bool)]
targets = s["targets"].astype(str)
L = np.nan_to_num(np.asarray(s["lfc"], float))
Q = np.asarray(s["q"], float)
lin = np.asarray(r["ref_lin"], float)
count = lin * np.median(np.asarray(r["libraries"], float)) / 1e6
ref = pl.read_parquet(W323 / "reference_de.parquet")
gi = {g: i for i, g in enumerate(genes)}
gates = []
for k, t in enumerate(targets):
    m = ref.filter(pl.col("target") == t)
    g = m["feature"].to_numpy().astype(str); q = m["p_adj"].to_numpy().astype(float)
    l = m["log2_fold_change"].to_numpy().astype(float)
    ok = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if ok.sum() >= 10:
        gates.append((k, np.array([gi[x] for x in g[ok]]), l[ok]))


def nmae_t(Lhat, sel=None):
    v = []
    for n, (k, idx, l) in enumerate(gates):
        if sel is not None and n not in sel:
            continue
        v.append(np.abs(Lhat[k, idx] - l).sum() / np.abs(l).sum())
    return float(np.mean(v))


sg = np.sign(L)
base = nmae_t(L)
out = {"base": base, "add": {}, "floor": {}, "expo": {}, "split": {}}
print(f"base {base:.5f}")
for c in (0.05, 0.075, 0.1, 0.125, 0.15):
    v = nmae_t(L + sg * c / np.sqrt(np.maximum(count, 1e-3))[None, :]); out["add"][c] = v
    print(f"  add   c {c:5.3f}: {v:.5f}  Δ {v - base:+.5f}")
for c in (0.15, 0.2, 0.25, 0.3):
    v = nmae_t(sg * np.maximum(np.abs(L), c / np.sqrt(np.maximum(count, 1e-3))[None, :])); out["floor"][c] = v
    print(f"  floor c {c:5.3f}: {v:.5f}  Δ {v - base:+.5f}")
for e in (0.0, 0.25, 0.5, 0.75):
    best = None
    for c in np.linspace(0.02, 0.3, 29):
        v = nmae_t(L + sg * c / np.maximum(count, 1e-3)[None, :] ** e)
        if best is None or v < best[1]:
            best = (c, v)
    out["expo"][e] = best
    print(f"  add exponent {e}: best c {best[0]:.3f} -> {best[1]:.5f}  Δ {best[1] - base:+.5f}")
rng = np.random.default_rng(0)
res = []
for rep in range(20):
    perm = rng.permutation(len(gates)); A, B = set(perm[: len(gates) // 2]), set(perm[len(gates) // 2:])
    cs = np.linspace(0.02, 0.3, 29)
    fit = min(cs, key=lambda c: nmae_t(L + sg * c / np.sqrt(np.maximum(count, 1e-3))[None, :], A))
    gainB = nmae_t(L + sg * fit / np.sqrt(np.maximum(count, 1e-3))[None, :], B) - nmae_t(L, B)
    res.append((fit, gainB))
res = np.array(res)
out["split"] = {"c_fit_mean": res[:, 0].mean(), "c_fit_sd": res[:, 0].std(), "heldout_gain_mean": res[:, 1].mean(),
                "heldout_gain_sd": res[:, 1].std(), "heldout_gain_max": res[:, 1].max()}
print(f"  split-half (20 reps): c fit {res[:, 0].mean():.3f} ± {res[:, 0].std():.3f}; held-out Δnmae {res[:, 1].mean():+.5f} ± {res[:, 1].std():.5f} (worst {res[:, 1].max():+.5f})")
# per-target: how many targets improve
v0 = np.array([np.abs(L[k, idx] - l).sum() / np.abs(l).sum() for k, idx, l in gates])
L1 = L + sg * 0.1 / np.sqrt(np.maximum(count, 1e-3))[None, :]
v1 = np.array([np.abs(L1[k, idx] - l).sum() / np.abs(l).sum() for k, idx, l in gates])
print(f"  per-target (c 0.1): improves {np.sum(v1 < v0)} / {len(gates)}; median Δ {np.median(v1 - v0):+.5f}")
(HERE / "work/mde425b.json").write_text(json.dumps(out, indent=1, default=float))
