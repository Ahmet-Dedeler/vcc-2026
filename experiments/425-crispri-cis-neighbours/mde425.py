"""425b: nmae is scored only on the reference's significant genes, so the L1-optimal prediction for gene g is the median
of l_{p,g} CONDITIONAL on g being significant. Conditional on significance, |l| is bounded below by the minimum detectable
effect at g's control expression (MDE ~ 1/sqrt(mean count)). With sign accuracy a, the L1 optimum sits at the
(a - 1/2)/a quantile of |l| given significance: close to the MDE when a ~ 0.56. Test on H1 (395 delivered stats):
magnitudes by expression, signs ours.
"""
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
L = np.asarray(s["lfc"], float)
lin = np.asarray(r["ref_lin"], float)                 # control arithmetic-mean CPM per kept gene (checked below)
libs = np.asarray(r["libraries"], float)
print("ref_lin: min %.3g median %.3g max %.3g ; sum %.4g" % (lin.min(), np.median(lin), lin.max(), lin.sum()))
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


def nmae(Lhat):
    v = []
    for k, idx, l in gates:
        p = np.nan_to_num(Lhat[k, idx], nan=0.0)
        v.append(np.abs(p - l).sum() / np.abs(l).sum())
    return float(np.mean(v))


base = nmae(L)
print(f"gated {len(gates)}  base nmae {base:.5f}")
# |l| of gate genes vs expression; our sign accuracy vs expression
mean_count = lin * np.median(libs) / 1e6            # expected counts per control cell at median depth
allk = np.concatenate([np.full(i.size, k) for k, i, _ in gates]); alli = np.concatenate([i for _, i, _ in gates])
alll = np.concatenate([l for _, _, l in gates])
ours = L[allk, alli]
edges = np.quantile(mean_count[alli], np.linspace(0, 1, 9))
print("\nexpr bin (mean count/cell)   n    med|l|  q10|l|   med|ours|  sign_acc  ours/|l| med")
for a, b in zip(edges[:-1], edges[1:]):
    m = (mean_count[alli] >= a) & (mean_count[alli] <= b)
    sa = np.mean(np.sign(ours[m]) == np.sign(alll[m]))
    print(f"  {a:8.3f} - {b:8.3f}  {m.sum():6d}  {np.median(np.abs(alll[m])):.3f}   {np.quantile(np.abs(alll[m]), .1):.3f}"
          f"    {np.median(np.abs(np.nan_to_num(ours[m]))):.3f}     {sa:.3f}   {np.median(np.abs(np.nan_to_num(ours[m])) / np.abs(alll[m])):.3f}")
out = {"base": base, "rows": {}}


def rec(name, Lh):
    v = nmae(Lh); out["rows"][name] = v; print(f"  {name:52s} nmae {v:.5f}  Δ {v - base:+.5f}")


sg = np.sign(np.nan_to_num(L))
inv = 1.0 / np.sqrt(np.maximum(mean_count, 1e-3))
print()
for c in (0.05, 0.1, 0.15, 0.2, 0.3):
    rec(f"sign(ours) * {c}/sqrt(count)", sg * c * inv[None, :])
for c in (0.05, 0.1, 0.15, 0.2):
    rec(f"sign(ours) * max(|ours|, {c}/sqrt(count))", sg * np.maximum(np.abs(np.nan_to_num(L)), c * inv[None, :]))
for c in (0.05, 0.1, 0.2):
    rec(f"ours + sign(ours)*{c}/sqrt(count)", np.nan_to_num(L) + sg * c * inv[None, :])
for k in (0.8, 1.25):
    rec(f"ours x {k}", L * k)
(HERE / "work/mde425.json").write_text(json.dumps(out, indent=1))

# --- decomposition: where does the MDE-floor gain come from? called vs uncalled genes, and by expression bin
Q = np.asarray(s["q"], float)
called = (Q[allk, alli] < 0.05)
print(f"\ngate genes {alli.size}: we call {called.mean():.3f} of them; sign acc called {np.mean(np.sign(ours[called]) == np.sign(alll[called])):.3f}"
      f" / uncalled {np.mean(np.sign(ours[~called]) == np.sign(alll[~called])):.3f}")
best = np.nan_to_num(L) + sg * 0.1 * inv[None, :]
for name, Lh in (("ours+0.1/sqrt", best),):
    for lab, msk in (("called", called), ("uncalled", ~called)):
        Lm = np.nan_to_num(L).copy()
        kk, ii = allk[msk], alli[msk]
        Lm[kk, ii] = Lh[kk, ii]
        v = nmae(Lm)
        print(f"  {name} applied only to gate genes we {lab:8s}: nmae {v:.5f}  Δ {v - base:+.5f}")
