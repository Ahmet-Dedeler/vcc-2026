"""425: what do the CRISPRi cis-neighbours buy nmae on H1, given our shipped (395) delivered statistics?

Base = 410/work/h1/stats-base395.npz (our delivered per-gene lfc/q/p on h1-test), exactly 412's nmae.
Rows: our delivered lfc at the cis pairs, the oracle upper bound (cis pairs set to the reference lfc), and oracle-free
rules that set the pair's lfc to a constant c (log2) when the neighbour's TSS is within D of the target's.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
W323 = ROOT / "experiments/323-current-family-validation/work"
r = np.load(W323 / "reference.npz", allow_pickle=True)
s = np.load(ROOT / "experiments/410-pertscale-revisit/work/h1/stats-base395.npz", allow_pickle=True)
genes = r["genes"].astype(str)[r["keep"].astype(bool)]
targets = s["targets"].astype(str)
ref = pl.read_parquet(W323 / "reference_de.parquet")
gi = {g: i for i, g in enumerate(genes)}
L = np.asarray(s["lfc"], float)
tss = pd.read_csv(HERE / "work/tss_v44.tsv", sep="\t").set_index("gene_name")

gates = []
for k, t in enumerate(targets):
    m = ref.filter(pl.col("target") == t)
    g = m["feature"].to_numpy().astype(str); q = m["p_adj"].to_numpy().astype(float)
    l = m["log2_fold_change"].to_numpy().astype(float)
    ok = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if ok.sum() >= 10:
        gates.append((k, np.array([gi[x] for x in g[ok]]), l[ok]))
REF = np.full(L.shape, np.nan)
for k, t in enumerate(targets):
    m = ref.filter(pl.col("target") == t)
    idx = np.array([gi.get(x, -1) for x in m["feature"].to_numpy().astype(str)])
    ok = idx >= 0
    REF[k, idx[ok]] = m["log2_fold_change"].to_numpy()[ok]


def nmae(Lhat):
    v = []
    for k, idx, l in gates:
        p = np.nan_to_num(Lhat[k, idx], nan=0.0)
        v.append(np.abs(p - l).sum() / np.abs(l).sum())
    return float(np.mean(v))


# cis pairs on the delivered gene axis
gt = tss.reindex(genes)
pairs = []
for k, t in enumerate(targets):
    if t not in tss.index:
        continue
    tt = tss.loc[t]
    same = (gt["chrom"].to_numpy() == tt["chrom"])
    d = np.where(same, np.abs(gt["tss"].to_numpy(dtype=float) - float(tt["tss"])), np.inf)
    for j in np.flatnonzero(d <= 5000):
        if genes[j] == t:
            continue
        pairs.append((k, j, d[j]))
print(f"targets {len(targets)}  gated {len(gates)}  cis pairs <=5kb on axis {len(pairs)}")
for k, j, d in sorted(pairs, key=lambda x: x[2]):
    print(f"  {targets[k]:10s} {genes[j]:10s} {int(d):6d}  ours {L[k, j]:+.3f}  ref {REF[k, j]:+.3f}")

base = nmae(L)
out = {"base": base, "rows": {}}
print(f"\nbase nmae {base:.5f} (412/410 base395: 0.95210)")


def rec(name, Lh):
    v = nmae(Lh); out["rows"][name] = v; print(f"  {name:44s} nmae {v:.5f}  Δ {v - base:+.5f}")


for D in (1000, 5000):
    sel = [(k, j) for k, j, d in pairs if d <= D]
    Lo = L.copy()
    for k, j in sel:
        if np.isfinite(REF[k, j]):
            Lo[k, j] = REF[k, j]
    rec(f"ORACLE cis<= {D}", Lo)
    for c in (-0.5, -1.0, -1.5, -2.0, -3.0):
        Lc = L.copy()
        for k, j in sel:
            Lc[k, j] = min(L[k, j], c) if np.isfinite(L[k, j]) else c
        rec(f"rule cis<= {D}: lfc = min(ours, {c})", Lc)
(HERE / "work/h1nmae425.json").write_text(json.dumps(out, indent=1))
