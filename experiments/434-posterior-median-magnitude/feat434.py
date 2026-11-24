"""434: build the per-(target, gene) feature table on h1-test's gates for a calibrated sign-confidence model.
Cells = every gene of every scored h1-test target (58 targets with gate >= 10). Columns: truth l (NaN off-gate),
our delivered lfc/q (426's base_s0_bps = the shipped 410b construction), six non-H1 donor deltas, control count/cell.
Writes work/cells434.parquet (gate cells) and work/allcells434.npz (every cell, for the analytic nmae)."""
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
(HERE / "work").mkdir(exist_ok=True)
W323 = ROOT / "experiments/323-current-family-validation/work"
r = np.load(W323 / "reference.npz", allow_pickle=True)
print(sorted(r.files))
genes = r["genes"].astype(str)[r["keep"].astype(bool)]; gi = {g: i for i, g in enumerate(genes)}
s = np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz"); targets = s["targets"].astype(str)
L = np.nan_to_num(s["lfc"]); Q = s["q"]
print("L", L.shape, "Q", Q.shape, "genes", len(genes), "targets", len(targets))
DON = ["k562_gwps", "k562_essential", "rpe1", "nadig_jurkat", "orion_hct116", "orion_hek293t"]
D = np.full((len(DON), len(targets), len(genes)), np.nan, dtype=np.float32)
for di, d in enumerate(DON):
    base = ROOT / "data/donors" / d
    P = [l.strip() for l in open(base / "perts.txt")]; G = [l.strip() for l in open(base / "genes.txt")]
    pi = {p: i for i, p in enumerate(P)}; gcol = {g: i for i, g in enumerate(G)}
    M = np.load(base / "delta.npy", mmap_mode="r")
    cols = np.array([gcol.get(g, -1) for g in genes]); okc = cols >= 0
    for k, t in enumerate(targets):
        if t in pi:
            D[di, k, okc] = np.asarray(M[pi[t]], dtype=np.float32)[cols[okc]]
ref = pl.read_parquet(W323 / "reference_de.parquet")
T = np.full((len(targets), len(genes)), np.nan)
for k, t in enumerate(targets):
    m = ref.filter(pl.col("target") == t)
    g = m["feature"].to_numpy().astype(str); q = m["p_adj"].to_numpy().astype(float); l = m["log2_fold_change"].to_numpy().astype(float)
    ok = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if ok.sum() >= 10:
        T[k, [gi[x] for x in g[ok]]] = l[ok]
np.savez_compressed(HERE / "work/allcells434.npz", L=L.astype(np.float32), Q=Q.astype(np.float32), D=D, T=T.astype(np.float32),
                    targets=targets, genes=genes, donors=np.array(DON))
print("gate cells", np.isfinite(T).sum(), "targets with gate", np.isfinite(T).any(1).sum())
