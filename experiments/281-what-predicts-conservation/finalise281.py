"""281 — add the balanced-core target and its out-of-fold prediction to work/pred281.csv."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

import fit281 as F

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent / "work"

df = pd.read_csv(WORK / "feat281.csv", low_memory=False)
held = df[df["set"] == "held"].reset_index(drop=True)
z = pd.read_csv(ROOT / "experiments/279-znf-state-dependence/work/znf279.csv",
                usecols=["gene", "c3_balanced_core"])
held = held.merge(z, on="gene", how="left")
block = np.where(held["gene_group1"].fillna("").to_numpy() != "",
                 held["gene_group1"].fillna("").to_numpy(),
                 "NOFAM:" + held["cytoband"].fillna("").astype(str).to_numpy())
c3 = held["c3_balanced_core"].to_numpy(float)
idx = np.flatnonzero(np.isfinite(c3))
X, y, bl = held[F.FEATS].to_numpy(float)[idx], c3[idx], block[idx]
oof = np.full(len(idx), np.nan)
for a, b in GroupKFold(n_splits=5).split(X, y, groups=bl):
    oof[b], _ = F.fit_ridge(X[a], y[a], X[b])
print(f"core OOF rho {F.spearman(y, oof):+.4f} on n {len(idx)}")

col_obs = np.full(len(held), np.nan)
col_obs[idx] = y
col_pred = np.full(len(held), np.nan)
col_pred[idx] = oof
core = pd.DataFrame(dict(gene=held["gene"], observed_core=col_obs,
                         pred_oof_core_ridge=col_pred, n_pairs=held["n_pairs"]))

pred = pd.read_csv(WORK / "pred281.csv")
pred = pred.merge(core, on="gene", how="left")
pred.to_csv(WORK / "pred281.csv", index=False)
print(f"work/pred281.csv now {pred.shape[0]} rows x {pred.shape[1]} cols: "
      f"{', '.join(pred.columns)}")
