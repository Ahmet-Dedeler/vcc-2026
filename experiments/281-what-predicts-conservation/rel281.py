"""281 — POST-HOC: block-bootstrap intervals for the reliability-corrected target `C_rel`, and the
feature-block ablation on the balanced core (so the attribution is read where the measurement-design
confound is gone).

Outputs: work/rel281.json, work/rel281.log
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

import boot281 as B
import fit281 as F

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent / "work"
SEED = F.SEED
T0 = time.time()
LOG: list[str] = []


def P(m: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {m}" if m else ""
    print(line, flush=True)
    LOG.append(line)


BLOCKS = {
    "essentiality": ["dm_mean_effect", "dm_sd_effect", "dm_frac_dependent", "dm_frac_strong",
                     "dm_min_effect", "dm_common_essential", "replogle_essential"],
    "expression": ["ccle_breadth", "ccle_median_tpm", "ccle_sd"],
    "regulon": ["log_collectri_regulon_size", "collectri_indegree", "is_tf_collectri",
                "log_dorothea_regulon_size", "is_tf_dorothea"],
    "complexes": ["humap_n_complexes", "humap_max_complex_size", "reactome_n_complexes",
                  "reactome_max_complex_size"],
    "evolution": ["log_n_paralogues", "ortholog_depth"],
}

df = pd.read_csv(WORK / "feat281.csv", low_memory=False)
held = df[df["set"] == "held"].reset_index(drop=True)
z = pd.read_csv(ROOT / "experiments/279-znf-state-dependence/work/znf279.csv",
                usecols=["gene", "c3_balanced_core"])
held = held.merge(z, on="gene", how="left")
block = np.where(held["gene_group1"].fillna("").to_numpy() != "",
                 held["gene_group1"].fillna("").to_numpy(),
                 "NOFAM:" + held["cytoband"].fillna("").astype(str).to_numpy())
rng = np.random.default_rng(SEED)
blocks = np.array(sorted(set(block)))
order = rng.permutation(len(blocks))
counts = pd.Series(block).value_counts()
cum, take = 0, []
for k in blocks[order]:
    if cum >= F.HELD_FRAC * len(held):
        break
    take.append(k)
    cum += int(counts[k])
te_all = pd.Series(block).isin(set(take)).to_numpy()
X = held[F.FEATS].to_numpy(float)
XO = held[F.RELIABILITY].to_numpy(float)
out: dict = {}

# ---- (1) the reliability-corrected target, block bootstrap -------------------------------------
y = held["y_rel"].to_numpy(float)
ok = np.isfinite(y)
te, tr = te_all & ok, (~te_all) & ok
P(f"C_rel: train {int(tr.sum())}, held-family test {int(te.sum())} "
  f"({len(set(block[te]))} blocks)")
pr, mr = F.fit_ridge(X[tr], y[tr], X[te])
pg, mg, _ = F.fit_gbm(X[tr], y[tr], X[te])
po, mo = F.fit_ridge(XO[tr], y[tr], XO[te])
preds = {"features_ridge": pr, "features_gbm": pg, "reliability_only": po,
         "mean": np.full(int(te.sum()), float(y[tr].mean()))}
res = B.block_boot(y[te], preds, block[te], np.random.default_rng(SEED + 501), n=2000)
for k in [k for k in res if not k.startswith("_")]:
    v = res[k]
    P(f"  {k:20s} rho {v['rho']['point']:+.4f} [{v['rho']['lo']:+.4f}, {v['rho']['hi']:+.4f}]"
      f"   R2 {v['r2']['point']:+.4f} [{v['r2']['lo']:+.4f}, {v['r2']['hi']:+.4f}]"
      f"   AUC10 {v['auc']['point']:.4f} [{v['auc']['lo']:.4f}, {v['auc']['hi']:.4f}]")
v = res["_paired_rho"]["features_ridge-reliability_only"]
P(f"  Δrho features - reliability {v['mean']:+.4f} [{v['lo']:+.4f}, {v['hi']:+.4f}] "
  f"P(Δ<=0) {v['p_le_0']:.3f}")
out["c_rel_block_boot"] = res

# ---- (2) feature-block ablation on the balanced core -------------------------------------------
c3 = held["c3_balanced_core"].to_numpy(float)
ok3 = np.isfinite(c3)
te3, tr3 = te_all & ok3, (~te_all) & ok3
P(f"\nfeature-block ablation on the BALANCED CORE: train {int(tr3.sum())}, test {int(te3.sum())}")
abl = {}
for label, cols in ([(f"only::{b}", c) for b, c in BLOCKS.items()]
                    + [(f"drop::{b}", [f for f in F.FEATS if f not in c]) for b, c in BLOCKS.items()]
                    + [("all", list(F.FEATS))]):
    M = held[cols].to_numpy(float)
    p, m = F.fit_ridge(M[tr3], c3[tr3], M[te3])
    pgb, mgb, _ = F.fit_gbm(M[tr3], c3[tr3], M[te3])
    abl[label] = dict(n_features=len(cols), ridge_rho=F.spearman(c3[te3], p),
                      ridge_r2=F.r2(c3[te3], p), gbm_rho=F.spearman(c3[te3], pgb),
                      gbm_r2=F.r2(c3[te3], pgb))
    a = abl[label]
    P(f"  {label:22s} ({a['n_features']:2d})  ridge rho {a['ridge_rho']:+.4f} R2 "
      f"{a['ridge_r2']:+.4f}   gbm rho {a['gbm_rho']:+.4f} R2 {a['gbm_r2']:+.4f}")
out["core_block_ablation"] = abl

out["log"] = LOG
(WORK / "rel281.json").write_text(json.dumps(out, indent=1, default=float))
P("\nwrote work/rel281.json")
