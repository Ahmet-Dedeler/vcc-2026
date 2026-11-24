"""281 — POST-HOC: does the model reproduce 279's FAMILY AXIS, and is the residual table an artefact?

Two things the residual tail cannot answer on its own:

 1. A group predicted at the top of the ranking can only have a non-positive rank residual, and one
    predicted at the bottom only a non-negative one. The family residual table of core281 must
    therefore be checked against exactly that ceiling/floor before any of it is read as biology.
 2. The scientifically interesting question is whether the build-time features reproduce the axis
    279 found (core machinery conserved -> cell-identity families state-dependent). That is a
    family-level comparison of OBSERVED vs OUT-OF-FOLD PREDICTED medians, not a residual.

Runs on the balanced core (identical three cross-line pairs for every gene).
Outputs: work/fam281.json, work/fam281.log
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import GroupKFold

import fit281 as F

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
SEED = F.SEED

T0 = time.time()
LOG: list[str] = []


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


def main() -> None:
    df = pd.read_csv(WORK / "feat281.csv", low_memory=False)
    held = df[df["set"] == "held"].reset_index(drop=True)
    z = pd.read_csv(ROOT / "experiments/279-znf-state-dependence/work/znf279.csv",
                    usecols=["gene", "c3_balanced_core"])
    held = held.merge(z, on="gene", how="left")
    block = np.where(held["gene_group1"].fillna("").to_numpy() != "",
                     held["gene_group1"].fillna("").to_numpy(),
                     "NOFAM:" + held["cytoband"].fillna("").astype(str).to_numpy())
    c3 = held["c3_balanced_core"].to_numpy(float)
    ok = np.isfinite(c3)
    idx = np.flatnonzero(ok)
    X, y, bl = held[F.FEATS].to_numpy(float)[idx], c3[idx], block[idx]
    sub = held.iloc[idx].copy()

    oof = np.full(len(idx), np.nan)
    for a, b in GroupKFold(n_splits=5).split(X, y, groups=bl):
        oof[b], _ = F.fit_ridge(X[a], y[a], X[b])
    P(f"OOF on the balanced core reproduced: rho {F.spearman(y, oof):+.4f} "
      f"(core281 recorded +0.3398)")
    sub["oof"] = oof
    sub["rank_resid"] = (stats.rankdata(y) - stats.rankdata(oof)) / len(y)

    fam = sub.loc[sub["gene_group1"].fillna("") != ""].groupby("gene_group1")
    tab = fam.agg(n=("gene", "size"), obs=("c3_balanced_core", "median"),
                  pred=("oof", "median"), resid=("rank_resid", "median"))
    tab = tab[tab["n"] >= 20].copy()
    tab["obs_rank"] = stats.rankdata(tab["obs"]) / len(tab)
    tab["pred_rank"] = stats.rankdata(tab["pred"]) / len(tab)
    rho_fam = float(stats.spearmanr(tab["obs"], tab["pred"]).statistic)
    rho_art = float(stats.spearmanr(tab["pred"], tab["resid"]).statistic)
    P(f"\n{len(tab)} HGNC families with >= 20 genes on the core")
    P(f"  Spearman(observed family median, PREDICTED family median) = {rho_fam:+.4f}"
      f"   <- does the model reproduce 279's family axis?")
    P(f"  Spearman(predicted family median, residual family median) = {rho_art:+.4f}"
      f"   <- the ceiling/floor artefact in the residual table")

    # a family-level bootstrap over families
    rng = np.random.default_rng(SEED + 401)
    vals = []
    for _ in range(2000):
        p = rng.integers(0, len(tab), len(tab))
        vals.append(stats.spearmanr(tab["obs"].to_numpy()[p], tab["pred"].to_numpy()[p]).statistic)
    ci = (float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5)))
    P(f"  family bootstrap 95% CI on that axis correlation: [{ci[0]:+.4f}, {ci[1]:+.4f}]")

    P("\n  the two ends of 279 §4(a2)'s axis, observed vs predicted (core, medians):")
    named = ["Spliceosomal C complex", "SSU processome", "Proteasome", "Ribosome biogenesis factors",
             "DEAD-box helicases", "RNA polymerase subunits", "AAA ATPases",
             "Zinc fingers C2H2-type", "KRAB domain containing", "Basic helix-loop-helix proteins",
             "PRD class homeoboxes and pseudogenes", "HOXL subclass homeoboxes",
             "NKL subclass homeoboxes and pseudogenes", "Forkhead boxes",
             "Olfactory receptors, family 2", "Interleukin receptors", "C-type lectin domain containing"]
    rows = []
    for nm in named:
        if nm in tab.index:
            r = tab.loc[nm]
            rows.append(dict(family=nm, n=int(r["n"]), obs=float(r["obs"]),
                             pred=float(r["pred"]), resid=float(r["resid"])))
            P(f"    {nm[:46]:46s} n {int(r['n']):4d}  obs {r['obs']:+.4f}  pred {r['pred']:+.4f}"
              f"  rankres {r['resid']:+.3f}")
    P("\n  full family table sorted by OBSERVED median (top 10 / bottom 10):")
    ts = tab.sort_values("obs")
    for _, r in ts.head(10).reset_index().iterrows():
        P(f"    LOW   {str(r['gene_group1'])[:44]:44s} n {int(r['n']):4d}  obs {r['obs']:+.4f}"
          f"  pred {r['pred']:+.4f}")
    for _, r in ts.tail(10).reset_index().iloc[::-1].iterrows():
        P(f"    HIGH  {str(r['gene_group1'])[:44]:44s} n {int(r['n']):4d}  obs {r['obs']:+.4f}"
          f"  pred {r['pred']:+.4f}")

    out = dict(seed=SEED, n_families=int(len(tab)), rho_obs_pred=rho_fam,
               rho_obs_pred_ci=list(ci), rho_pred_resid_artefact=rho_art,
               named=rows, table=tab.reset_index().to_dict("records"), log=LOG)
    (WORK / "fam281.json").write_text(json.dumps(out, indent=1, default=float))
    P("\nwrote work/fam281.json")


if __name__ == "__main__":
    main()
