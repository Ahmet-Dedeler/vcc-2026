#!/usr/bin/env python3
"""700 (coordinator [400b6e]'s ask, 10-14 23:20Z): what predicts the reference's lfc SIGN on our called-but-NOT-significant
genes (721's a0, ~89 % of our calls on h1-test)? A cheap re-check of 651/610's closures on the current call set and the
production emitter (435 seed 0 = lib521 state), h1-test truth.

Set: our prescribed call set (head | CPM block), adjudicable (finite, non-zero truth lfc), NOT truth-significant (q >= 0.05),
own target excluded. Score: a0 = share of genes whose source sign equals sign(truth lfc). Sources (sign per (target, gene)):
  clean4          today's block prior (P_down > 0.5 -> down)
  emitted         our production emitter's delivered lfc sign (what fid actually reads)
  key             our design/bulk key sign (435's cap15 effect, built without H1 donors)
  panel_mean      the panel-mean of our truth-free key (one shared programme per gene)
  pie             PIE lfc_pred sign (xdataset checkpoint, cloud run 20261014T225548Z)
  lib_oracle      sign(beta_g * D_t): beta_g = slope of log1p(CP10k) on log library in the CONTROL cells; D_t = log median
                  library of target t's REAL cells minus the controls' (an ORACLE shift: the day would need a predicted one)
  cpm_rank        no sign by itself: enters only the regression
  logit_cv        cross-fitted logistic regression (5 target folds) on [clean4 P-0.5, key z, panel_mean z, pie lfc z,
                  beta*D (oracle), log control CPM]; reported with and without the oracle term
usage: uv run python experiments/700-pie-transfer/a0sources700.py <pred_h1test.parquet>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402


def main():
    S = L.load_state(0)
    targets, genes = list(S["targets"]), list(S["genes"])
    T, G = len(targets), len(genes)
    # PIE
    t = pq.read_table(sys.argv[1]); pg = json.loads(t.schema.metadata[b"pie"])["genes"]["h1test"]; df = t.to_pandas()
    pos = {g: i for i, g in enumerate(pg)}; have = np.array([g in pos for g in genes]); ix = np.array([pos.get(g, 0) for g in genes])
    plfc = np.zeros((T, G))
    for _, r in df.iterrows():
        plfc[targets.index(r["perturbation"])] = np.where(have, np.asarray(r["lfc_pred"])[ix], 0)
    # CLEAN4
    c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=False)
    cp = {g: i for i, g in enumerate(map(str, c4["genes"]))}
    pdn = np.array([c4["P_down"][cp[g]] if g in cp else np.nan for g in genes])
    # cells: control library slopes and per-target library shifts
    import anndata as ad
    a = ad.read_h5ad(ROOT / "data/harness/h1-test.h5ad")
    vpos = {g: i for i, g in enumerate(a.var_names)}
    cols = np.array([vpos[g] for g in genes])
    X = a.X.tocsr()
    tg = a.obs["target_gene"].to_numpy()
    lib = np.asarray(X.sum(1)).ravel().astype(float)
    ctrl = tg == "non-targeting"
    Xc = X[ctrl][:, cols].toarray().astype(np.float64)
    llc = np.log(lib[ctrl])
    E = np.log1p(1e4 * Xc / lib[ctrl][:, None])
    llc0 = llc - llc.mean()
    beta = (llc0 @ (E - E.mean(0))) / (llc0 @ llc0)            # slope per gene
    cpm = 1e6 * Xc.sum(0) / lib[ctrl].sum()
    D = np.array([np.log(np.median(lib[tg == x])) - np.log(np.median(lib[ctrl])) for x in targets])
    del Xc, E
    own = S["own"]; rl = S["rl"]
    adj = np.isfinite(rl) & (rl != 0) & ~own
    sig = (S["rq"] < L.ALPHA) & np.isfinite(S["rq"]) & ~own
    B = S["C"] & adj & ~sig
    y = np.sign(rl)
    key = S["key"]
    pm = key.mean(0)
    src = {"clean4": np.broadcast_to(np.where(np.isfinite(pdn), np.where(pdn > 0.5, -1.0, 1.0), 0.0), (T, G)),
           "emitted": np.sign(np.nan_to_num(S["ol"])), "key": np.sign(key), "panel_mean": np.broadcast_to(np.sign(pm), (T, G)),
           "pie": np.sign(plfc), "lib_oracle": np.sign(beta[None, :] * D[:, None])}
    res = {"n_B": int(B.sum()), "library_shift_D": {"median": float(np.median(D)), "q10": float(np.quantile(D, .1)),
                                                     "q90": float(np.quantile(D, .9))}, "a0": {}}
    for k, s in src.items():
        m = B & (s != 0)
        res["a0"][k] = dict(a0=round(float((s[m] == y[m]).mean()), 4), coverage=round(float(m.sum() / B.sum()), 4))
    # cross-fitted logistic regression over target folds
    from sklearn.linear_model import LogisticRegression
    def z(v):
        sd = np.nanstd(v)
        return np.nan_to_num((v - np.nanmean(v)) / sd if sd > 0 else v * 0)
    feats_all = {
        "clean4": np.broadcast_to(np.nan_to_num(pdn - 0.5), (T, G)),
        "key": z(key), "panel_mean": np.broadcast_to(z(pm), (T, G)), "pie": z(plfc),
        "lib_oracle": z(beta[None, :] * D[:, None]), "logcpm": np.broadcast_to(z(np.log1p(cpm)), (T, G)),
    }
    fold = np.arange(T) % 5
    for name, use in (("logit_cv_all", list(feats_all)), ("logit_cv_no_oracle", [f for f in feats_all if f != "lib_oracle"]),
                      ("logit_cv_clean4+logcpm", ["clean4", "logcpm"]), ("logit_cv_clean4+lib_oracle", ["clean4", "lib_oracle"])):
        correct = total = 0
        for k in range(5):
            tr = (fold != k)[:, None] & B
            te = (fold == k)[:, None] & B
            Xtr = np.stack([feats_all[f][tr] for f in use], 1); ytr = y[tr] > 0
            Xte = np.stack([feats_all[f][te] for f in use], 1); yte = y[te] > 0
            m = LogisticRegression(max_iter=200).fit(Xtr, ytr)
            correct += int((m.predict(Xte) == yte).sum()); total += int(te.sum())
        res["a0"][name] = dict(a0=round(correct / total, 4), coverage=1.0, features=use)
    out = HERE / "work/a0sources700.json"
    out.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
