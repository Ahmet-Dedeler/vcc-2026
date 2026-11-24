"""284 — the preregistered tables: S1 (primary), S2, 281's own split and the random reference, on
`C_held`, `C_rel` and the balanced core; component block bootstrap (PREREG §5).

  PYTHONDONTWRITEBYTECODE=1 uv run python fit284.py

Outputs: work/fit284.json, work/fit284.log, work/preds284.npz
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor

import lib284 as L

T0 = time.time()
LOG: list[str] = []


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


def arms(held: pd.DataFrame, y: np.ndarray, te_all: np.ndarray) -> tuple[dict, dict]:
    """281's arm set, every hyperparameter chosen inside TRAIN (fit281.run + boot281 preds)."""
    ok = np.isfinite(y)
    te, tr = te_all & ok, (~te_all) & ok
    X = held[L.FEATS].to_numpy(float)
    XR = held[L.FEATS + L.RELIABILITY].to_numpy(float)
    XO = held[L.RELIABILITY].to_numpy(float)
    j = L.FEATS.index("replogle_essential")
    meta: dict = {"n_train": int(tr.sum()), "n_test": int(te.sum())}
    preds: dict = {}
    preds["features_ridge"], meta["ridge"] = L.fit_ridge(X[tr], y[tr], X[te])
    preds["features_gbm"], meta["gbm"], _ = L.fit_gbm(X[tr], y[tr], X[te])
    preds["reliability_only"], meta["rel_ridge"] = L.fit_ridge(XO[tr], y[tr], XO[te])
    preds["reliability_only_gbm"], meta["rel_gbm"], _ = L.fit_gbm(XO[tr], y[tr], XO[te])
    preds["features_plus_reliability"], meta["fr_ridge"] = L.fit_ridge(XR[tr], y[tr], XR[te])
    preds["features_plus_reliability_gbm"], meta["fr_gbm"], _ = L.fit_gbm(XR[tr], y[tr], XR[te])
    preds["best_single_calibrated"], meta["best_single"] = L.fit_ridge(
        X[:, [j]][tr], y[tr], X[:, [j]][te])
    preds["mean"] = np.full(int(te.sum()), float(y[tr].mean()))
    return preds, dict(meta=meta, te=te, tr=tr)


def show(res: dict, label: str) -> None:
    P(f"  --- {label} ---")
    for k in [k for k in res if not k.startswith("_")]:
        v = res[k]
        P(f"    {k:30s} rho {v['rho']['point']:+.4f} [{v['rho']['lo']:+.4f}, {v['rho']['hi']:+.4f}]"
          f" (sd {v['rho']['sd']:.4f})  R2 {v['r2']['point']:+.4f} [{v['r2']['lo']:+.4f}, "
          f"{v['r2']['hi']:+.4f}]  AUC10 {v['auc']['point']:.4f} [{v['auc']['lo']:.4f}, "
          f"{v['auc']['hi']:.4f}]")
    for k in ("features_ridge-reliability_only", "features_gbm-reliability_only_gbm",
              "features_plus_reliability-reliability_only",
              "features_plus_reliability_gbm-reliability_only_gbm",
              "features_ridge-best_single_calibrated"):
        if k in res["_paired_rho"]:
            v = res["_paired_rho"][k]
            P(f"    Δrho {k:56s} {v['mean']:+.4f} [{v['lo']:+.4f}, {v['hi']:+.4f}] "
              f" P(Δ<=0) {v['p_le_0']:.3f}")


def main() -> None:
    h = hashlib.sha256((L.HERE / "PREREGISTRATION.md").read_bytes()).hexdigest()
    frozen = (L.WORK / "PREREGISTRATION.sha256").read_text().split()[0]
    if h != frozen:
        sys.exit("preregistration digest mismatch — aborting")
    P(f"prereg {h} MATCH")

    held = L.load_held()
    sp = np.load(L.WORK / "splits284.npz", allow_pickle=True)
    te281, teS1, teS2, teRND = sp["te281"], sp["teS1"], sp["teS2"], sp["teRND"]
    comp, block = sp["comp"], np.asarray(sp["block"]).astype(str)
    out: dict = {"prereg_sha256": h, "results": {}}
    store: dict = {}

    targets = [("C_held", held["y_raw"].to_numpy(float)),
               ("C_rel", held["y_rel"].to_numpy(float)),
               ("core", held["c3_balanced_core"].to_numpy(float))]
    splits = [("S1_all_membership", teS1), ("S2_component", teS2),
              ("split281_first_group", te281), ("random_20260907", teRND)]

    for tname, y in targets:
        for sname, te_all in splits:
            key = f"{tname}::{sname}"
            P(f"\n=== {key} ===")
            preds, info = arms(held, y, te_all)
            te = info["te"]
            P(f"  train {info['meta']['n_train']}, test {info['meta']['n_test']}; "
              f"ridge alpha {info['meta']['ridge']['alpha']:.4g}; "
              f"gbm {info['meta']['gbm']['params']['learning_rate']}/"
              f"{info['meta']['gbm']['params']['max_iter']}; test spans "
              f"{len(set(comp[te].tolist()))} components / {len(set(block[te].tolist()))} "
              f"first-family blocks")
            res_c = L.block_boot(y[te], preds, comp[te], np.random.default_rng(L.SEED + 11))
            show(res_c, "component bootstrap (frozen, PREREG §5)")
            entry = dict(meta=info["meta"], n_components=len(set(comp[te].tolist())),
                         n_first_family_blocks=len(set(block[te].tolist())),
                         boot_component=res_c)
            # 281's first-family block bootstrap, for the comparison the mandate asks for
            if tname == "C_held":
                res_b = L.block_boot(y[te], preds, block[te], np.random.default_rng(L.SEED + 101))
                show(res_b, "281's first-family block bootstrap (comparison)")
                entry["boot_first_family"] = res_b
            out["results"][key] = entry
            if tname == "C_held" and sname in ("S1_all_membership", "S2_component"):
                store[f"{sname}::idx"] = np.flatnonzero(te)
                for k, v in preds.items():
                    store[f"{sname}::{k}"] = np.asarray(v)

    # ---- KA-4: shuffled training target on S1, 200 independent shuffles (281's corrected form) --
    P("\n=== KA-4: shuffled training target on S1, 200 shuffles (ridge, alpha fixed at the "
      "real-data choice) ===")
    y = held["y_raw"].to_numpy(float)
    te, tr = teS1, ~teS1
    X = held[L.FEATS].to_numpy(float)
    alpha = out["results"]["C_held::S1_all_membership"]["meta"]["ridge"]["alpha"]
    st = L.Standardiser().fit(X[tr])
    from sklearn.linear_model import Ridge
    Xtr, Xte = st.transform(X[tr]), st.transform(X[te])
    rng = np.random.default_rng(L.SEED + 77)
    null = []
    for _ in range(200):
        ys = y[tr].copy()
        rng.shuffle(ys)
        null.append(L.spearman(y[te], Ridge(alpha=alpha).fit(Xtr, ys).predict(Xte)))
    null = np.asarray(null)
    obs = out["results"]["C_held::S1_all_membership"]["boot_component"]["features_ridge"]["rho"][
        "point"]
    ka4 = dict(null_mean=float(null.mean()), null_sd=float(null.std()), observed=float(obs),
               n_reaching=int((null >= obs).sum()), n_draws=200,
               z=float((obs - null.mean()) / null.std()),
               pass_=bool((null >= obs).sum() == 0))
    P(f"  null {ka4['null_mean']:+.4f} ± {ka4['null_sd']:.4f}; observed {obs:+.4f}; "
      f"{ka4['n_reaching']}/200 draws reach it; z {ka4['z']:+.1f} -> "
      f"{'PASS' if ka4['pass_'] else 'FAIL'}")
    out["ka4"] = ka4

    np.savez_compressed(L.WORK / "preds284.npz", **store)
    out["log"] = LOG
    (L.WORK / "fit284.json").write_text(json.dumps(out, indent=1, default=float))
    (L.WORK / "fit284.log").write_text("\n".join(LOG))
    P("wrote work/fit284.json, work/preds284.npz")


if __name__ == "__main__":
    main()
