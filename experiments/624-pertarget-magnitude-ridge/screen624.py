#!/usr/bin/env python3
"""624 — a learned per-target bulk MAGNITUDE (Testh-style ridge on truth-free features) vs our LODO key + percentile cap.

Council 600 bulk item ([4d0207] -> [7772dc]). Surfaces: 438's October-like panels (efp, efr, real) x the three cancer
folds = 9 folds (522's dumps + saved arms). Stage 1 (this script): an APPROXIMATE screen from 522's saved per-target
arrays. A per-target weight scales a target's row, so x_t(w) ~ x_t(w0) w_t / w0_t and a_t(w) ~ a_t(w0) (w_t / w0_t)^2,
and the aggregate cosine per frame is sum_t x_t / sqrt(sum_t a_t * sum_t d_t) (averaged over frames A/B/C). Stage 2
(exact fold522.score on the best arms) follows only if stage 1 shows a gain.

Learning target (truth-derived, training folds only): the ideal relative weight along our direction on the pctPROD arm,
w*_t = w0_t * (x_t / a_t) / (X / A) (522's u x w), clipped to [0.1, 20], in logs.
Truth-free features per target (from the dump, available on the day):
  log key/median, log LODO terms (min, max, |diff|), cross-donor cosine cos_dd, log donor norms (min, max),
  log Fnorm/median, pnorm, cosPd (min, max), zsrc (518's strength, via 522's z), carried.
Fit: ridge (alpha grid chosen inside the training panels by leave-one-fold-out), leave-one-PANEL-out (train on the other
two panels' 6 folds, predict the held-out panel's 3 folds), so no target or panel leaks.
Arms vs the percentile cap (522's October form), relative aggregate-cosine change in % on each held-out fold:
  ridge    w = clip(exp(pred) / median_carried, 0.5, hi_pct)            (replacement for the key)
  tilt th  w = w_pct * exp(th * zres), zres = ridge prediction residualised on log key, standardised, |z| <= 2
Analysis only, < 2 GB.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W522 = ROOT / "experiments/522-october-bulk-key/work"
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/522-october-bulk-key"))
PANELS = ("efp", "efr", "real")
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t")
P_PROD = 96.64697522769474
EPS = 1e-12


def load(panel, fold):
    d = np.load(W522 / f"dump522-{panel}-{fold}.npz")
    a = np.load(W522 / f"arms522-{panel}-{fold}.npz")
    kv, car = d["kv"].astype(float), d["carried"].astype(bool)
    med = float(np.median(kv[car]))
    V = d["V"].astype(float)          # 2 x n LODO terms
    dn = d["dnorm"].astype(float); cpd = d["cosPd"].astype(float); pn = d["pnorm"].astype(float)
    av = d["avail"].astype(bool)
    with np.errstate(divide="ignore", invalid="ignore"):
        Vn = np.where(av, V / np.median(V[av]), np.nan)
        dnn = np.where(av, dn / np.median(dn[av]), np.nan)
    def slog(x):
        return np.log(np.clip(x, 1e-3, None))
    F = {
        "lkey": slog(np.where(car, kv / med, np.nan)),
        "lV_min": slog(np.nanmin(Vn, 0)), "lV_max": slog(np.nanmax(Vn, 0)),
        "lV_absdiff": np.abs(slog(Vn[0]) - slog(Vn[1])),
        "cos_dd": d["cos_dd"].astype(float),
        "ldn_min": slog(np.nanmin(dnn, 0)), "ldn_max": slog(np.nanmax(dnn, 0)),
        "lF": slog(d["Fnorm"] / np.median(d["Fnorm"])),
        "pn_min": np.nanmin(np.where(av, pn, np.nan), 0), "pn_max": np.nanmax(np.where(av, pn, np.nan), 0),
        "cosPd_min": np.nanmin(np.where(av, cpd, np.nan), 0), "cosPd_max": np.nanmax(np.where(av, cpd, np.nan), 0),
        "zsrc": a["z"].astype(float), "carried": car.astype(float),
    }
    X = np.stack([F[k] for k in F], 1)
    X = np.where(np.isfinite(X), X, np.nan)
    w0 = a["w__pctPROD"].astype(float)
    x0, a0, d0 = a["pctPROD__x"].astype(float), a["pctPROD__a"].astype(float), a["pctPROD__d"].astype(float)
    xs, as_ = x0.sum(0), a0.sum(0)
    Xs, As = xs.sum(), as_.sum()
    u = (xs / np.maximum(as_, EPS)) / (Xs / As)
    wstar = np.clip(w0 * u, 0.1, 20.0)
    return dict(panel=panel, fold=fold, X=X, names=list(F), y=np.log(wstar), car=car, kv=kv, med=med, w0=w0,
                x0=x0, a0=a0, d0=d0, hi=float(np.percentile(kv[car] / med, P_PROD)))


def agg_cos(S, w):
    r = np.asarray(w, float) / np.maximum(S["w0"], EPS)
    out = []
    for f in range(S["x0"].shape[0]):
        x = (S["x0"][f] * r).sum(); a = (S["a0"][f] * r * r).sum(); dd = S["d0"][f].sum() if S["d0"].ndim > 1 else S["d0"].sum()
        out.append(x / np.sqrt(a * dd))
    return float(np.mean(out))


def ridge_fit(X, y, alpha):
    mu = np.nanmean(X, 0); sd = np.nanstd(X, 0); sd[sd == 0] = 1
    Z = np.nan_to_num((X - mu) / sd)
    A = Z.T @ Z + alpha * np.eye(Z.shape[1])
    beta = np.linalg.solve(A, Z.T @ (y - y.mean()))
    return dict(mu=mu, sd=sd, beta=beta, b0=float(y.mean()))


def ridge_pred(m, X):
    return m["b0"] + np.nan_to_num((X - m["mu"]) / m["sd"]) @ m["beta"]


def main() -> None:
    S = {(p, f): load(p, f) for p in PANELS for f in FOLDS}
    names = S[(PANELS[0], FOLDS[0])]["names"]
    alphas = (0.1, 1.0, 10.0, 100.0, 1000.0)
    out = {"features": names, "rows": []}
    for hp in PANELS:
        train = [k for k in S if k[0] != hp]
        # alpha by leave-one-fold-out inside the training panels (Spearman of pred vs y on carried)
        def score_alpha(al):
            rs = []
            for k in train:
                tr = [j for j in train if j != k]
                X = np.concatenate([S[j]["X"][S[j]["car"]] for j in tr]); y = np.concatenate([S[j]["y"][S[j]["car"]] for j in tr])
                m = ridge_fit(X, y, al)
                p = ridge_pred(m, S[k]["X"][S[k]["car"]])
                from scipy.stats import spearmanr
                rs.append(spearmanr(p, S[k]["y"][S[k]["car"]]).correlation)
            return float(np.mean(rs))
        al_scores = {al: score_alpha(al) for al in alphas}
        al = max(al_scores, key=al_scores.get)
        X = np.concatenate([S[j]["X"][S[j]["car"]] for j in train]); y = np.concatenate([S[j]["y"][S[j]["car"]] for j in train])
        m = ridge_fit(X, y, al)
        for f in FOLDS:
            T = S[(hp, f)]
            car = T["car"]
            pred = ridge_pred(m, T["X"])
            from scipy.stats import spearmanr
            rho = spearmanr(pred[car], T["y"][car]).correlation
            rho_key = spearmanr(T["X"][car, 0], T["y"][car]).correlation
            base = agg_cos(T, T["w0"])
            res = {"panel": hp, "fold": f, "alpha": al, "alpha_cv": al_scores, "spearman_pred_vs_ideal": float(rho),
                   "spearman_key_vs_ideal": float(rho_key), "base_cos": base}
            # ridge replacement
            wr = np.ones_like(T["w0"]); e = np.exp(pred[car]); wr[car] = np.clip(e / np.median(e), 0.5, T["hi"])
            res["ridge_vs_pct_pct"] = 100 * (agg_cos(T, wr) / base - 1)
            # tilt on the percentile cap: residual of pred on log key
            lk = T["X"][:, 0]
            ok = car & np.isfinite(lk)
            A_ = np.vstack([np.ones(ok.sum()), lk[ok]]).T
            coef, *_ = np.linalg.lstsq(A_, pred[ok], rcond=None)
            res_ = np.zeros_like(pred); res_[ok] = pred[ok] - A_ @ coef
            z = np.zeros_like(pred); z[ok] = np.clip(res_[ok] / (res_[ok].std() + EPS), -2, 2)
            for th in (0.1, 0.25, 0.5):
                wt = T["w0"] * np.exp(th * z)
                res[f"tilt{th}_vs_pct_pct"] = 100 * (agg_cos(T, wt) / base - 1)
            # 522's strength rider (pct x exp(0.25 z_src) on targets above the floor) and the ridge tilt ON TOP of it
            a522 = np.load(W522 / f"arms522-{hp}-{f}.npz")
            zsrc = a522["z"].astype(float)
            above = car & (T["kv"] / T["med"] > 0.5)
            wS = T["w0"] * np.exp(0.25 * np.where(above, zsrc, 0.0))
            res["S522_0.25_vs_pct_pct"] = 100 * (agg_cos(T, wS) / base - 1)
            for th in (0.1, 0.25):
                res[f"S522+tilt{th}_vs_pct_pct"] = 100 * (agg_cos(T, wS * np.exp(th * z)) / base - 1)
            # oracle bounds for scale: ideal weights (truth) on this fold
            wo = np.where(car, np.exp(T["y"]), 1.0)
            res["oracle_vs_pct_pct"] = 100 * (agg_cos(T, wo) / base - 1)
            out["rows"].append(res)
            print(f"{hp:4s} {f:14s} α {al:g}  ρ(pred, ideal) {rho:+.3f} (key {rho_key:+.3f}) | ridge {res['ridge_vs_pct_pct']:+.2f}%  "
                  f"tilt .1 {res['tilt0.1_vs_pct_pct']:+.2f}%  .25 {res['tilt0.25_vs_pct_pct']:+.2f}%  .5 {res['tilt0.5_vs_pct_pct']:+.2f}%  "
                  f"| oracle {res['oracle_vs_pct_pct']:+.1f}%", flush=True)
        out[f"coef_{hp}"] = dict(zip(names, map(float, m["beta"])))
    for k in ("ridge_vs_pct_pct", "tilt0.1_vs_pct_pct", "tilt0.25_vs_pct_pct", "tilt0.5_vs_pct_pct", "S522_0.25_vs_pct_pct",
              "S522+tilt0.1_vs_pct_pct", "S522+tilt0.25_vs_pct_pct"):
        v = [r[k] for r in out["rows"]]
        print(f"{k:22s} mean {np.mean(v):+.2f}%  folds > 0: {sum(x > 0 for x in v)}/9")
    (WORK / "screen624.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
