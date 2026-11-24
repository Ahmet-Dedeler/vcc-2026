"""234 — the COMBINATION key, and the two mandatory nulls.

COMBINATION.  A single fit over every key we hold, trained LEAVE-ONE-FOLD-OUT: the model
sees the other two cancer folds' (key, correctness) pairs and never the fold it is scored
on.  Two variants, because the honest choice is not obvious:
  `pool`  fitted on the reference DE pool only -- the coordinates that are actually scored
  `all`   fitted on a subsample of every scored coordinate -- no knowledge of the pool
and two heads, logistic on rank-transformed features and a plain |t|-weighted linear
discriminant, so the answer does not depend on one estimator's inductive bias.
IN-FOLD is reported alongside OUT-OF-FOLD, because 224 §7 showed a panel-fitted quantity
reads ~0 from power alone and 195 cleared on held perturbations before failing at the line.

NULL 1 -- RANDOM KEYS.  The Phase-1 table is a maximum over 108 (key x region x orientation)
arms.  108 PURE-NOISE keys, swept identically, calibrate exactly how much of the winner is
selection.  A random key selects a random count-matched subset of each row, so this is also
the mandate's sign-shuffle null in its count-matched form.

NULL 2 -- |t|-MASS-MATCHED SHUFFLE.  For each arm, a random subset of the row's POOL
coordinates carrying the same |t| mass the arm flipped.  Its expectation is exactly
`mass * (1 - 2 q_pool)`, so it is the sharpest statement of what a key has to beat.

    uv run python combo234.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common234 import CANCER, QW_NEEDED, RATE_PER_QW, SEED, WORK  # noqa: E402
from screen234 import FRACS, REGIONS, Fold, curve  # noqa: E402

ALL_FOLDS = CANCER + ("kolf2",)
FEATURES = ["absH", "absfin", "abspre", "key189", "theta", "expr", "cpm", "zpred",
            "absH_rank", "whitenkept", "rownorm", "ncov", "margin_maj", "margin_z",
            "coh", "tdonor", "minz", "sumz"]
N_RANDOM_KEYS = 108          # = 18 keys x 3 regions x 2 orientations, the real search size
SUBSAMPLE = 1_500_000


def rank_features(f: Fold) -> np.ndarray:
    """Per-fold rank-transform of every feature to [0, 1], so a model fitted on one fold's
    scale transfers to another.  Shape (T, G, F), float32."""
    T, G = f.truth.shape
    X = np.empty((T * G, len(FEATURES)), dtype=np.float32)
    for c, name in enumerate(FEATURES):
        v = f.keys[name].reshape(-1)
        o = np.argsort(v, kind="stable")
        r = np.empty(v.size, dtype=np.float32)
        r[o] = np.arange(v.size, dtype=np.float32) / v.size
        X[:, c] = r
    return X


def fit_head(X, y, w, kind: str, rng):
    if kind == "logistic":
        from sklearn.linear_model import LogisticRegression
        m = LogisticRegression(max_iter=300, C=1.0, solver="lbfgs")
        m.fit(X, y, sample_weight=w)
        return lambda Z: m.decision_function(Z)
    # |t|-weighted linear discriminant: the weighted mean difference, whitened
    w1 = w * y
    w0 = w * (1 - y)
    m1 = (X * w1[:, None]).sum(0) / max(w1.sum(), 1e-12)
    m0 = (X * w0[:, None]).sum(0) / max(w0.sum(), 1e-12)
    Xc = X - (X * w[:, None]).sum(0) / max(w.sum(), 1e-12)
    S = (Xc * w[:, None]).T @ Xc / max(w.sum(), 1e-12)
    beta = np.linalg.solve(S + 1e-6 * np.eye(S.shape[0]), m1 - m0)
    return lambda Z: Z @ beta


def main() -> None:
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    folds = {h: Fold(h) for h in ALL_FOLDS}
    feats = {h: rank_features(folds[h]) for h in CANCER}
    print(f"  features built ({time.time()-t0:.0f}s)", flush=True)
    audit = {"features": FEATURES, "n_random_keys": N_RANDOM_KEYS}

    # ------------------------------------------------------------------ the combination key
    combo_rows = []
    for train_on in ("pool", "all"):
        for kind in ("logistic", "lda"):
            keys_lofo, keys_infold = {}, {}
            for h in CANCER:
                Xs, ys, ws = [], [], []
                for o in CANCER:
                    if o == h:
                        continue
                    g = folds[o]
                    sel = g.pool.reshape(-1) if train_on == "pool" \
                        else np.ones(g.pool.size, dtype=bool)
                    idx = np.flatnonzero(sel)
                    if idx.size > SUBSAMPLE:
                        idx = rng.choice(idx, SUBSAMPLE, replace=False)
                    Xs.append(feats[o][idx])
                    ys.append(g.correct.reshape(-1)[idx].astype(np.float64))
                    ww = g.w.reshape(-1)[idx] if train_on == "pool" \
                        else np.abs(g.truth).reshape(-1)[idx]
                    ws.append(ww / max(ww.mean(), 1e-30))
                sc = fit_head(np.concatenate(Xs), np.concatenate(ys),
                              np.concatenate(ws), kind, rng)
                keys_lofo[h] = sc(feats[h]).reshape(folds[h].truth.shape).astype(np.float32)
                # in-fold: the same head fitted ON the evaluation fold itself
                g = folds[h]
                sel = g.pool.reshape(-1) if train_on == "pool" \
                    else np.ones(g.pool.size, dtype=bool)
                idx = np.flatnonzero(sel)
                if idx.size > SUBSAMPLE:
                    idx = rng.choice(idx, SUBSAMPLE, replace=False)
                ww = (g.w.reshape(-1)[idx] if train_on == "pool"
                      else np.abs(g.truth).reshape(-1)[idx])
                sc2 = fit_head(feats[h][idx], g.correct.reshape(-1)[idx].astype(np.float64),
                               ww / max(ww.mean(), 1e-30), kind, rng)
                keys_infold[h] = sc2(feats[h]).reshape(g.truth.shape).astype(np.float32)

            for region in REGIONS:
                for tag, kk in (("OUT-of-fold", keys_lofo), ("in-fold", keys_infold)):
                    cur = {h: curve(folds[h], kk[h], region, +1) for h in CANCER}
                    cm = np.mean([cur[h] for h in CANCER], axis=0)
                    i = int(np.argmax(cm))
                    combo_rows.append(dict(
                        head=kind, train_on=train_on, region=region, fit=tag,
                        best_frac=FRACS[i], dqw=float(cm[i]),
                        per_fold={h: float(cur[h][i]) for h in CANCER},
                        overall=float(cm[i]) * RATE_PER_QW["k1.384"],
                        curve={h: [float(v) for v in cur[h]] for h in CANCER}))
            print(f"  combo {train_on}/{kind} done ({time.time()-t0:.0f}s)", flush=True)
    audit["combination"] = combo_rows

    hdr = (f"{'head':<10}{'train':<7}{'region':<11}{'fit':<13}{'f*':>7}{'dqw':>11}"
           f"{'overall':>10}{'k+':>5}")
    print(f"\n{'='*len(hdr)}\nTHE COMBINATION KEY — one fit over all "
          f"{len(FEATURES)} keys, leave-one-fold-out\n{'='*len(hdr)}")
    print(hdr)
    print("-" * len(hdr))
    for r in sorted(combo_rows, key=lambda x: -x["dqw"]):
        npos = sum(1 for v in r["per_fold"].values() if v > 0)
        print(f"{r['head']:<10}{r['train_on']:<7}{r['region']:<11}{r['fit']:<13}"
              f"{r['best_frac']:>7.4f}{r['dqw']:>+11.5f}{r['overall']:>+10.5f}{npos:>4d}/3")

    # ------------------------------------------------------- NULL 1: 108 pure-noise keys
    nulls = []
    for n in range(N_RANDOM_KEYS):
        region = REGIONS[n % len(REGIONS)]
        rk = {h: rng.random(folds[h].truth.shape).astype(np.float32) for h in CANCER}
        cur = {h: curve(folds[h], rk[h], region, +1) for h in CANCER}
        cm = np.mean([cur[h] for h in CANCER], axis=0)
        i = int(np.argmax(cm))
        lofo = []
        for h in CANCER:
            other = np.mean([cur[o] for o in CANCER if o != h], axis=0)
            lofo.append(cur[h][int(np.argmax(other))])
        nulls.append(dict(region=region, best_frac=FRACS[i], dqw_infold=float(cm[i]),
                          dqw_lofo=float(np.mean(lofo))))
    ni = np.array([r["dqw_infold"] for r in nulls])
    nl = np.array([r["dqw_lofo"] for r in nulls])
    audit["null_random_keys"] = dict(
        rows=nulls, infold_max=float(ni.max()), infold_mean=float(ni.mean()),
        infold_p95=float(np.percentile(ni, 95)), infold_sd=float(ni.std(ddof=1)),
        lofo_max=float(nl.max()), lofo_mean=float(nl.mean()),
        lofo_p95=float(np.percentile(nl, 95)), lofo_sd=float(nl.std(ddof=1)))
    print(f"\nNULL 1 — {N_RANDOM_KEYS} PURE-NOISE keys, swept identically "
          "(= the count-matched sign-shuffle, and the multiple-comparison calibration)")
    print(f"  in-fold best-dose d q_w   max {ni.max():+.5f}   p95 {np.percentile(ni,95):+.5f}"
          f"   mean {ni.mean():+.5f}   sd {ni.std(ddof=1):.5f}")
    print(f"  LOFO-dose      d q_w      max {nl.max():+.5f}   p95 {np.percentile(nl,95):+.5f}"
          f"   mean {nl.mean():+.5f}   sd {nl.std(ddof=1):.5f}")

    # ------------------------------------ NULL 2: |t|-mass-matched random pool selections
    screen = json.loads((WORK / "screen234.json").read_text())
    top = sorted(screen["phase1"], key=lambda r: -r["dqw_infold"])[:6]
    mass_rows = []
    for r in top:
        kn, region, sg = r["key"], r["region"], (1 if r["orientation"] == "low_first" else -1)
        frac = r["best_frac"]
        obs, nullv, mass = [], [], []
        for h in CANCER:
            f = folds[h]
            A = f.regions[region]
            k = np.where(A, sg * f.keys[kn].astype(np.float64), np.inf)
            order = np.argsort(k, axis=1, kind="stable")
            m = A.sum(1)
            kk = np.minimum(np.rint(frac * m).astype(int), m)
            sel = np.zeros(f.truth.shape, dtype=bool)
            for i in range(f.truth.shape[0]):
                sel[i, order[i, :kk[i]]] = True
            obs.append(float(np.mean((f.g * sel).sum(1)[f.live])))
            wsel = (f.w * sel).sum(1)
            mass.append(float(np.mean(wsel[f.live] / f.den[f.live])))
            draws = []
            for d in range(200):
                sn = np.zeros(f.truth.shape, dtype=bool)
                for i in np.flatnonzero(f.live):
                    idx = np.flatnonzero(f.pool[i])
                    if idx.size == 0 or wsel[i] <= 0:
                        continue
                    p = rng.permutation(idx)
                    cw = np.cumsum(f.w[i, p])
                    n_take = int(np.searchsorted(cw, wsel[i]) + 1)
                    sn[i, p[:min(n_take, idx.size)]] = True
                draws.append(float(np.mean((f.g * sn).sum(1)[f.live])))
            nullv.append(draws)
        nv = np.array(nullv)                                   # (3 folds, 200 draws)
        cm_null = nv.mean(0)
        mass_rows.append(dict(
            key=kn, region=region, orientation=r["orientation"], frac=frac,
            dqw=float(np.mean(obs)), pool_mass=float(np.mean(mass)),
            null_mean=float(cm_null.mean()), null_sd=float(cm_null.std(ddof=1)),
            null_p95=float(np.percentile(cm_null, 95)),
            z=float((np.mean(obs) - cm_null.mean()) / max(cm_null.std(ddof=1), 1e-12)),
            p_null_beats=float((cm_null >= np.mean(obs)).mean())))
        print(f"  {kn}/{region}@{frac:g}: d q_w {np.mean(obs):+.5f}   pool mass "
              f"{np.mean(mass):.4f}   mass-matched null {cm_null.mean():+.5f} "
              f"+- {cm_null.std(ddof=1):.5f}   z {mass_rows[-1]['z']:+.1f}   "
              f"P(null >= arm) {mass_rows[-1]['p_null_beats']:.3f}", flush=True)
    audit["null_mass_matched"] = mass_rows
    print(f"\nNULL 2 — |t|-mass-matched random POOL selections, 200 draws per fold "
          "(printed above)")

    audit["seconds"] = time.time() - t0
    (WORK / "combo234.json").write_text(json.dumps(audit, indent=1))
    print(f"\n-> {WORK/'combo234.json'}  ({audit['seconds']:.0f}s)   "
          f"gate needs d q_w >= +{QW_NEEDED['k1.384']:.5f}")


if __name__ == "__main__":
    main()
