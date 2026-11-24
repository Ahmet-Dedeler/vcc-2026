"""192 Phase 2 — the joint (`fid`, `jac`, `reach`) frontier over the call set.

191 proved the production emitter realises ANY prescribed call set at bit-identical per-gene
panel totals (900/900 digests), so the call set is a free design variable.  This script prices
it.

Conventions, all forced by the mandate and by AGENTS.md rule 9 / harness trap 2:
  * `n_real` is a DISTRIBUTION.  Its LEVEL is unknown on the challenge, its SHAPE is h1-test's.
    Every sweep scales each perturbation's pool PROPORTIONALLY (`PowerSurface`), never
    `Surface(n_real_target=N)`, which flattens it and deletes the small-`n_real` stratum.
  * `reach` is reported as a CHECK.  171 proved it is a no-op in the call set; the order is
    `|H|` inside the called set and `|H|` inside its complement (191's `rc_hybrid` law).
  * The emitted `pred_lfc` is a15's own -- the emitter preserves the delivered profile, so the
    SIGNS do not move with the call set.  Only which genes are called moves.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
E185 = os.path.join(REPO, "experiments/185-rank-controlled-emission")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, HERE)
from debench import Surface                             # noqa: E402
from keys192 import gene_stats, h_to_lfc                # noqa: E402

A15 = os.path.join(E185, "work/de185.a15.npz")
ALPHA = 0.05
# board-fitted 1/span (185/181; inside AGENTS.md §4's published spans)
INV = dict(reach=1.12927, fid=3.32482, jac=2.72081, nmae=1.64787)


class PowerSurface(Surface):
    """157's PROPORTIONAL power matching: pool_p = the reference's top round(f*n_real_p)."""

    def __init__(self, path, factor=1.0):
        super().__init__(path)
        z = np.load(path, allow_pickle=True)
        rp = np.where(np.isfinite(z["real_p"]), z["real_p"], np.inf)
        base = (z["real_padj"] < ALPHA) & np.isfinite(z["real_padj"])
        sig = np.zeros(rp.shape, bool)
        for i in range(rp.shape[0]):
            k = int(round(factor * int(base[i].sum())))
            k = max(0, min(k, int(np.isfinite(z["real_p"][i]).sum())))
            if k:
                sig[i, np.argsort(rp[i], kind="stable")[:k]] = True
        self.real_sig = sig & np.isfinite(z["real_p"])
        for p, j in enumerate(self.self_idx):
            if j >= 0:
                self.real_sig[p, j] = False


def topk_call(key, n_per_row):
    out = np.zeros(key.shape, bool)
    for i in range(key.shape[0]):
        n = int(min(max(n_per_row[i], 0), key.shape[1]))
        if n:
            out[i, np.argpartition(-key[i], n - 1)[:n]] = True
    return out


def order_key(call, absH):
    """191's `rc_hybrid` delivery law as a debench ordering key: called genes first (a
    p-threshold block), `|H|` descending inside the block and inside its complement."""
    r = absH / (absH.max(1, keepdims=True) + 1e-12)
    return call.astype(float) * 2.0 + r


def build_keys():
    E = np.load(os.path.join(E169, "work/h1test-effect.npz"), allow_pickle=True)
    pos = {str(g): i for i, g in enumerate(E["genes"])}
    tpos = {str(t): i for i, t in enumerate(E["targets"])}
    s = Surface(A15)
    z15 = np.load(A15, allow_pickle=True)
    gi = np.array([pos[g] for g in s.genes])
    ti = np.array([tpos[t] for t in s.perts])
    cpm = E["cpm"][gi].astype(np.float64)
    Hm = E["H"][np.ix_(ti, gi)].astype(np.float64)
    lfcP = h_to_lfc(Hm, cpm)
    gs = gene_stats(os.path.join(E185, "work/ctrl_cs.npy"),
                    os.path.join(WORK, "genestats.npz"))
    sig = s.real_sig & np.isfinite(s.real_lfc)

    from sklearn.linear_model import LogisticRegression
    X = np.stack([np.broadcast_to(np.log(cpm + 1)[None, :], Hm.shape),
                  np.broadcast_to(np.log(gs["nz"][:] + 1e-3)[None, :], Hm.shape),
                  np.broadcast_to(np.log(gs["cv2"][:] + 1e-3)[None, :], Hm.shape)], -1)
    mu, sd = X.reshape(-1, 3).mean(0), X.reshape(-1, 3).std(0)
    Xz = (X - mu) / sd
    rngf = np.random.default_rng(7)
    fit = np.zeros(Hm.shape)
    coefs = []
    for p in range(len(s.perts)):
        tr = np.setdiff1d(np.arange(len(s.perts)), [p])
        Xt = Xz[tr].reshape(-1, 3)
        yt = sig[tr].reshape(-1)
        idx = rngf.choice(Xt.shape[0], size=min(300_000, Xt.shape[0]), replace=False)
        m = LogisticRegression(max_iter=400).fit(Xt[idx], yt[idx])
        fit[p] = m.decision_function(Xz[p])
        coefs.append(m.coef_[0])
    keys = {
        "|H| (185's key)": np.abs(Hm),
        "a15 p-order (incumbent)": -np.log(np.maximum(z15["pred_p"], 1e-300)),
        "|lfc_pred|*sqrt(cpm)": np.abs(lfcP) * np.sqrt(cpm)[None, :],
        "detectability (cpm)": np.broadcast_to(cpm[None, :], Hm.shape).copy(),
        "FIT (cpm,nz,cv2) LOTO": fit,
        "ORACLE |lfc_real|*sqrt(cpm)": np.abs(np.nan_to_num(s.real_lfc)) * np.sqrt(cpm)[None, :],
        "random": np.random.default_rng(0).random(Hm.shape),
    }
    return s, keys, np.abs(Hm), sig, np.array(coefs).mean(0), (mu, sd)


def main() -> None:
    s, keys, absH, sig, coef, _ = build_keys()
    nreal = sig.sum(1)
    T, G = absH.shape
    print(f"fitted gene-only logistic, mean LOTO coefficients (standardised): "
          f"log cpm {coef[0]:+.3f}, log detection-rate {coef[1]:+.3f}, log cv2 {coef[2]:+.3f}")

    # ---------------------------------------------------------------- n_real_p predictor
    print("\n" + "=" * 100)
    print("SIZING — can we predict n_real_p at build time?  (leave-one-target-out)")
    hn = np.linalg.norm(absH, axis=1)
    feats = np.stack([np.log(hn), np.log(np.abs(absH).sum(1)),
                      np.log((absH > np.quantile(absH, 0.99)).sum(1) + 1)], 1)
    y = np.log(nreal + 1.0)
    pred = np.zeros(T)
    for p in range(T):
        tr = np.setdiff1d(np.arange(T), [p])
        A = np.c_[feats[tr], np.ones(T - 1)]
        beta, *_ = np.linalg.lstsq(A, y[tr], rcond=None)
        pred[p] = np.r_[feats[p], 1.0] @ beta
    ok = nreal > 0
    from scipy.stats import pearsonr, spearmanr
    r = pearsonr(pred[ok], y[ok])
    sp = spearmanr(pred[ok], y[ok])
    resid = y[ok] - pred[ok]
    print(f"  corr(predicted, true log n_real_p): pearson {r[0]:+.4f} (p {r[1]:.3g}), "
          f"spearman {sp.correlation:+.4f}")
    print(f"  residual sd in log space {resid.std():.4f}  -> multiplicative error factor "
          f"exp(sd) = {np.exp(resid.std()):.2f}x")
    print(f"  a constant predictor (median n_real) would have residual sd "
          f"{(y[ok] - np.median(y[ok])).std():.4f}")
    npred_hat = np.exp(pred) - 1.0

    # ---------------------------------------------------------------- the frontier
    Ks = [25, 50, 100, 200, 400, 800, 1600, 3200, 6400, 10779]
    res = {}
    print("\n" + "=" * 100)
    print("THE FRONTIER — global top-K per perturbation.  fid / jac / reach at proportional "
          "power f")
    for f in (0.5, 1.0, 2.0, 3.5):
        sf = PowerSurface(A15, f)
        nr = np.array([int(sf.real_sig[p].sum()) for p in range(T)])
        f0, j0, r0 = sf.fid(), sf.jac(), sf.reach()
        print(f"\n  f = {f}  (median n_real {np.median(nr):.0f}, mean {nr.mean():.0f});  "
              f"a15 REALISED: fid {f0:.5f} jac {j0:.5f} reach {r0:.5f}")
        print(f"  {'key':<28}" + "".join(f"{'K=' + str(K):>22}" for K in Ks))
        for kname, kv in keys.items():
            cells = []
            for K in Ks:
                call = topk_call(kv, np.full(T, K))
                fv, jv = sf.fid(call=call), sf.jac(call=call)
                rv = sf.reach(key=order_key(call, absH))
                res[f"{f}|{kname}|{K}"] = (fv, jv, rv)
                cells.append(f"{fv:>6.4f}/{jv:<6.4f}/{rv:<5.3f}")
            print(f"  {kname:<28}" + "".join(f"{c:>22}" for c in cells))
        # per-perturbation sizing, oracle and predicted, with a margin sweep
        print(f"  {'--- per-perturbation sizing (margin m: n_pred_p = m * n_hat_p) ---':<28}")
        for lab, base in (("ORACLE n_real_p", nr.astype(float)),
                          ("PREDICTED n_hat_p", npred_hat * f)):
            for m in (0.5, 1.0, 1.5, 2.0, 3.0, 5.0):
                for kname in ("FIT (cpm,nz,cv2) LOTO", "a15 p-order (incumbent)"):
                    call = topk_call(keys[kname], np.round(m * base).astype(int))
                    fv, jv = sf.fid(call=call), sf.jac(call=call)
                    rv = sf.reach(key=order_key(call, absH))
                    res[f"{f}|SZ:{lab}|m{m}|{kname}"] = (fv, jv, rv)
                    d = ((fv - f0) * INV["fid"] + (jv - j0) * INV["jac"]) / 6
                    print(f"    {lab:<18} m={m:<4} {kname:<26} fid {fv:.5f} jac {jv:.5f} "
                          f"reach {rv:.5f}   d(fid+jac) overall {d:+.5f}")

    with open(os.path.join(WORK, "frontier192.json"), "w") as fh:
        json.dump(dict(frontier=res, npred_hat=npred_hat.tolist(),
                       nreal=nreal.tolist(), coef=coef.tolist()), fh, indent=1, default=float)
    print(f"\nwrote {WORK}/frontier192.json")


if __name__ == "__main__":
    main()
