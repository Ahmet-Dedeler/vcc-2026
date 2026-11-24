"""521 steps 3-6 — per-target call-set policies driven by the strength key, priced exactly from table521 on all four DE
members, at each power regime, with target-bootstrap CIs and the top-3 share.

Policies (every one a vector K_p over the 85 H1 targets; 'base' = today's call set, 2,600 prescribed):
  ORACLE     K_p = best in the grid for target p (per-target oracle, the ceiling at this content)
  NEG-ORACLE n_real_p < 20 -> K (which targets are negligible, known)
  RANK(q,K)  the bottom q % of targets by strength -> K; the rest keep base          (no fitted parameter)
  TWO(q1,K1,q2,K2) bottom q1 % -> K1, the next q2 % -> K2
  BAYES      the Bayes action: LOTO Gaussian predictive distribution of log1p n_real (linear in the strength rank, sigma
             from the other targets' residuals) x a LOTO kernel regression of each action's per-target utility on
             log1p n_real; K_p = argmax expected utility (base has 0). Everything that sets target p's action is
             fitted without target p.
Predictors: the strength ranks (S4 thr4 primary = the H1-available part of 518's meta key; S3 = nreal521's; K562 alone;
S4 thr-avg), and DEGRADED versions (Gaussian noise on the ranks, 200 draws) that bring the rank correlation with
log1p n_real down to the cancer folds' level (r 0.40 / 0.25; pred521) — the board's predictor is closer to those.
Units: board overall (192's INV; nmae lower is better); pds/mse do not move in this family.

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/policy521.py [seed]
"""
from __future__ import annotations

import json
import sys

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402
from scipy.stats import norm  # noqa: E402

import lib521 as L  # noqa: E402

MEM = ("fid", "jac", "reach", "nmae")
SGN = np.array([L.INV["fid"], L.INV["jac"], L.INV["reach"], -L.INV["nmae"]]) / 6.0
QS = (5, 10, 15, 20, 25, 30, 40, 50)
KS = (0, 10, 25, 50, 100, 200, 400, 800, 1200, 1600, 2000, 2200)
NB = 4000


class Table:
    def __init__(self, seed=0):
        z = np.load(L.WORK / f"table521_s{seed}.npz")
        self.rows = {str(k): r for k, r in zip(z["keys"], z["rows"])}
        self.nreal = {str(r): n for r, n in zip(z["regimes"], z["nreal"])}
        self.kgrid = [int(k) for k in z["kgrid"]]

    def block(self, reg, kind, var, K):
        return self.rows[f"{reg}|base"] if K is None else self.rows[f"{reg}|{kind}|{var}|{K}"]

    def assemble(self, reg, kind, var, Kvec):
        """rows (4 x 85) of a policy; Kvec entries None = base."""
        out = self.rows[f"{reg}|base"].copy()
        for i, K in enumerate(Kvec):
            if K is not None:
                out[:, i] = self.rows[f"{reg}|{kind}|{var}|{K}"][:, i]
        return out


def means(R, idx=None):
    if idx is not None:
        R = R[:, idx]
    return np.array([np.nanmean(R[m]) for m in range(4)])


def value(R, R0, idx=None):
    return float(SGN @ (means(R, idx) - means(R0, idx)))


def boot(R, R0, nb=NB, seed=521):
    rng = np.random.default_rng(seed)
    n = R.shape[1]
    v = np.empty(nb)
    for b in range(nb):
        idx = rng.integers(0, n, n)
        v[b] = value(R, R0, idx)
    return float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975)), float(np.mean(v <= 0))


def per_target_contrib(R, R0):
    """each target's contribution to the overall delta (cohort sizes of the base)."""
    c = np.zeros(R.shape[1])
    for m in range(4):
        ok0 = np.isfinite(R0[m])
        N = ok0.sum()
        d = np.where(np.isfinite(R[m]) & ok0, R[m] - R0[m], 0.0)
        if m == 0:   # a target leaving fid's cohort (n_pred = n_real = 0) moves the mean by (mean - row)/N
            gone = ok0 & ~np.isfinite(R[m])
            d = np.where(gone, np.nanmean(R0[m]) - R0[m], d)
        c += SGN[m] * d / N
    return c


def top3_share(c):
    tot = c.sum()
    if tot <= 0:
        return np.nan
    return float(np.sort(c)[::-1][:3].sum() / tot)


def summarize(T, reg, kind, var, Kvec, label, out, verbose=True):
    R0 = T.rows[f"{reg}|base"]
    R = T.assemble(reg, kind, var, Kvec)
    v = value(R, R0)
    lo, hi, p0 = boot(R, R0)
    dm = means(R) - means(R0)
    c = per_target_contrib(R, R0)
    rec = dict(value=v, ci=[lo, hi], p_le0=p0, d=dict(zip(MEM, dm.tolist())), top3=top3_share(c),
               n_shrunk=int(sum(K is not None for K in Kvec)))
    out[f"{reg}|{kind}|{var}|{label}"] = rec
    if verbose:
        print(f"    {label:30s} {v:+.5f} [{lo:+.5f}, {hi:+.5f}] P<=0 {p0:.3f}  fid {dm[0]:+.5f} jac {dm[1]:+.5f} "
              f"reach {dm[2]:+.5f} nmae {dm[3]:+.5f}  shrunk {rec['n_shrunk']:2d}  top3 {rec['top3']:.2f}")
    return rec


def rank_policy(z, q, K):
    thr = np.quantile(z, q / 100.0)
    n = z.size
    k = int(round(q / 100.0 * n))
    o = np.argsort(z, kind="stable")[:k]
    Kv = [None] * n
    for i in o:
        Kv[i] = K
    return Kv


def two_level(z, q1, K1, q2, K2):
    n = z.size
    o = np.argsort(z, kind="stable")
    k1, k2 = int(round(q1 / 100 * n)), int(round((q1 + q2) / 100 * n))
    Kv = [None] * n
    for i in o[:k1]:
        Kv[i] = K1
    for i in o[k1:k2]:
        Kv[i] = K2
    return Kv


def bayes_policy(T, reg, kind, var, z, h=0.5, margin=0.0):
    """LOTO Bayes action over {base} U KS."""
    n = T.nreal[reg].astype(float)
    y = np.log1p(n)
    R0 = T.rows[f"{reg}|base"]
    U = {}
    for K in KS:
        R = T.assemble(reg, kind, var, [K] * 85)
        U[K] = per_target_contrib(R, R0)          # utility of shrinking target q to K (alone), per target
    grid = np.linspace(0, max(y.max(), 1.0) + 2.0, 161)
    Kv = [None] * 85
    eu_rec = []
    for p in range(85):
        m = np.arange(85) != p
        b, a = np.polyfit(z[m], y[m], 1)
        sd = float(np.std(y[m] - (a + b * z[m]), ddof=2))
        mu = a + b * z[p]
        w = norm.pdf(grid, mu, sd)
        w /= w.sum()
        best, bestK = margin, None
        for K in KS:
            kw = norm.pdf((grid[:, None] - y[m][None, :]) / h)
            uhat = (kw * U[K][m][None, :]).sum(1) / np.maximum(kw.sum(1), 1e-300)
            eu = float((w * uhat).sum())
            if eu > best:
                best, bestK = eu, K
        Kv[p] = bestK
        eu_rec.append(best)
    return Kv


def degrade(z, y, r_target, rng):
    """rank-noise so that corr(rank(z + e), y) ~ r_target (noise scale solved on this draw's noise shape)."""
    zz = (z - z.mean()) / z.std()
    e = rng.standard_normal(z.size)
    lo, hi = 0.0, 20.0
    for _ in range(40):
        s = 0.5 * (lo + hi)
        zn = zz + s * e
        r = np.corrcoef(np.argsort(np.argsort(zn)), y)[0, 1]
        if r > r_target:
            lo = s
        else:
            hi = s
    zn = zz + 0.5 * (lo + hi) * e
    return (np.argsort(np.argsort(zn)) / (z.size - 1) - 0.5)


def main(seed=0):
    T = Table(seed)
    P = np.load(L.WORK / "pred521_h1.npz")
    feats = {"S4": P["S4__S3_p_vip_thr4"], "S3": P["S3__k562phct116phek293t_nreal521_thr4"],
             "K562": P["k562_t4"], "S4avg": P["S4_thr-avg"]}
    # K562 alone lacks 3 targets: give them the S4 rank (the build would do the same)
    k = feats["K562"].copy()
    kr = np.full(85, np.nan)
    ok = np.isfinite(k)
    kr[ok] = np.argsort(np.argsort(k[ok])) / (ok.sum() - 1) - 0.5
    kr[~ok] = feats["S4"][~ok]
    feats["K562"] = kr
    out = {}
    for reg in ("1.0", "0.5", "0.35", "cmp"):
        n = T.nreal[reg]
        print(f"\n=== regime {reg}: n_real median {np.median(n):.0f}, <20 {int((n < 20).sum())}, zero {int((n == 0).sum())}")
        for var in ("b", "a"):
            kind = "headcpm"
            print(f"  -- content {kind}, variant {var}")
            # oracles
            R0 = T.rows[f"{reg}|base"]
            best = []
            for i in range(85):
                vals = [0.0] + [per_target_contrib(T.assemble(reg, kind, var, [K if j == i else None for j in range(85)]), R0)[i]
                                for K in KS]
                j = int(np.argmax(vals))
                best.append(None if j == 0 else KS[j - 1])
            summarize(T, reg, kind, var, best, "ORACLE per-target best K", out)
            for K in (10, 25, 50, 100, 400):
                summarize(T, reg, kind, var, [K if n[i] < 20 else None for i in range(85)], f"NEG-ORACLE n<20 -> {K}", out)
            summarize(T, reg, kind, var, [0 if n[i] == 0 else None for i in range(85)], "NULL-ORACLE n=0 -> 0", out)
            for fname in ("S4", "S3", "K562", "S4avg"):
                z = feats[fname]
                verbose = fname == "S4"
                if verbose:
                    print(f"   {fname}: RANK(q, K)")
                grid = {}
                R0g = T.rows[f"{reg}|base"]
                for q in QS:
                    for K in KS:
                        grid[(q, K)] = value(T.assemble(reg, kind, var, rank_policy(z, q, K)), R0g)
                out[f"{reg}|{kind}|{var}|{fname} RANK grid"] = {f"q{q}_K{K}": v for (q, K), v in grid.items()}
                bestqk = max(grid, key=grid.get)
                if verbose:
                    hdr = "      q\\K " + "".join(f"{K:>9d}" for K in KS)
                    print(hdr)
                    for q in QS:
                        print(f"      {q:3d}  " + "".join(f"{grid[(q, K)]:+9.5f}" for K in KS))
                rec = summarize(T, reg, kind, var, rank_policy(z, *bestqk), f"{fname} RANK best q{bestqk[0]} K{bestqk[1]}", out)
                for q, K in ((10, 25), (20, 100), (20, 400), (30, 400)):
                    summarize(T, reg, kind, var, rank_policy(z, q, K), f"{fname} RANK q{q} K{K}", out, verbose=verbose)
                for (q1, K1, q2, K2) in ((10, 25, 20, 400), (10, 50, 30, 800), (5, 0, 20, 100)):
                    summarize(T, reg, kind, var, two_level(z, q1, K1, q2, K2), f"{fname} TWO {q1}/{K1}+{q2}/{K2}", out,
                              verbose=verbose)
                summarize(T, reg, kind, var, bayes_policy(T, reg, kind, var, z), f"{fname} BAYES (LOTO)", out)
            # degraded predictors: the board-like predictor quality
            y = np.log1p(T.nreal["1.0"].astype(float))
            for r_t in (0.40, 0.25):
                rng = np.random.default_rng(5210 + int(100 * r_t))
                vals = {"RANK q10 K25": [], "RANK q20 K100": [], "RANK q20 K400": [], "RANK q30 K400": [], "BAYES": [],
                        "TWO 10/25+20/400": []}
                for d in range(200):
                    zd = degrade(feats["S4"], y, r_t, rng)
                    R0 = T.rows[f"{reg}|base"]
                    for lab, Kv in (("RANK q10 K25", rank_policy(zd, 10, 25)), ("RANK q20 K100", rank_policy(zd, 20, 100)),
                                    ("RANK q20 K400", rank_policy(zd, 20, 400)), ("RANK q30 K400", rank_policy(zd, 30, 400)),
                                    ("TWO 10/25+20/400", two_level(zd, 10, 25, 20, 400))):
                        vals[lab].append(value(T.assemble(reg, kind, var, Kv), R0))
                    if d < 40:   # the Bayes action is slower; 40 draws
                        vals["BAYES"].append(value(T.assemble(reg, kind, var, bayes_policy(T, reg, kind, var, zd)), R0))
                for lab, v in vals.items():
                    v = np.asarray(v)
                    out[f"{reg}|{kind}|{var}|DEGRADED r{r_t} {lab}"] = dict(mean=float(v.mean()), q05=float(np.quantile(v, 0.05)),
                                                                          q95=float(np.quantile(v, 0.95)), p_le0=float(np.mean(v <= 0)))
                    print(f"    DEGRADED r {r_t:.2f} {lab:18s} mean {v.mean():+.5f}  [5-95 % over draws {np.quantile(v, 0.05):+.5f}, "
                          f"{np.quantile(v, 0.95):+.5f}]  P(<=0) {np.mean(v <= 0):.2f}")
    (L.WORK / f"policy521_s{seed}.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
