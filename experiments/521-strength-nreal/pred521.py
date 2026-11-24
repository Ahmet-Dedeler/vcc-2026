"""521 step 2 — the n_real predictor, out of sample, with its full predictive distribution.

Strength statistic (518's `strength_for`, parameterised): per bank, genes > 10 CPM in the bank's control minus the 300
A/B/C panel genes and the scored targets; per-gene noise scale = 1.4826 MAD over the bank's A/B/C-panel rows; per-target
noise = median |z| / 0.6745; strength = log1p #{g : |z| > thr * noise}. Variants: thr 3 / 4 / 5; CD4 T cells as the
robust statistic and as DESeq2's own nDE count (adj p < 0.05, mean log1p over 3 conditions).

Surfaces (truth = the scorer's rules; H1 data never a source when H1 is the truth):
  H1       85 targets, 323's reference (400 cells, A/B/C power). Sources: K562-GWPS, Orion HCT116, Orion HEK293T,
           VIPerturb K562 (+ essential screens, which cover only H1's essential targets). Also the 67 non-essential
           targets alone (the A/B/C panel is essential-free).
  folds    480's gates (Wilcoxon vs the line's own controls, BH, 150 cells, ~12k UMI: BELOW A/B/C power):
           K562-GWPS (267 panel targets), Orion HCT116 thinned (132), Orion HEK293T thinned (132).
           Sources exclude the held-out line and any same-line bank (VIPerturb on K562); 'xs' also drops the Orion
           sibling (same study, shared guides).
Reported per feature: coverage, Pearson/Spearman with log1p n_real, the LOTO linear calibration's residual sd and its
multiplicative error exp(sd), LOTO R^2 vs a LOTO constant, AUC(negligible: n_real < 20) and AUC at the surface's own
41st percentile (H1's negligible base rate), and interval coverage of the LOTO Gaussian predictive distribution.

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/pred521.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402
from scipy.stats import norm, spearmanr  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
PANEL = [l.split(",")[0] for l in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()]
PSET = set(PANEL)
BANKS = {"k562": ROOT / "data/donors/k562_gwps", "hct116": ROOT / "data/donors/orion_hct116",
         "hek293t": ROOT / "data/donors/orion_hek293t", "vip": ROOT / "data/donors/viperturb_k562",
         "kolf2": ROOT / "data/external/kolf2", "k562ess": ROOT / "data/donors/k562_essential",
         "rpe1": ROOT / "data/donors/rpe1", "jurkat": ROOT / "data/donors/nadig_jurkat",
         "hepg2": ROOT / "experiments/166-donor-acquisition/donors/nadig_hepg2"}
CD4 = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
CD4C = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
GATES = {"k562": ROOT / "experiments/480-nmae-instrument/work/gate-k562_gwps.npz",
         "hct116": ROOT / "experiments/480-nmae-instrument/work/gate-orion_hct116.npz",
         "hek293t": ROOT / "experiments/480-nmae-instrument/work/gate-orion_hek293t.npz"}
REF = ROOT / "experiments/323-current-family-validation/work/reference.npz"


def read_list(p):
    return [l.strip() for l in p.read_text().splitlines() if l.strip()]


def strength(path, targets, thr=4.0):
    """518's statistic with its threshold as a parameter; NaN where the bank does not carry the target."""
    perts = read_list(path / "perts.txt")
    genes = read_list(path / "genes.txt")
    b = np.load(path / "b_ctrl.npy").astype(np.float64)
    ex = np.expm1(b)
    cpm = 1e6 * ex / ex.sum()
    tset = set(targets)
    keep = (cpm > 10.0) & np.asarray([g not in PSET and g not in tset for g in genes])
    pi = {p: i for i, p in enumerate(perts)}
    delta = np.load(path / "delta.npy", mmap_mode="r")
    prow = sorted(pi[t] for t in PANEL if t in pi)
    if len(prow) < 50:                       # an essential screen: no panel rows; scale over its first 300 rows
        prow = list(range(min(300, len(perts))))
    P = np.asarray(delta[np.asarray(prow)], np.float64)[:, keep]
    sc = 1.4826 * np.median(np.abs(P - np.median(P, 0)), 0)
    ok = sc > 0
    cov = [t for t in targets if t in pi]
    out = {}
    if cov:
        rows = np.asarray([pi[t] for t in cov])
        o = np.argsort(rows)
        R = np.empty((len(cov), int(keep.sum())))
        R[o] = np.asarray(delta[rows[o]], np.float64)[:, keep]
        Z = R[:, ok] / sc[ok]
        nz = np.median(np.abs(Z), 1) / 0.6745
        cnt = (np.abs(Z) > thr * nz[:, None]).sum(1)
        out = dict(zip(cov, np.log1p(cnt).tolist()))
    return np.asarray([out.get(t, np.nan) for t in targets])


def cd4_feats(targets):
    rob, nde = [], []
    for c in CD4C:
        d = CD4 / c
        rob.append(strength(d, targets))
        perts, genes = read_list(d / "perts.txt"), read_list(d / "genes.txt")
        nonp = np.asarray([g not in PSET for g in genes])
        pi = {p: i for i, p in enumerate(perts)}
        adj = np.load(d / "adj_p_value.npy", mmap_mode="r")
        v = np.full(len(targets), np.nan)
        for k, t in enumerate(targets):
            if t in pi:
                v[k] = np.log1p(float((np.asarray(adj[pi[t]], float)[nonp] < 0.05).sum()))
        nde.append(v)
    with np.errstate(all="ignore"):
        return np.nanmean(np.vstack(rob), 0), np.nanmean(np.vstack(nde), 0)


def crank(v):
    v = np.asarray(v, float)
    ok = np.isfinite(v)
    out = np.full(v.size, np.nan)
    if ok.sum() > 1:
        out[ok] = np.argsort(np.argsort(v[ok])) / (ok.sum() - 1) - 0.5
    return out


def meanrank(cols):
    with np.errstate(all="ignore"):
        return np.nanmean(np.vstack([crank(c) for c in cols]), 0)


def auc(x, pos):
    """P(x of a positive > x of a negative), ties 1/2; pos = the NON-negligible class."""
    xp, xn = x[pos], x[~pos]
    if xp.size == 0 or xn.size == 0:
        return np.nan
    return float(np.mean((xp[:, None] > xn[None, :]) + 0.5 * (xp[:, None] == xn[None, :])))


def loto(x, y):
    """leave-one-target-out linear calibration y ~ a + b x; returns (pred, sd_loo) on the covered targets."""
    n = x.size
    pred = np.empty(n)
    for i in range(n):
        m = np.arange(n) != i
        b, a = np.polyfit(x[m], y[m], 1)
        pred[i] = a + b * x[i]
    return pred


def evaluate(x, nreal, label):
    ok = np.isfinite(x)
    x, n = x[ok], nreal[ok]
    y = np.log1p(n)
    r = float(np.corrcoef(x, y)[0, 1])
    rs = float(spearmanr(x, y).correlation)
    pred = loto(x, y)
    res = y - pred
    const = np.array([(y.sum() - y[i]) / (y.size - 1) for i in range(y.size)])
    sd, sd0 = float(res.std(ddof=1)), float((y - const).std(ddof=1))
    r2 = float(1 - (res ** 2).sum() / ((y - const) ** 2).sum())
    q41 = np.quantile(n, 0.41)
    # LOTO Gaussian predictive distribution: coverage of central intervals (sigma from the other targets' residuals)
    cov = {}
    for lev in (0.5, 0.8, 0.9):
        zc = norm.ppf(0.5 + lev / 2)
        hits = []
        for i in range(y.size):
            s_i = np.sqrt(np.delete(res, i).var(ddof=1))
            hits.append(abs(y[i] - pred[i]) <= zc * s_i)
        cov[str(lev)] = float(np.mean(hits))
    # heteroscedasticity: residual sd in the bottom / top half of the prediction
    lo = pred <= np.median(pred)
    row = dict(label=label, n=int(ok.sum()), pearson=r, spearman=rs, loto_resid_sd=sd, const_sd=sd0,
               mult_err=float(np.exp(sd)), loto_r2=r2, auc_neg20=auc(x, n >= 20), auc_q41=auc(x, n > q41),
               n_neg20=int((n < 20).sum()), cover=cov, sd_lowhalf=float(res[lo].std(ddof=1)),
               sd_highhalf=float(res[~lo].std(ddof=1)))
    print(f"  {label:34s} n {row['n']:3d}  r {r:+.3f}  rho {rs:+.3f}  LOTO sd {sd:.3f} (const {sd0:.3f}; x{np.exp(sd):.1f})"
          f"  R2 {r2:+.3f}  AUC<20 {row['auc_neg20']:.3f}  AUCq41 {row['auc_q41']:.3f}  cover50/80/90 "
          f"{cov['0.5']:.2f}/{cov['0.8']:.2f}/{cov['0.9']:.2f}  sd lo/hi {row['sd_lowhalf']:.2f}/{row['sd_highhalf']:.2f}")
    return row


def h1_surface():
    ref = np.load(REF)
    T = [str(t) for t in ref["targets"]]
    keep = ref["keep"].astype(bool)
    genes = ref["genes"].astype(str)[keep]
    q = np.asarray(ref["q"], float)
    nreal = np.array([int(((q[k] < 0.05) & (genes != t)).sum()) for k, t in enumerate(T)])
    ess = set(read_list(BANKS["k562ess"] / "perts.txt"))
    return T, nreal, np.array([t in ess for t in T])


def fold_surface(name):
    z = np.load(GATES[name])
    T = [str(t) for t in z["targets"]]
    g = z["genes"].astype(str)
    q = np.asarray(z["q"], float)
    nreal = np.array([int(((q[k] < 0.05) & (g != t)).sum()) for k, t in enumerate(T)])
    return T, nreal


def main():
    out = {"h1": {}, "folds": {}}
    cache = {}
    # ---------------------------------------------------------------- H1
    T, nreal, ess = h1_surface()
    print(f"H1: 85 targets, n_real median {np.median(nreal):.0f}, <20 {int((nreal < 20).sum())}, 0 {int((nreal == 0).sum())}, "
          f"essential {int(ess.sum())}")
    F = {}
    for thr in (3.0, 4.0, 5.0):
        for b in ("k562", "hct116", "hek293t", "vip", "k562ess", "rpe1", "jurkat", "hepg2"):
            F[f"{b}@{thr:g}"] = strength(BANKS[b], T, thr)
    cache["h1"] = {k: v.tolist() for k, v in F.items()}
    rows = []
    print("single banks:")
    for thr in (3.0, 4.0, 5.0):
        for b in ("k562", "hct116", "hek293t", "vip"):
            rows.append(evaluate(F[f"{b}@{thr:g}"], nreal, f"{b} thr{thr:g}"))
    for b in ("k562ess", "rpe1", "jurkat", "hepg2"):
        rows.append(evaluate(F[f"{b}@4"], nreal, f"{b} thr4 (essential screen)"))
    print("combinations (mean centred rank over the banks carrying a target):")
    combos = {"S3 = k562+hct116+hek293t (nreal521)": ("k562", "hct116", "hek293t"),
              "S4 = S3 + vip": ("k562", "hct116", "hek293t", "vip"),
              "k562+vip (two K562 labs)": ("k562", "vip"),
              "S4 + essential screens": ("k562", "hct116", "hek293t", "vip", "k562ess", "rpe1", "jurkat", "hepg2")}
    Z = {}
    for thr in (3.0, 4.0, 5.0):
        for nm, bs in combos.items():
            z = meanrank([F[f"{b}@{thr:g}"] for b in bs])
            Z[f"{nm} thr{thr:g}"] = z
            rows.append(evaluate(z, nreal, f"{nm[:22]} thr{thr:g}"))
    # the mean of the thresholds (a threshold-robust statistic)
    z = meanrank([meanrank([F[f"{b}@{t:g}"] for t in (3.0, 4.0, 5.0)]) for b in ("k562", "hct116", "hek293t", "vip")])
    Z["S4 thr-avg"] = z
    rows.append(evaluate(z, nreal, "S4 thr-avg (3/4/5)"))
    # K562-weighted: the best single source counted twice
    z = meanrank([F["k562@4"], F["k562@4"], F["hct116@4"], F["hek293t@4"], F["vip@4"]])
    Z["S4 k562x2 thr4"] = z
    rows.append(evaluate(z, nreal, "S4 with K562 x2, thr4"))
    print("non-essential H1 targets only (the A/B/C panel is essential-free):")
    ne = ~ess
    for nm in ("S3 = k562+hct116+hek293t (nreal521) thr4", "S4 = S3 + vip thr4", "S4 thr-avg"):
        zz = Z[nm].copy()
        zz[~ne] = np.nan
        rows.append(evaluate(zz, nreal, f"[non-ess] {nm[:24]}"))
    zz = F["k562@4"].copy(); zz[~ne] = np.nan
    rows.append(evaluate(zz, nreal, "[non-ess] k562 thr4"))
    out["h1"]["rows"] = rows
    np.savez(WORK / "pred521_h1.npz", targets=np.array(T), nreal=nreal, ess=ess,
             **{k.replace(" ", "_").replace("=", "").replace("+", "p").replace("(", "").replace(")", ""): v for k, v in Z.items()},
             **{k.replace("@", "_t"): v for k, v in F.items()})
    # ---------------------------------------------------------------- cancer folds
    src_all = {"k562": ("hct116", "hek293t", "kolf2", "cd4rob", "cd4nde"),
               "hct116": ("k562", "vip", "hek293t", "kolf2", "cd4rob", "cd4nde"),
               "hek293t": ("k562", "vip", "hct116", "kolf2", "cd4rob", "cd4nde")}
    src_xs = {"k562": ("hct116", "hek293t", "kolf2", "cd4rob", "cd4nde"),
              "hct116": ("k562", "vip", "kolf2", "cd4rob", "cd4nde"),
              "hek293t": ("k562", "vip", "kolf2", "cd4rob", "cd4nde")}
    for fold in ("k562", "hct116", "hek293t"):
        Tf, nf = fold_surface(fold)
        print(f"\nfold {fold}: {len(Tf)} targets, n_real median {np.median(nf):.0f}, <20 {int((nf < 20).sum())}, "
              f"0 {int((nf == 0).sum())}, q41 {np.quantile(nf, 0.41):.0f}")
        G = {}
        for b in ("k562", "hct116", "hek293t", "vip", "kolf2"):
            G[b] = strength(BANKS[b], Tf, 4.0)
        G["cd4rob"], G["cd4nde"] = cd4_feats(Tf)
        cache[f"fold_{fold}"] = {k: v.tolist() for k, v in G.items()}
        frows = []
        for b in src_all[fold]:
            frows.append(evaluate(G[b], nf, f"{b} alone"))
        # the H1-validated combination's analogue: production sources available on this fold (+ vip)
        prod = [b for b in ("k562", "hct116", "hek293t", "vip") if b in src_all[fold]]
        frows.append(evaluate(meanrank([G[b] for b in prod]), nf, "prod sources (S4 minus held-out)"))
        frows.append(evaluate(meanrank([G[b] for b in src_all[fold]]), nf, "all sources (meta)"))
        frows.append(evaluate(meanrank([G[b] for b in src_xs[fold]]), nf, "cross-study only (xs)"))
        frows.append(evaluate(meanrank([G[b] for b in ("kolf2", "cd4rob", "cd4nde")]), nf, "kolf2 + CD4 (no cancer line)"))
        out["folds"][fold] = frows
    (WORK / "pred521.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "strength521_cache.json").write_text(json.dumps(cache))


if __name__ == "__main__":
    main()
