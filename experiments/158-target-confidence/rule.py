"""158 — does any build-time key predict `s_p*`, and what does using it actually buy?

WHAT IS MEASURED
----------------
`sse(s) = sum_p ||s_p d_p - delta_p||^2 / sum_p ||delta_p||^2` is SEPARABLE in `s_p`, so the
unconstrained optimum is `s_p* = <d_p, delta_p>/||d_p||^2` per row and the constrained-positive
optimum is `max(s_p*, 0)`. Every rule below produces a positive weight vector `w`, which is
then multiplied by ONE globally re-optimised scale `g` (free, and pds-neutral), and scored.

Because `g` is re-optimised, only the SHAPE of `w` matters -- a key is being asked to order
targets by how much of their predicted direction is real, nothing more.

LEAKAGE DISCIPLINE
------------------
Two surfaces. Every reported number is produced by a rule whose key choice, sign, gamma and
floor were selected on the OTHER surface. Both directions are reported. The `--shrink-pert-*`
weight form is production's `pert_shrink_vector`, reimplemented here only so the NaN policy is
explicit; `check_prod_form.py` asserts the two agree bit for bit.
"""
from __future__ import annotations

import itertools
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from oracle import (SURFACES, cos_agg, delivered, global_star, pds_of,  # noqa: E402
                    per_target_cos, per_target_star, sse_of)

WORK = os.path.join(HERE, "work")
OTHER = {"hct116": "hek293t", "hek293t": "hct116"}


# ------------------------------------------------------------------- weight form
def rank01(v: np.ndarray) -> np.ndarray:
    """097's ascending mid-rank on (0, 1); the best target -> ~1, the worst -> 0.5/n."""
    order = np.argsort(np.argsort(v, kind="stable"), kind="stable").astype(np.float64)
    # ties share a midrank, as `_descending_midrank` does
    out = np.empty_like(order)
    for val in np.unique(v):
        sel = v == val
        out[sel] = order[sel].mean()
    return (out + 0.5) / v.size


def weight(key: np.ndarray, gamma: float, floor: float, sign: float = 1.0) -> np.ndarray:
    """`z_p = max(rank01(sign*key), floor)**gamma`, production's form, `missing='neutral'`.

    `rank01` never returns 0, so `z_p > 0` for every finite floor -- which is exactly the
    condition under which the scaling is pds-exact (verified in `oracle.py`).
    """
    k = sign * np.asarray(key, dtype=np.float64)
    ok = np.isfinite(k)
    z = np.ones(k.size)
    if ok.sum() < 2 or np.unique(k[ok]).size < 2:
        return z
    z[ok] = np.maximum(rank01(k[ok]), floor) ** gamma
    z[~ok] = float(np.median(z[ok]))          # absence of evidence is not its opposite
    return z


def quantile_shrink(key: np.ndarray, q: float, floor: float, sign: float = 1.0) -> np.ndarray:
    """The cheap discrete version: the worst `q` fraction of targets keep `floor` of their
    amplitude, everybody else keeps 1. `floor > 0` keeps pds exact."""
    k = sign * np.asarray(key, dtype=np.float64)
    ok = np.isfinite(k)
    z = np.ones(k.size)
    if ok.sum() < 2:
        return z
    r = np.full(k.size, np.nan)
    r[ok] = rank01(k[ok])
    z[ok] = np.where(r[ok] <= q, floor, 1.0)
    return z


# ------------------------------------------------------------------- scoring one w
def score(w: np.ndarray, P: np.ndarray, T: np.ndarray) -> dict:
    Pw = w[:, None] * P
    g = global_star(Pw, T)
    s = g * w
    u = sse_of(s, P, T)
    return dict(sse=u, cos=cos_agg(u), g=g)


# --------------------------------------------------------------------- correlations
def spearman(a: np.ndarray, b: np.ndarray) -> float:
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 5:
        return np.nan
    return float(np.corrcoef(rank01(a[ok]), rank01(b[ok]))[0, 1])


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 5:
        return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


# ---------------------------------------------------------------- cluster bootstrap
def clusters(P: np.ndarray, n_clusters: int = 20, seed: int = 0) -> np.ndarray:
    """Group targets by the similarity of OUR OWN predicted direction (build-time, no truth).

    Targets whose predictions point the same way share whatever systematic error the donors
    carry, so resampling individual targets would understate the variance. Ward linkage on
    the unit-normalised delivered rows; deterministic, no truth touched.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    U = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    Z = linkage(U, method="ward")
    return fcluster(Z, t=n_clusters, criterion="maxclust") - 1


def boot_gain(w: np.ndarray, P: np.ndarray, T: np.ndarray, lab: np.ndarray,
              n: int = 2000, seed: int = 1) -> tuple[float, float, float]:
    """Paired cluster bootstrap of `sse(global) - sse(rule)`; both terms recomputed with the
    global scale re-optimised INSIDE each resample, so `g`'s own sampling error is included."""
    rng = np.random.default_rng(seed)
    uniq = np.unique(lab)
    idx = [np.flatnonzero(lab == c) for c in uniq]
    gains = np.empty(n)
    for i in range(n):
        pick = np.concatenate([idx[j] for j in rng.integers(0, len(idx), len(idx))])
        Pb, Tb, wb = P[pick], T[pick], w[pick]
        base = score(np.ones(len(pick)), Pb, Tb)["sse"]
        gains[i] = base - score(wb, Pb, Tb)["sse"]
    return float(np.percentile(gains, 2.5)), float(np.median(gains)), float(np.percentile(gains, 97.5))


def perm_null(key: np.ndarray, gamma: float, floor: float, sign: float,
              P: np.ndarray, T: np.ndarray, n: int = 2000, seed: int = 2,
              discrete: bool = False) -> tuple[float, np.ndarray]:
    """Shuffle the KEY across targets, holding the weight form and every hyperparameter fixed.
    Returns `(one-sided p, null gains)`."""
    rng = np.random.default_rng(seed)
    base = score(np.ones(len(key)), P, T)["sse"]
    obs = base - score(quantile_shrink(key, gamma, floor, sign) if discrete
                       else weight(key, gamma, floor, sign), P, T)["sse"]
    null = np.empty(n)
    for i in range(n):
        k = key[rng.permutation(len(key))]
        w = quantile_shrink(k, gamma, floor, sign) if discrete else weight(k, gamma, floor, sign)
        null[i] = base - score(w, P, T)["sse"]
    return float((np.sum(null >= obs) + 1) / (n + 1)), null


# ---------------------------------------------------------------------------- main
GAMMAS = (0.25, 0.5, 1.0, 2.0, 4.0)
FLOORS = (0.0, 0.1, 0.2, 0.3, 0.5)
QS = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50)
QFLOORS = (0.05, 0.15, 0.30, 0.50, 0.70)


def main() -> None:
    with open(os.path.join(WORK, "keys.json")) as fh:
        K = json.load(fh)

    S = {}
    for surface in SURFACES:
        P, T, ptar, aux = delivered(surface)
        assert [str(x) for x in ptar] == K[surface]["targets"]
        keys = {k: np.asarray(v, dtype=np.float64) for k, v in K[surface]["keys"].items()}
        keys = {k: v for k, v in keys.items() if np.unique(v[np.isfinite(v)]).size >= 2}
        S[surface] = dict(P=P, T=T, keys=keys,
                          s=per_target_star(P, T), cos=per_target_cos(P, T),
                          base=score(np.ones(P.shape[0]), P, T),
                          lab=clusters(P))

    names = sorted(set(S["hct116"]["keys"]) & set(S["hek293t"]["keys"]))
    print(f"{len(names)} keys carried by both surfaces\n")

    # ---------------------------------------------------------------- ceilings
    print("=" * 108)
    print("CEILINGS (what any positive per-target scale could reach)")
    print(f"  {'surface':10s} {'sse global':>11s} {'sse s*>=0':>11s} {'sse oracle':>11s} "
          f"{'cos glob':>9s} {'cos s*>=0':>10s} {'cos oracle':>10s}")
    ceil = {}
    for surface in SURFACES:
        d = S[surface]
        s = d["s"]
        pos = sse_of(np.maximum(s, 0.0), d["P"], d["T"])
        orc = sse_of(s, d["P"], d["T"])
        ceil[surface] = dict(base=d["base"]["sse"], pos=pos, oracle=orc)
        print(f"  {surface:10s} {d['base']['sse']:11.6f} {pos:11.6f} {orc:11.6f} "
              f"{cos_agg(d['base']['sse']):9.5f} {cos_agg(pos):10.5f} {cos_agg(orc):10.5f}")

    # ------------------------------------------------------ correlations to s_p*
    print("\n" + "=" * 108)
    print("KEY -> s_p* ASSOCIATION (both surfaces reported; nothing selected here)")
    print(f"  {'key':24s} {'sp(s*) hct':>11s} {'sp(s*) hek':>11s} "
          f"{'sp(cos) hct':>12s} {'sp(cos) hek':>12s} {'pe(s*) hct':>11s} {'pe(s*) hek':>11s}")
    corr = {}
    for k in names:
        row = {}
        for surface in SURFACES:
            d = S[surface]
            row[surface] = dict(sp_s=spearman(d["keys"][k], d["s"]),
                                sp_cos=spearman(d["keys"][k], d["cos"]),
                                pe_s=pearson(d["keys"][k], d["s"]))
        corr[k] = row
        print(f"  {k:24s} {row['hct116']['sp_s']:+11.3f} {row['hek293t']['sp_s']:+11.3f} "
              f"{row['hct116']['sp_cos']:+12.3f} {row['hek293t']['sp_cos']:+12.3f} "
              f"{row['hct116']['pe_s']:+11.3f} {row['hek293t']['pe_s']:+11.3f}")

    # ------------------------------------------------- cross-surface rule fitting
    print("\n" + "=" * 108)
    print("CROSS-SURFACE CONTINUOUS RULE  z = max(rank01(sign*key), floor)^gamma")
    print("  hyperparameters (sign, gamma, floor) chosen on the FIT surface, reported on EVAL")
    print(f"  {'key':24s} {'fit':>8s} {'eval':>8s} {'sign':>5s} {'gam':>5s} {'flr':>4s} "
          f"{'eval sse':>10s} {'d sse':>10s} {'eval cos':>9s} {'d cos %':>8s} {'own-best':>10s}")
    grid = list(itertools.product((1.0, -1.0), GAMMAS, FLOORS))
    table = {}
    for k in names:
        rec = {}
        for fit in SURFACES:
            ev = OTHER[fit]
            fd, ed = S[fit], S[ev]
            best = max(grid, key=lambda hp: fd["base"]["sse"]
                       - score(weight(fd["keys"][k], hp[1], hp[2], hp[0]), fd["P"], fd["T"])["sse"])
            sign, gam, flr = best
            r = score(weight(ed["keys"][k], gam, flr, sign), ed["P"], ed["T"])
            own_best = max(ed["base"]["sse"]
                           - score(weight(ed["keys"][k], hp[1], hp[2], hp[0]),
                                   ed["P"], ed["T"])["sse"] for hp in grid)
            rec[ev] = dict(sign=sign, gamma=gam, floor=flr, sse=r["sse"], cos=r["cos"],
                           d_sse=ed["base"]["sse"] - r["sse"],
                           d_cos_pct=100 * (r["cos"] / cos_agg(ed["base"]["sse"]) - 1),
                           own_best=own_best)
            print(f"  {k:24s} {fit:>8s} {ev:>8s} {sign:+5.0f} {gam:5.2f} {flr:4.1f} "
                  f"{r['sse']:10.6f} {ed['base']['sse'] - r['sse']:+10.6f} {r['cos']:9.5f} "
                  f"{100 * (r['cos'] / cos_agg(ed['base']['sse']) - 1):+7.1f}% {own_best:+10.6f}")
        table[k] = rec

    # ------------------------------------------------------- discrete quantile rule
    print("\n" + "=" * 108)
    print("CROSS-SURFACE DISCRETE RULE  worst q fraction -> floor, rest -> 1")
    print(f"  {'key':24s} {'fit':>8s} {'eval':>8s} {'sign':>5s} {'q':>5s} {'flr':>5s} "
          f"{'eval sse':>10s} {'d sse':>10s} {'eval cos':>9s} {'d cos %':>8s} {'own-best':>10s}")
    dgrid = list(itertools.product((1.0, -1.0), QS, QFLOORS))
    dtable = {}
    for k in names:
        rec = {}
        for fit in SURFACES:
            ev = OTHER[fit]
            fd, ed = S[fit], S[ev]
            best = max(dgrid, key=lambda hp: fd["base"]["sse"]
                       - score(quantile_shrink(fd["keys"][k], hp[1], hp[2], hp[0]),
                               fd["P"], fd["T"])["sse"])
            sign, q, flr = best
            r = score(quantile_shrink(ed["keys"][k], q, flr, sign), ed["P"], ed["T"])
            own_best = max(ed["base"]["sse"]
                           - score(quantile_shrink(ed["keys"][k], hp[1], hp[2], hp[0]),
                                   ed["P"], ed["T"])["sse"] for hp in dgrid)
            rec[ev] = dict(sign=sign, q=q, floor=flr, sse=r["sse"], cos=r["cos"],
                           d_sse=ed["base"]["sse"] - r["sse"],
                           d_cos_pct=100 * (r["cos"] / cos_agg(ed["base"]["sse"]) - 1),
                           own_best=own_best)
            print(f"  {k:24s} {fit:>8s} {ev:>8s} {sign:+5.0f} {q:5.2f} {flr:5.2f} "
                  f"{r['sse']:10.6f} {ed['base']['sse'] - r['sse']:+10.6f} {r['cos']:9.5f} "
                  f"{100 * (r['cos'] / cos_agg(ed['base']['sse']) - 1):+7.1f}% {own_best:+10.6f}")
        dtable[k] = rec

    out = dict(ceiling=ceil, corr=corr, continuous=table, discrete=dtable)
    with open(os.path.join(WORK, "rule.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {os.path.join(WORK, 'rule.json')}")

    # ------------------------------------------------- worst-surface ranking, both forms
    print("\n" + "=" * 108)
    print("WORST-SURFACE d sse (the repo's selection rule). Positive = mse improves.")
    print(f"  {'key':24s} {'cont worst':>11s} {'cont hct':>10s} {'cont hek':>10s} "
          f"{'disc worst':>11s} {'disc hct':>10s} {'disc hek':>10s}")
    rank = []
    for k in names:
        c = [table[k][s]["d_sse"] for s in SURFACES]
        d = [dtable[k][s]["d_sse"] for s in SURFACES]
        rank.append((max(min(c), min(d)), k, c, d))
    for _, k, c, d in sorted(rank, reverse=True):
        print(f"  {k:24s} {min(c):+11.6f} {c[0]:+10.6f} {c[1]:+10.6f} "
              f"{min(d):+11.6f} {d[0]:+10.6f} {d[1]:+10.6f}")


if __name__ == "__main__":
    main()
