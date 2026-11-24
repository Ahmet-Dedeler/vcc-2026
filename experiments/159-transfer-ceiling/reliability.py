"""159 §2 -- split-half reliability of a pseudobulk effect, as a function of cells per target.

DEFINITIONS (all in the scorer's own space: b = log1p(TS*P/sum P), TS = 5e4)

  full effect      d_t   = b(all cells of t)      - b(all control cells)
  half effect      e_t^i = b(m cells of t, set i) - b(control half i)
  split-half r     r(m)  = cos(e^1, e^2), averaged over disjoint pairings and over targets

Under d_obs = d_true + eps with eps independent of d_true and of the other half,
  E[r(m)] ~= ||d||^2 / (||d||^2 + sigma^2(m))   =: rho(m)
so r(m) IS the squared correlation of an m-cell effect with the truth, and rho at any other
depth follows from the fitted noise-vs-depth law (see `fit`).

TWO CONTROL VARIANTS, because they bracket a real bias:
  `disjoint`  each half scored against a DISJOINT half of the control -> independent noise;
              this is the estimator the disattenuation formula assumes.
  `shared`    both halves scored against the FULL control -> shared control noise is
              common-mode and inflates r.  Reported so the inflation is visible, not assumed.

TWO SPACES, because they answer different questions:
  `full`      the raw effect vector -- what `pds_cosine` and `expr_mse` actually read.
  `resid`     the effect minus the line's across-target mean response at the SAME depth --
              the target-SPECIFIC component, the only part that can earn score (003: the mean
              response IS the metric's zero).  Much lower SNR, so much more attenuated.

KNOWN-ANSWER ROW: h1-test / h1-full carry exactly 400 cells per target -- the challenge's own
count -- so the extrapolation the 150-cell Orion surfaces need can be validated against a
directly measured value (fit on m <= half the range, predict the largest measured m).
"""
from __future__ import annotations

import itertools
import json
import os
import sys

import numpy as np

TS = 5.0e4
SEED = 20260902
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")


def prof(P):
    """counts -> b = log1p(TS * P / sum P), row-wise."""
    P = np.asarray(P, dtype=np.float64)
    if P.ndim == 1:
        return np.log1p(TS * P / P.sum())
    return np.log1p(TS * P / P.sum(1, keepdims=True))


def unit(A):
    return A / np.maximum(np.linalg.norm(A, axis=-1, keepdims=True), 1e-300)


def load(name):
    z = np.load(os.path.join(WORK, f"shards-{name}.npz"), allow_pickle=True)
    targets = np.asarray([str(x) for x in z["targets"]])
    genes = np.asarray([str(x) for x in z["genes"]])
    gpos = {g: i for i, g in enumerate(genes)}
    drop = {gpos[t] for t in targets if t in gpos}
    mask = np.asarray([i for i in range(len(genes)) if i not in drop])
    return dict(targets=targets, genes=genes, mask=mask,
                tsum=z["tsum"].astype(np.float64), tcell=z["tcell"],
                csum=z["csum"].astype(np.float64), ccell=z["ccell"],
                kt=int(z["kt"]), kc=int(z["kc"]))


def pairings(k, j, n_pair, rng):
    """n_pair disjoint (j-shard, j-shard) index pairs out of k shards."""
    all_sets = list(itertools.combinations(range(k), j))
    out, seen = [], set()
    tries = 0
    while len(out) < n_pair and tries < 60 * n_pair:
        tries += 1
        s1 = all_sets[rng.integers(len(all_sets))]
        rest = [x for x in range(k) if x not in s1]
        if len(rest) < j:
            continue
        s2 = tuple(sorted(rng.choice(rest, size=j, replace=False).tolist()))
        key = tuple(sorted([s1, s2]))
        if key in seen:
            continue
        seen.add(key)
        out.append((list(s1), list(s2)))
    return out


def curve(d, n_pair=24, verbose=True):
    """r(m) for both control variants and both spaces, per target and aggregated."""
    rng = np.random.default_rng(SEED + 1)
    tsum, csum, mask = d["tsum"], d["csum"], d["mask"]
    kt, kc = d["kt"], d["kc"]
    nt = tsum.shape[0]
    cells_per_shard = float(d["tcell"].mean())

    b_c_full = prof(csum.sum(0))[mask]
    b_c1 = prof(csum[: kc // 2].sum(0))[mask]
    b_c2 = prof(csum[kc // 2:].sum(0))[mask]

    d_full = prof(tsum.sum(1))[:, mask] - b_c_full[None, :]
    m_line_full = d_full.mean(0)
    ctrl_gap = float(np.linalg.norm(b_c1 - b_c2) ** 2)

    # DEPTH BIAS.  log1p is concave, so E[b(m cells)] != b(infinite cells): a shallow
    # pseudobulk carries a systematic offset that is a function of m and is nearly the SAME
    # vector for every target.  Both halves carry it, so it is common-mode and inflates the
    # split-half cosine -- far more in the low-SNR residual space than in the full space.  A
    # first pass that centred on the FULL-depth mean made resid r NON-MONOTONE in m, which is
    # that bias's signature.  Centring each side on the mean over targets of THAT side at THAT
    # depth removes it exactly (the bias is common across targets) and shares no noise between
    # sides (the two means come from disjoint cells).
    rows = []
    for j in range(1, kt // 2 + 1):
        m = j * cells_per_shard
        pp = pairings(kt, j, n_pair, rng)
        acc = {k: np.zeros(nt) for k in ("full_disjoint", "full_shared",
                                         "resid_disjoint", "resid_shared")}
        gap = np.zeros(nt)
        bias = np.zeros(len(mask))
        for s1, s2 in pp:
            p1 = prof(tsum[:, s1].sum(1))[:, mask]
            p2 = prof(tsum[:, s2].sum(1))[:, mask]
            gap += np.linalg.norm(p1 - p2, axis=1) ** 2
            for var, c1, c2 in (("disjoint", b_c1, b_c2), ("shared", b_c_full, b_c_full)):
                e1, e2 = p1 - c1[None, :], p2 - c2[None, :]
                acc[f"full_{var}"] += (unit(e1) * unit(e2)).sum(1)
                r1 = e1 - e1.mean(0, keepdims=True)
                r2 = e2 - e2.mean(0, keepdims=True)
                acc[f"resid_{var}"] += (unit(r1) * unit(r2)).sum(1)
                if var == "disjoint":
                    bias += 0.5 * (e1.mean(0) + e2.mean(0))
        for k in acc:
            acc[k] /= len(pp)
        gap /= len(pp)
        bias /= len(pp)
        rows.append(dict(j=j, m=m, n_pair=len(pp),
                         per_target={k: v.copy() for k, v in acc.items()},
                         target_gap=gap.copy(),
                         mean_norm=float(np.linalg.norm(bias)),
                         mean_cos_to_full=float(
                             bias @ m_line_full / max(np.linalg.norm(bias) *
                                                      np.linalg.norm(m_line_full), 1e-300))))
        if verbose:
            print(f"   m={m:6.1f}/side  pairs={len(pp):3d}  "
                  f"full/disj {acc['full_disjoint'].mean():.4f}  "
                  f"full/shar {acc['full_shared'].mean():.4f}  "
                  f"resid/disj {acc['resid_disjoint'].mean():.4f}  "
                  f"resid/shar {acc['resid_shared'].mean():.4f}  "
                  f"|mean| {rows[-1]['mean_norm']:6.3f} "
                  f"(cos to full-depth mean {rows[-1]['mean_cos_to_full']:.4f})")
    return rows, dict(cells_per_shard=cells_per_shard, ctrl_gap=ctrl_gap,
                      ctrl_cells_half=float(d["ccell"][: kc // 2].sum()),
                      n_targets=nt, n_scored_genes=int(len(mask)),
                      mean_norm_full_depth=float(np.linalg.norm(m_line_full)),
                      d_full_norm_median=float(np.median(np.linalg.norm(d_full, axis=1))),
                      resid_norm_median=float(np.median(
                          np.linalg.norm(d_full - m_line_full[None, :], axis=1))))


def fit(ms, rs, model="power"):
    """Fit the noise-vs-depth law to a measured split-half curve.  Returns (par, R2, f(n)).

    `linear` : 1/r - 1 = a/m + c   -- what a pure 1/n Poisson argument predicts.
    `power`  : 1/r - 1 = a/m^p     -- what the data shows.  log1p is not variance-stabilising
               at these depths, so the effective exponent is 0.72 on H1 and ~0.96-1.05 on
               Orion.  The extrapolation check in `main` decides which is used.
    """
    ms, rs = np.asarray(ms, float), np.asarray(rs, float)
    ok = (rs > 1e-6) & (rs < 1 - 1e-9)
    y = 1.0 / rs[ok] - 1.0
    if model == "linear":
        X = np.column_stack([1.0 / ms[ok], np.ones(ok.sum())])
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        pred = X @ coef
        par = {"model": "linear", "a": float(coef[0]), "c": float(coef[1])}

        def f(n):
            return 1.0 / (1.0 + par["a"] / np.asarray(n, float) + 0.5 * par["c"])
    else:
        X = np.column_stack([np.ones(ok.sum()), np.log(ms[ok])])
        coef, *_ = np.linalg.lstsq(X, np.log(y), rcond=None)
        pred = np.exp(X @ coef)
        par = {"model": "power", "a": float(np.exp(coef[0])), "p": float(-coef[1])}

        def f(n):
            return 1.0 / (1.0 + par["a"] * np.asarray(n, float) ** (-par["p"]))
    ss = 1.0 - ((y - pred) ** 2).sum() / max(((y - y.mean()) ** 2).sum(), 1e-30)
    return par, float(ss), f


KEYS = ("full_disjoint", "full_shared", "resid_disjoint", "resid_shared")


def main() -> None:
    names = sys.argv[1:] or ["hct116", "hek293t", "h1test"]
    out = {}
    for name in names:
        p = os.path.join(WORK, f"shards-{name}.npz")
        if not os.path.exists(p):
            print(f"[{name}] no shards, skipping")
            continue
        d = load(name)
        print(f"\n=== {name}: {d['tsum'].shape[0]} targets, {len(d['mask'])} scored genes, "
              f"{d['kt']} shards of {d['tcell'].mean():.1f} cells, "
              f"control {d['ccell'].sum():.0f} cells in {d['kc']} shards")
        rows, meta = curve(d)
        ms = [r["m"] for r in rows]
        rec = {"meta": meta, "seed": SEED, "m": ms,
               "mean_norm": [r["mean_norm"] for r in rows]}

        # ---- extrapolation check first: it decides which law we trust ------------------
        half = max(2, len(ms) // 2)
        print(f"\n   EXTRAPOLATION CHECK  (fit on m <= {ms[half - 1]:.0f}, "
              f"predict the measured m = {ms[-1]:.0f})")
        rec["extrap_check"] = {}
        for key in ("full_disjoint", "resid_disjoint"):
            agg = np.array([r["per_target"][key] for r in rows]).mean(1)
            line = []
            for model in ("linear", "power"):
                _, _, f = fit(ms[:half], agg[:half], model)
                pr = float(f(ms[-1]))
                rec["extrap_check"].setdefault(key, {})[model] = {
                    "predicted": pr, "measured": float(agg[-1]), "ratio": pr / agg[-1]}
                line.append(f"{model} {pr:.4f} (x{pr / agg[-1]:.3f})")
            print(f"     {key:<16} measured {agg[-1]:.4f}   " + "   ".join(line))

        print(f"\n   {'space/variant':<16} {'model':>6} {'a':>9} {'p|c':>8} {'R2':>7} "
              f"{'r(meas)':>8} {'rho150':>8} {'rho400':>8} {'rho_inf':>8}")
        for key in KEYS:
            per = np.array([r["per_target"][key] for r in rows])     # (n_m, n_targets)
            agg = per.mean(1)
            rec[key] = {"r_mean": agg.tolist(),
                        "r_median": np.median(per, axis=1).tolist(),
                        "r_p10": np.percentile(per, 10, axis=1).tolist(),
                        "r_p90": np.percentile(per, 90, axis=1).tolist(),
                        "per_target_largest_m": per[-1].tolist()}
            for model in ("linear", "power"):
                par, r2, f = fit(ms, agg, model)
                rec[key][model] = {**par, "r2": r2,
                                   "rho_150": float(f(150)), "rho_400": float(f(400)),
                                   "rho_inf": float(f(1e12))}
                sec = par.get("p", par.get("c"))
                print(f"   {key:<16} {model:>6} {par['a']:9.2f} {sec:8.4f} {r2:7.4f} "
                      f"{agg[-1]:8.4f} {f(150):8.4f} {f(400):8.4f} {f(1e12):8.4f}")
        out[name] = rec

    fn = os.path.join(WORK, "reliability.json")
    prev = {}
    if os.path.exists(fn):
        with open(fn) as fh:
            prev = json.load(fh)
    prev.update(out)
    with open(fn, "w") as fh:
        json.dump(prev, fh, indent=1)
    print(f"\nwrote {fn}")


if __name__ == "__main__":
    main()
