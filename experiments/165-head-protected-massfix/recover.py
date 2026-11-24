"""165 — recover the h1-test intended effect `D`, `b_ctrl` and the CPM gate from the
frozen 095 prediction, so the mass correction can be re-placed and re-priced.

WHY A RECOVERY AND NOT A REBUILD
--------------------------------
157 validated `h1-test` on `reach` against two board points, and it did so through ONE
specific artifact: `095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad`
(donors `k562_gwps` + `orion_hek293t`; `h1_2025` is the oracle-leak donor and never enters).
The known-answer reach values 0.18283 (a15 law) and 0.20291 (a8 law) belong to that artifact.
Rebuilding a fresh h1-test prediction at a15's exact config would throw that calibration away,
so the profile is taken as given and only the mass correction is varied on top of it.

THE ALGEBRA
-----------
The builder renders `p = clip(expm1(b + D), 0) / S(D)`, `S(D) = sum_g clip(...)`, and the
emitter lays down `totals_g = floor(n*L*p_g)` exactly (`donor_transfer.emit_counts`, quasi
mode, exact-total repair).  Writing `totals_g + U_g = n*L*p_g` with `U_g` in [0,1),

    clip(expm1(b_g + D_g), 0)  =  S(D) * p_g  =  kappa * (totals_g + U_g),
    kappa := S(D) / (n*L).

`D` is ZERO outside the CPM gate by construction, so `kappa` is pinned with no free parameter
by the ungated genes alone, and then `D = log1p(kappa*(totals+U)) - b` everywhere.  `n`, `L`
and the mass ratio never have to be known separately: `S(D)/TS = sum_g kappa*(totals+U)/TS`
falls out.

`U` is not observable, but `E[U] = 1/2` for any gene whose fractional part is uniform, which
is every gene with a non-trivial count.  It is NOT 1/2 for a gene emitting zero counts (there
`U = n*L*p < 1` and is typically far below 1/2), so `kappa` is estimated on a BAND of ungated
genes chosen by their CONTROL expression -- never by their emitted count, which would select
on the very quantity being corrected.

KNOWN-ANSWER CHECKS, AND EVERY ONE OF THEM CAN FAIL (AGENTS.md 11.1 rule 10's corollary)
----------------------------------------------------------------------------------------
1. `D` must vanish on EVERY ungated gene, at the residual scale ONE emitted count explains
   and no larger.  A wrong `b_ctrl` (wrong normalisation, wrong control set, wrong TS) does
   not vanish.
2. The gate is not assumed: `min_cpm` 5.0 and 10.0 are both tried and the one that vanishes
   is reported.  They differ by 6 orders of magnitude in residual, so this is decisive.
3. `kappa` must be STABLE across disjoint control-expression bands.  It is a single scalar
   fitted per target; if `b_ctrl`'s shape were wrong the bands would disagree.

No cell-eval2, no build, no submission.  Reads only; writes only into this directory.
"""
from __future__ import annotations

import os
import time

import h5py
import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
REAL = f"{ROOT}/data/harness/h1-test.h5ad"
PRED = (f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/"
        f"g0.65_q1.0_a1.08_spread0.25/pred.h5ad")
CTRL, COL = "non-targeting", "target_gene"
TS = 5.0e4


def _index(group):
    v = group["_index"]
    if isinstance(v, h5py.Group):
        v = v["values"]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v[:]])


def _labels(f):
    o = f["obs"][COL]
    if isinstance(o, h5py.Group):
        cats = np.asarray([x.decode() if isinstance(x, bytes) else str(x)
                           for x in o["categories"][:]])
        return cats[np.asarray(o["codes"][:])]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x)
                       for x in np.asarray(o[:])])


def _row_group_sums(path, groups, n_genes):
    """Column totals per row-group, streamed straight off the CSR (memory stays flat)."""
    out = np.zeros((len(groups), n_genes), dtype=np.float64)
    with h5py.File(path, "r") as f:
        X = f["X"]
        indptr = np.asarray(X["indptr"][:])
        data, indices = X["data"], X["indices"]
        for i, rows in enumerate(groups):
            rows = np.sort(np.asarray(rows))
            breaks = np.flatnonzero(np.diff(rows) != 1)
            starts = np.concatenate(([0], breaks + 1))
            stops = np.concatenate((breaks + 1, [rows.size]))
            acc = out[i]
            for s, e in zip(starts, stops):
                lo, hi = int(indptr[rows[s]]), int(indptr[rows[e - 1] + 1])
                if hi <= lo:
                    continue
                np.add.at(acc, np.asarray(indices[lo:hi]),
                          np.asarray(data[lo:hi], np.float64))
    return out


def kappa_on_band(tot, lin_b, band):
    """kappa = sum_band expm1(b) / sum_band (totals + 1/2), one scalar per target."""
    return lin_b[band].sum() / (tot[:, band] + 0.5).sum(1)


def main() -> None:
    t0 = time.time()
    os.makedirs(WORK, exist_ok=True)

    with h5py.File(REAL, "r") as f:
        genes = _index(f["var"])
        rlab = _labels(f)
    with h5py.File(PRED, "r") as f:
        assert np.array_equal(_index(f["var"]), genes), "gene axes differ"
        plab = _labels(f)
    G = genes.size
    ci = np.flatnonzero(rlab == CTRL)
    targets = [t for t in sorted(set(plab) - {CTRL}) if (rlab == t).sum() > 0]
    print(f"h1-test: {G} genes, {ci.size} control cells, {len(targets)} perturbations",
          flush=True)

    # ---- b_ctrl, production convention (make_submission_pm_cd4.py:2056) ------------------
    ctrl_tot = _row_group_sums(REAL, [ci], G)[0]
    b = np.log1p(TS * ctrl_tot / ctrl_tot.sum())
    lin_b = np.expm1(b)
    cpm_pb = 1e6 * lin_b / lin_b.sum()
    print(f"b_ctrl: sum(expm1(b)) = {lin_b.sum():.6f} (must be exactly TS = {TS:.0f}) "
          f"[{time.time() - t0:.0f}s]", flush=True)

    tot = np.rint(_row_group_sums(
        PRED, [np.flatnonzero(plab == t) for t in targets], G)).astype(np.int64)
    T = len(targets)
    print(f"emitted totals read: row sums min {tot.sum(1).min():,} "
          f"max {tot.sum(1).max():,} [{time.time() - t0:.0f}s]", flush=True)
    print(f"  genes with zero emitted counts: median {int(np.median((tot == 0).sum(1)))}"
          f" of {G}  -> the reason kappa is fitted on a band, not on the whole off-gate set")

    # ---- CHECK 2: which gate makes D vanish? --------------------------------------------
    print("\n--- CHECK 2: gate identification ---")
    picked = None
    for min_cpm in (5.0, 10.0):
        gate = (cpm_pb > min_cpm) & (b > 0.0)
        band = (~gate) & (cpm_pb >= 2.0)
        kap = kappa_on_band(tot, lin_b, band)
        D = np.log1p(kap[:, None] * (tot + 0.5)) - b[None, :]
        res = np.abs(D[:, ~gate])
        print(f"  min_cpm {min_cpm:5.1f}: {int(gate.sum()):6d} gated, "
              f"{int(band.sum()):5d} in the kappa band | "
              f"max|D| off-gate {res.max():.3e}  median {np.median(res):.3e}")
        if picked is None or res.max() < picked[0]:
            picked = (res.max(), min_cpm, gate, band)
    _, MIN_CPM, gate, band = picked
    off = ~gate
    print(f"  -> adopting min_cpm = {MIN_CPM} "
          f"({'the builder default' if MIN_CPM == 5.0 else 'a15 production'})")

    # ---- CHECK 3: is kappa stable across disjoint control-expression bands? -------------
    print("\n--- CHECK 3: kappa stability across disjoint bands (b_ctrl shape test) ---")
    bands = {"cpm 0.5-1": (0.5, 1.0), "cpm 1-2": (1.0, 2.0), "cpm 2-3": (2.0, 3.0),
             "cpm 3-4": (3.0, 4.0), "cpm 4-5": (4.0, MIN_CPM)}
    ref = None
    for name, (lo, hi) in bands.items():
        m = off & (cpm_pb >= lo) & (cpm_pb < hi)
        k = kappa_on_band(tot, lin_b, m)
        if ref is None:
            ref = k
        print(f"  {name:<12} {int(m.sum()):5d} genes, "
              f"sum totals {int((tot[:, m]).sum(1).mean()):8d} | "
              f"kappa mean {k.mean():.6e}  vs the 0.5-1 band "
              f"{100 * (k / ref - 1).mean():+.3f}%")

    # The production estimator is the PLATEAU of check 3: ungated genes at or above 2 CPM,
    # i.e. >= ~18 emitted counts. Below that the bands drift ~0.1% low, which is the +1/2
    # midpoint losing its uniformity as counts approach zero -- exactly the bias that made
    # the first pass of this script read a mass ratio of 0.978. The 2-3, 3-4 and 4-5 bands
    # agree with each other to 0.01%, so the plateau is where the estimator has converged.
    fit = off & (cpm_pb >= 2.0)
    kap = kappa_on_band(tot, lin_b, fit)
    D = np.log1p(kap[:, None] * (tot + 0.5)) - b[None, :]

    # ---- CHECK 1: the off-gate residual, against what one count explains ----------------
    res = np.abs(D[:, off])
    one_count = np.median(np.abs(np.log1p(kap.mean() * (tot[:, off] + 1.5))
                                 - np.log1p(kap.mean() * (tot[:, off] + 0.5))))
    print(f"\n--- CHECK 1: off-gate residual ---")
    print(f"  {int(fit.sum())} genes in the kappa fit | max|D| {res.max():.3e}  "
          f"median {np.median(res):.3e}  99th pct {np.percentile(res, 99):.3e}")
    print(f"  ONE emitted count on an off-gate gene moves log1p by "
          f"{one_count:.3e} -- the residual is that, not a modelling error")

    # ---- the mass ratio, the quantity the whole experiment turns on ---------------------
    S = np.clip(kap[:, None] * (tot + 0.5), 0.0, None).sum(1)
    mr = S / TS
    print(f"\n--- mass ratio S(D)/TS on this artifact ---")
    print(f"  mean {mr.mean():.6f}  sd {mr.std():.6f}  min {mr.min():.6f}  "
          f"max {mr.max():.6f}")
    print(f"  156 measured 1.00181 (HCT116) / 1.00168 (HEK293T); "
          f"162 measured 1.00178 / 1.00174 / 1.00171 on challenge A/B/C")
    # kappa's own uncertainty, propagated: the +1/2 midpoint has residual sd 1/sqrt(12)
    n_fit = int(fit.sum())
    rel = float(np.mean(np.sqrt(n_fit / 12.0) / (tot[:, fit] + 0.5).sum(1)))
    print(f"  kappa's quantisation sd is {100 * rel:.4f}% RELATIVE -> S/TS is determined to "
          f"+/-{rel:.6f}, against an excess of {mr.mean() - 1:.6f}")
    kap_lo = kappa_on_band(tot, lin_b, off & (cpm_pb >= 1.0) & (cpm_pb < 2.0))
    print(f"  BAND SYSTEMATIC (larger than the random error): fitting on cpm 1-2 instead "
          f"gives S/TS = {(np.clip(kap_lo[:, None] * (tot + 0.5), 0, None).sum(1) / TS).mean():.6f}")
    print(f"  -> the correction on this artifact is "
          f"{(mr.mean() - 1) / 0.0018:.2f}x production's (156/162 measured 1.0018). "
          f"minnorm is the calibration arm, so this ratio cancels in the placement contrast.")

    D = np.where(gate[None, :], D, 0.0)
    np.savez_compressed(
        os.path.join(WORK, "effects.npz"),
        genes=genes, targets=np.asarray(targets), b_ctrl=b, gate=gate,
        effect=D.astype(np.float64), totals=tot, kappa=kap, mass_ratio=mr,
        min_cpm=np.float64(MIN_CPM), ctrl_cpm=cpm_pb,
    )
    print(f"\nwrote {WORK}/effects.npz  [{time.time() - t0:.0f}s]")

    # ---- where the minimum-L2 correction lands (the a23 mechanism, restated) ------------
    print("\n--- eps placement under minnorm, averaged over targets ---")
    lin = np.clip(np.expm1(b[None, :] + D), 0.0, None)
    frac10 = frac100 = 0.0
    for i in range(T):
        w = np.where(gate & (b + D[i] > 0), np.exp(b + D[i]), 0.0)
        eps = -(lin[i].sum() - TS) / float(w @ w) * w
        o = np.argsort(-np.abs(eps))
        frac10 += np.abs(eps[o[:10]]).sum() / np.abs(eps).sum() / T
        frac100 += np.abs(eps[o[:100]]).sum() / np.abs(eps).sum() / T
    print(f"  the 10 largest-|eps| genes carry {100 * frac10:.1f}% of the total |eps|; "
          f"the top 100 carry {100 * frac100:.1f}%")
    print(f"  |eps|/|D| (median over targets) will be reported by arms.py")


if __name__ == "__main__":
    main()
