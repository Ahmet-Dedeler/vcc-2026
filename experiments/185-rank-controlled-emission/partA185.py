"""185 Part A — controllability: the reachable Wilcoxon-z interval per gene.

For every h1-test perturbation, take a15's SHIPPED per-gene totals M_g and the real control
cells, and ask: over legal allocations of M_g counts to 400 cells at fixed per-cell total
T = S/400 (level cap T per cell), what z can gene g reach?

Exact objects per gene (`eme185.envelope_all`):
  z_even   the maximally-spread allocation (floor/ceil of M/400 in every cell)
  z_conc   the most concentrated allocation (cells filled to T, the rest at zero)
  z_min    the minimum-Fbar member of the zero-inflated family (lower convex envelope)
  z_max    the maximum-Fbar member of the hull family (upper concave envelope)
  z_floor  the member closest to the null (the reachable |z| floor)

Reported by decile of control mean CPM, and for one headline perturbation in full.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness185 as H  # noqa: E402
from eme185 import EME  # noqa: E402

n = 400


def main() -> None:
    t0 = time.time()
    C = H.load_controls()
    cs, tie, keep, cpm = C["cs"], C["tie"], C["keep"], C["cpm"]
    G, m = cs.shape
    TOT = H.load_totals()
    tmask = TOT["targets"] != "non-targeting"
    targets = TOT["targets"][tmask]
    totals = TOT["totals"][tmask]
    Hd = H.load_H()
    assert np.array_equal(Hd["genes"], C["genes"])
    hpos = {t: i for i, t in enumerate(Hd["targets"])}
    eme = EME(cs, tie, 20000, n)
    zero_frac = (cs == 0).mean(1)                      # control zero fraction per gated gene
    dec = np.digitize(cpm[keep], np.quantile(cpm[keep], np.linspace(0, 1, 11)[1:-1]))

    Z = np.zeros((targets.size, G, 6), np.float32)
    absH = np.zeros((targets.size, G), np.float32)
    for i, t in enumerate(targets):
        Mfull = totals[i]
        T = int(Mfull.sum()) // n
        Z[i] = eme.envelope(Mfull[keep], T)
        absH[i] = np.abs(Hd["H"][hpos[t]][keep])
        if (i + 1) % 20 == 0:
            print(f"  {i+1}/{targets.size} [{time.time()-t0:.0f}s]", flush=True)
    np.savez_compressed(os.path.join(H.WORK, "partA.npz"), targets=targets, genes=C["genes"][keep],
                        Z=Z, absH=absH, zero_frac=zero_frac, cpm=cpm[keep], decile=dec)

    z_even, z_conc, z_min, z_max, z_floor = [Z[:, :, k].astype(np.float64) for k in range(5)]
    hi = np.maximum(np.abs(z_min), np.abs(z_max))          # reachable |z| ceiling
    lo = np.abs(z_floor)                                   # reachable |z| floor
    width = z_max - z_min
    null_in = (z_min < 0) & (z_max > 0)

    def row(mask, label):
        return (f"{label:<12}{np.median(width[mask]):>9.1f}{np.median(z_min[mask]):>9.1f}"
                f"{np.median(z_even[mask]):>9.2f}{np.median(z_max[mask]):>9.2f}"
                f"{np.median(hi[mask]):>9.1f}{np.median(lo[mask]):>9.3f}"
                f"{100*null_in[mask].mean():>8.1f}%{100*(lo[mask] < 1.0).mean():>8.1f}%"
                f"{100*(hi[mask] > 5.0).mean():>8.1f}%{np.median(np.broadcast_to(zero_frac, mask.shape)[mask]):>8.2f}")

    hdr = (f"{'decile':<12}{'width':>9}{'z_min':>9}{'z_even':>9}{'z_max':>9}{'max|z|':>9}"
           f"{'min|z|':>9}{'null in':>9}{'floor<1':>9}{'ceil>5':>9}{'ctrl0':>8}")
    print("\nPART A — reachable z interval per gene, ALL 85 perturbations pooled "
          "(medians over (perturbation, gene) pairs; decile 1 = lowest control CPM)")
    print(hdr)
    for d in range(10):
        mask = np.broadcast_to(dec == d, z_min.shape)
        lab = f"d{d+1} {np.quantile(cpm[keep][dec == d], 0.5):.0f}cpm"
        print(row(mask, lab))
    print(row(np.ones_like(null_in), "ALL"))

    # the pool that reach walks: genes with large |H| vs small |H|
    q = np.quantile(absH, [0.5, 0.9, 0.99], axis=1)
    print("\nby |H| stratum (per perturbation quantiles of |H| over the gated axis):")
    print(hdr)
    for lab, mask in (("|H| top1%", absH >= q[2][:, None]), ("|H| top10%", absH >= q[1][:, None]),
                      ("|H| below med", absH < q[0][:, None])):
        print(row(mask, lab))

    # headline perturbation: the one with median n_real is not known here; use the first
    i = 0
    t = targets[i]
    print(f"\nHEADLINE — {t}: S = {totals[i].sum()}, T = {totals[i].sum()//n}, gated {G}")
    print(hdr)
    for d in range(10):
        mask = np.zeros_like(null_in)
        mask[i] = dec == d
        print(row(mask, f"d{d+1}"))
    m_all = np.zeros_like(null_in)
    m_all[i] = True
    print(row(m_all, "ALL"))
    # a15's own realised |z| on this target, for scale: from 171's table if present
    p171 = f"{H.ROOT}/experiments/171-reach-prefix/work/de171.a15.npz"
    if os.path.exists(p171):
        from scipy.stats import norm
        z171 = np.load(p171, allow_pickle=True)
        j = list(np.asarray(z171["perts"]).astype(str)).index(t)
        pa = z171["pred_p"][j]
        za = norm.isf(pa / 2)
        print(f"\na15's realised |z| on {t} (171's table): median {np.median(za):.2f}, "
              f"90th pct {np.quantile(za, 0.9):.2f}, max {za.max():.1f}; "
              f"fraction of genes whose a15 |z| lies inside [min|z|, max|z|] reachable here: "
              f"{np.mean((za >= lo[i]) & (za <= hi[i])):.3f}")
    print(f"\ndone [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
