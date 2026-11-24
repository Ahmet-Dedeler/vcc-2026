"""Did a12's tiered emission move the PROFILE, or only the emitted cells?

`A`, `B` and `X` are functionals of the group SUMS. The 096 model assumes emission is invariant
for them (the quasi emitter repairs each gene's panel total exactly), so a12's profile should be
exactly a10's re-rendered at `s = 1.0108/1.0796`. If the `tail-spread 1.0` policy broke the
repair, `X` would move without `A` giving it away -- the direction can rotate while the norm is
preserved, and RESEARCH-METHOD's "validation statistic at the wrong scale" note says a cosine on
the whole vector cannot see it.

Two things measured here, and the second is the point:
  * per-row cosine and norm ratio of the shipped profile against the re-rendered one -- a
    NECESSARY condition only;
  * the residual's size relative to `X` itself, which is what decides whether the cosine check
    has any power at all.

    uv run python mech_a12_profile.py
"""
from __future__ import annotations

import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from mech_alpha_chal import CTX, TS, load, render, terms  # noqa: E402

W = f"{ROOT}/experiments/096-mse-floor/work"
ALPHA_PROD = 1.0796


def compare(base_stats, tgt_stats, s, label):
    Db = load(base_stats)
    Dt = load(tgt_stats)
    print(f"\n===== {label}: base re-rendered at s={s:.6f} vs the shipped target =====")
    print(f"  {'ctx':>4}{'A(render)':>12}{'A(ship)':>12}{'ratio':>8}{'B(render)':>12}"
          f"{'B(ship)':>12}{'ratio':>8}{'cos min':>10}{'cos med':>10}"
          f"{'|dres|/|dh|':>13}")
    agg = []
    for c in CTX:
        db, dt = Db[c], Dt[c]
        sums = render(db["bctrl"], db["dh"], db["S"], s)
        Ar, Br, _ = terms(db, sums)
        # the re-rendered effect, in the same coordinates the metric uses
        from cell_eval2.prep import bulk_lognorm_means
        dh_r = bulk_lognorm_means(np.asarray(sums, dtype=np.float64), TS) - db["bctrl"][None, :]
        dh_s = dt["dh"]
        for k, col in db["cols"].items():
            dh_r[k, col] = 0.0
        for k, col in dt["cols"].items():
            dh_s[k, col] = 0.0
        nr = np.linalg.norm(dh_r, axis=1)
        ns = np.linalg.norm(dh_s, axis=1)
        cos = np.einsum("ij,ij->i", dh_r, dh_s) / np.maximum(nr * ns, 1e-300)
        res = np.linalg.norm(dh_s - dh_r, axis=1) / np.maximum(nr, 1e-300)
        At, Bt, _ = terms(dt, dt["sums"] if "sums" in dt else None) if False else (None, None, None)
        # A/B of the shipped target come straight from its own means
        a = np.einsum("ij,ij->i", dh_s, dh_s)
        from cell_eval2.metrics.delta import _across_pert_budget, _row_weights
        w = _row_weights(len(dt["tgt"]), dt["G"], dt["cols"])
        As = float((a * w).sum())
        Bs = _across_pert_budget(dt["bctrl"][None, :] + dh_s, dt["cols"], w)
        print(f"  {c:>4}{Ar:>12.6f}{As:>12.6f}{As/Ar:>8.4f}{Br:>12.6f}{Bs:>12.6f}"
              f"{Bs/Br:>8.4f}{cos.min():>10.6f}{np.median(cos):>10.6f}{res.mean():>13.6f}")
        agg.append((cos, res, ns / np.maximum(nr, 1e-300)))
    cos = np.concatenate([a[0] for a in agg])
    res = np.concatenate([a[1] for a in agg])
    rat = np.concatenate([a[2] for a in agg])
    print(f"  all rows: cos min {cos.min():.6f}  p5 {np.percentile(cos,5):.6f}  "
          f"median {np.median(cos):.6f}   ||d_res||/||dh|| mean {res.mean():.5f}  "
          f"max {res.max():.5f}   norm ratio median {np.median(rat):.5f}")
    return cos, res


def main():
    # a11 = a10's profile at VCC_ALPHA 0.838 (uniform emission) -- the CONTROL for this test
    compare(f"{W}/a10-stats.npz", f"{W}/a11-stats.npz", 0.838 / ALPHA_PROD,
            "a10 -> a11  (amplitude only, emission unchanged)")
    # a12 = the same profile at 1.0108 but with the tiered spread policy -- the TEST
    compare(f"{W}/a10-stats.npz", f"{W}/a12-stats.npz", 1.0108 / ALPHA_PROD,
            "a10 -> a12  (amplitude + tiered spread/overdispersion)")

    print("\n===== how much X a given angular error can move =====")
    print("  X = sum_p w_p <dh_p, d_p>.  A residual `r` with ||r||/||dh|| = e can move X by at")
    print("  most e*sqrt(A*sumD) if it is perfectly aligned with the truth, i.e.")
    A, sD, X = 0.248147, 0.2319, 0.018753
    for e in (0.0001, 0.001, 0.01, 0.02):
        print(f"    ||r||/||dh|| = {e:7.4f}  ->  |dX| <= {e*np.sqrt(A*sD):.6f}  "
              f"= {100*e*np.sqrt(A*sD)/X:6.1f}% of X")
    print("  So a cosine of 0.9999 on the whole vector still permits a 12% swing in X.")


if __name__ == "__main__":
    main()
