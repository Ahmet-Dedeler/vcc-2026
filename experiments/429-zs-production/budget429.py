"""429 diagnostic: why zs keeps only part of its gain on the MDE carriers. On the dumped seed-0 blocks (the same 6 targets
with base and MDE carriers) compare (a) the carriers' sub-gate holdings (the zs pass's library-compensation budget),
(b) the achieved A-channel move with the real budget vs a counterfactual unlimited budget (the vendored loop with room =
BIG; not a valid block, only the achieved |dA|), and (c) the delivered |lfc| of the genes zs can move.
usage: uv run python budget429.py"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for _p in (str(ROOT / "experiments/323-current-family-validation"), str(ROOT / "src"), str(HERE)):
    sys.path.insert(0, _p)
import evaluate323 as E  # noqa: E402
import src_rank_mag as RM  # noqa: E402


def achieved(X, n_car, cc, l, unlimited):
    """The vendored zs loop's per-gene achieved dA (no block written back)."""
    n = X.shape[0]
    keep = np.asarray(cc.keep, bool)
    gidx, sidx = np.flatnonzero(keep), np.flatnonzero(~keep)
    s = X.sum(1) / 1e6
    car = np.arange(n_car)
    Xk = X[:, gidx]
    A = (Xk / s[:, None]).mean(0)
    ref = np.asarray(cc.ref_lin, float)
    ok = np.isfinite(l) & (A > 0) & (l != 0)
    touch = RM.select("dn", l, ok)
    A_t = (ref + 1e-9) * np.power(2.0, 2.5 * l) - 1e-9
    dA = np.where(touch, A_t - A, 0.0)
    order = np.flatnonzero(touch)[np.argsort(-np.abs(l[touch]), kind="stable")]
    sub_avail = X[np.ix_(car, sidx)].sum(1).astype(np.int64)
    sub_delta = np.zeros(n_car, dtype=np.int64)
    cw = 1.0 / (n * s[car])
    groups, gw = RM._levels(cw)
    ach = np.zeros(gidx.size)
    for j in order:
        x0 = Xk[car, j].copy()
        room0 = np.full(n_car, RM.BIG) if unlimited else np.maximum(sub_avail + sub_delta, 0)
        xl, xh, *_ = RM.cell_bounds(x0, s[car], cc.cs[j])
        x = RM._zs_step(x0, xl, xh, room0, cw, dA[j], groups, gw)
        dx = x - x0
        if not dx.any():
            continue
        A_new = A[j] + (dx * cw).sum()
        if np.sign(np.log2((A_new + 1e-9) / (ref[j] + 1e-9))) != np.sign(l[j]):
            continue
        sub_delta -= dx
        ach[j] = (dx * cw).sum()
    return ach, dA, touch, sub_avail


def main():
    r, cc, _ = E.load_reference()
    del r
    rows = []
    for p in sorted((HERE / "work/blocks").glob("base_s0_bps_*.npz")):
        z = np.load(p)
        X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
        n_car = int(z["n_car"])
        l = RM.block_lfc(X, cc)
        a_real, dA, touch, sub = achieved(X, n_car, cc, l, False)
        a_inf, _, _, _ = achieved(X, n_car, cc, l, True)
        moved = a_real != 0
        arm = "MDE" if "mde" in p.stem else "base"
        rows.append((arm, str(z["target"]), np.median(sub), float(np.mean(sub < 50)), int(moved.sum()),
                     float(np.abs(a_real).sum()), float(np.abs(a_inf).sum()),
                     float(np.median(np.abs(l[moved]))) if moved.any() else np.nan))
    print("arm   target   carrier sub-gate median  frac<50  genes moved  sum|dA| real  sum|dA| unlimited  median |l| moved")
    for rw in rows:
        print(f"{rw[0]:5s} {rw[1]:8s} {rw[2]:10.0f} {rw[3]:16.3f} {rw[4]:11d} {rw[5]:14.2f} {rw[6]:18.2f} {rw[7]:16.3f}")
    for arm in ("base", "MDE"):
        R = [rw for rw in rows if rw[0] == arm]
        print(f"{arm}: sum|dA| real {sum(r_[5] for r_ in R):.1f}, unlimited budget {sum(r_[6] for r_ in R):.1f} "
              f"(budget binds on {1 - sum(r_[5] for r_ in R) / sum(r_[6] for r_ in R):.1%}); genes moved {sum(r_[4] for r_ in R)}")


if __name__ == "__main__":
    main()
