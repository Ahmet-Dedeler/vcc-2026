"""Is the residual between surrogate and realized DETERMINISTIC (the floor) or stochastic?

`emit_counts(mode="quasi")` sets `totals = floor(n * L * profile)` and then repairs every gene's
total exactly, so the scored group sums are that floor -- no emission noise reaches the
pseudobulk at all.  The only stochastic input is `L = int(median(libs[pick]))`, one integer per
perturbation.  If that is the whole story, then

    dh_floor = log1p(TS * floor(n*L*p) / sum_g floor(n*L*p)) - b_ctrl

reproduces the realized `dh` bit-for-bit, and the "surrogate" is not a surrogate at all.
"""
from __future__ import annotations

import numpy as np
from wkit import MIN_CPM, ROOT, TS, Ctx, uniform


def recover_L(total: int, prof: np.ndarray, n: int = 400) -> int:
    """The integer `L` the emitter used, from the perturbation's observed panel total."""
    lo = int(total / n)
    for step in range(0, 3000):
        for L in (lo + step, lo - step):
            if L > 0 and int(np.floor(n * L * prof).sum()) == total:
                return L
    return -1


def main():
    for name in ("h1-test", "orion"):
        c = Ctx(name)
        E = c.effects(uniform)
        obs = c.t["pm"][c.t["keep"]]                     # scorer's own log1p profile per row
        # observed group totals, from the cached prediction's group sums
        from msekit import group_stats
        pp, ps, pc, pj, pg = group_stats(str(__import__("wkit").HARNESS[name]["pred"]))
        tot = {str(p): int(ps[i].sum()) for i, p in enumerate(pp)}

        dh_floor = np.empty_like(c.dh_real)
        Ls, exact = [], 0
        for i, p in enumerate(c.perts):
            prof = c.model.profile_from_effect(E[i], c.b_ctrl, min_cpm=MIN_CPM,
                                               gate_unexpressed=True)
            L = recover_L(tot[p], prof)
            Ls.append(L)
            t = np.floor(400 * L * prof) if L > 0 else 400 * 1 * prof
            dh_floor[i] = np.log1p(TS * t / t.sum()) - c.b_ctrl
            exact += int(L > 0)
        for k, col in c.t["cols"].items():
            dh_floor[k, col] = 0.0

        r = np.abs(dh_floor - c.dh_real)
        Xf, Af, Wf = c.XAW(dh_floor)
        Xr, Ar, Wr = c.XAW(c.dh_real)
        print(f"\n########## {name} ##########")
        print(f"  L recovered exactly for {exact}/{len(c.perts)} perturbations; "
              f"L range {min(Ls)}..{max(Ls)}")
        print(f"  max |dh_floor - dh_realized| = {r.max():.3e}   "
              f"mean = {r.mean():.3e}")
        print(f"  X  floor {Xf:.8g}  realized {Xr:.8g}  ratio {Xf/Xr:.8f}")
        print(f"  A  floor {Af:.8g}  realized {Ar:.8g}  ratio {Af/Ar:.8f}")
        print(f"  W  floor {Wf:.8g}  realized {Wr:.8g}  sqrt-ratio {np.sqrt(Wf/Wr):.8f}")
        np.savez(ROOT / f"experiments/097-direction/work/{name}-floorL.npz",
                 L=np.asarray(Ls), Xf=Xf, Af=Af, Wf=Wf, Xr=Xr, Ar=Ar, Wr=Wr)
        del c


if __name__ == "__main__":
    main()
