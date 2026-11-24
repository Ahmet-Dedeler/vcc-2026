"""Is `W` really amplitude-invariant through this pipeline?

The whole objective rests on `W = X^2/A` being the only thing an arm can move, because the
amplitude is re-optimised onto the rho cliff afterwards.  That argument is exact when `dh` is
LINEAR in the effect (`X -> sX`, `A -> s^2 A`).  It is not: `dh = log1p(TS*floor(n*L*p)) -
b_ctrl` with `p` normalised and CPM-gated, which compresses low-expression genes.

So before any arm can be read, measure `W(s)` for a pure rescale of `E`.  Whatever `W(s)`
varies by across the plausible amplitude range is the part of an arm's `W` change that an
amplitude move could have produced on its own, and it must be subtracted from the claim.
"""
from __future__ import annotations

import sys

import numpy as np
from wkit import Ctx, uniform

S = (0.25, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0)


def main():
    for name in sys.argv[1:] or ["h1-test", "orion"]:
        c = Ctx(name)
        b = c.run(uniform)
        Wb = b["W"]
        print(f"\n########## {name}   W(1.0) = {Wb:.8g}   "
              f"median ||E||={b['Enorm']:.4f}  median ||dh||={b['dhnorm']:.5f}")
        print(f"  {'s':>6}{'X':>13}{'A':>13}{'W':>13}{'sqrt(W/Wb)':>13}"
              f"{'X/s':>13}{'A/s^2':>13}{'med||dh||':>12}")
        for s in S:
            dh = c.dh_exact(b["E"] * s)
            X, A, W = c.XAW(dh)
            print(f"  {s:>6.2f}{X:>13.6g}{A:>13.6g}{W:>13.6g}{np.sqrt(W/Wb):>13.4f}"
                  f"{X/s:>13.6g}{A/s**2:>13.6g}"
                  f"{float(np.median(np.linalg.norm(dh, axis=1))):>12.5f}")
        del c


if __name__ == "__main__":
    main()
