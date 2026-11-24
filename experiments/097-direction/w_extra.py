"""STEP 3 -- the honesty pass: noise floor, fine grids, cross-surface fits, and the
direction/magnitude decomposition of the one arm that looked positive on both surfaces.

Three questions the arm sweep leaves open:

  1. What is the NOISE FLOOR on `sqrt(W/W_base)`?  The reconstruction is deterministic given
     the emitter's per-perturbation `L`, so the floor is exactly the spread of `W` over
     alternative `L` draws.  Anything smaller than that is unreadable.
  2. Where is the argmax of each tuned family, and what does the parameter fitted on ONE
     harness buy on the OTHER?  An in-sample argmax is not a result.
  3. `a6 normrel^gam` reweights donors by the INVERSE of their effect norm, which changes two
     things at once: the mixture DIRECTION of each row, and the row's own MAGNITUDE relative
     to the other rows (a per-perturbation shrinkage, which 096 (c) measured as nil out of
     sample).  Splitting them says which mechanism is actually paying.
"""
from __future__ import annotations

import json
import sys

import numpy as np
from w_arms import key_normrel, key_selfkd, soft
from wkit import ROOT, Ctx, fixed, uniform

GAMS = (-3.0, -2.5, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.5, 1.0)
SHARES = tuple(np.round(np.arange(0.05, 1.0, 0.05), 2))


def main():
    name = sys.argv[1]
    c = Ctx(name)
    out = {"harness": name}
    base = c.run(uniform)
    Wb = base["W"]
    nb = np.linalg.norm(base["E"], axis=1)
    print(f"########## {name}   W_base = {Wb:.8g}")

    # -- 1. noise floor: re-draw the emitter's L ------------------------------------------
    print("\n--- noise floor on sqrt(W/W_base) from the emitter's only stochastic input, L")
    mults = []
    for seed in range(8):
        Lr = c.resample_L(seed)
        X, A, W = c.XAW(c.dh_exact(base["E"], Lr))
        mults.append(np.sqrt(W / Wb))
    mults = np.array(mults)
    print(f"  L-redraw multipliers: {np.round(mults, 5).tolist()}")
    print(f"  mean {mults.mean():.5f}  sd {mults.std(ddof=1):.5f}  "
          f"range {mults.min():.5f}..{mults.max():.5f}")
    print(f"  => anything inside +/-{2*mults.std(ddof=1):.4f} of 1.000 is emitter noise")
    out["noise_floor_sd"] = float(mults.std(ddof=1))
    out["noise_floor_mults"] = mults.tolist()

    # -- 2. fine a3 grid ------------------------------------------------------------------
    print("\n--- (a3 fine) global k562 share")
    a3 = {}
    for s in SHARES:
        X, A, W = c.XAW(c.dh_exact(c.effects(
            fixed({"k562_gwps": float(s), "orion_hek293t": 1.0 - float(s)}))))
        a3[float(s)] = dict(X=X, A=A, W=W, mult=float(np.sqrt(W / Wb)))
        print(f"  share {s:.2f}   X={X:>11.6g} A={A:>11.6g} W={W:>11.6g} "
              f"mult={np.sqrt(W/Wb):.4f}")
    best = max(a3, key=lambda k: a3[k]["mult"])
    print(f"  in-sample argmax share = {best:.2f} at {a3[best]['mult']:.4f}")
    out["a3"] = a3

    # -- 3. fine a6 normrel grid, plus the direction / magnitude split ---------------------
    for kname, keys in (("normrel", key_normrel(c)), ("selfkd", key_selfkd(c))):
        print(f"\n--- (a6 fine) {kname}: full / direction-only / magnitude-only")
        rec = {}
        print(f"  {'gam':>7}{'W':>13}{'mult':>9}{'dir-only':>11}{'mag-only':>11}"
              f"{'med rownorm ratio':>20}")
        for g in GAMS:
            E = c.effects(soft(c, keys, g))
            n = np.linalg.norm(E, axis=1)
            Edir = E * (nb / np.maximum(n, 1e-30))[:, None]      # arm's mixture, base's norms
            Emag = base["E"] * (n / np.maximum(nb, 1e-30))[:, None]   # base mixture, arm norms
            W = c.XAW(c.dh_exact(E))[2]
            Wd = c.XAW(c.dh_exact(Edir))[2]
            Wm = c.XAW(c.dh_exact(Emag))[2]
            rec[g] = dict(W=W, mult=float(np.sqrt(W / Wb)),
                          mult_dir=float(np.sqrt(Wd / Wb)), mult_mag=float(np.sqrt(Wm / Wb)))
            print(f"  {g:>7.2f}{W:>13.6g}{np.sqrt(W/Wb):>9.4f}{np.sqrt(Wd/Wb):>11.4f}"
                  f"{np.sqrt(Wm/Wb):>11.4f}{float(np.median(n/nb)):>20.4f}")
        out[f"a6_{kname}"] = {str(k): v for k, v in rec.items()}

    # -- 4. the pure per-perturbation magnitude controls -----------------------------------
    print("\n--- per-row magnitude controls (base mixture direction, re-scaled rows)")
    med = float(np.median(nb))
    for lbl, n_new in (("equal row norms", np.full_like(nb, med)),
                       ("rows ^0.5 of base", med * (nb / med) ** 0.5),
                       ("rows ^1.5 of base", med * (nb / med) ** 1.5),
                       ("rows ^2.0 of base", med * (nb / med) ** 2.0)):
        E = base["E"] * (n_new / np.maximum(nb, 1e-30))[:, None]
        W = c.XAW(c.dh_exact(E))[2]
        print(f"  {lbl:<22} W={W:>11.6g}  mult={np.sqrt(W/Wb):.4f}")
        out.setdefault("magnitude_controls", {})[lbl] = float(np.sqrt(W / Wb))

    (ROOT / f"experiments/097-direction/work/{name}-extra.json").write_text(
        json.dumps(out, indent=2))
    print(f"\nwrote work/{name}-extra.json")


if __name__ == "__main__":
    main()
