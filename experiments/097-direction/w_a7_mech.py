"""Where does the third donor's gain actually come from?

Adding a donor that carries 18 of 85 targets buys the same `sqrt(W)` as one that carries 84,
which is not what "more information per perturbation" predicts.  Three candidate mechanisms,
all separable:

  (i)  the rows the new donor CARRIES get a better mixture;
  (ii) the rows it does NOT carry move anyway, because `MultiDonorModel` re-derives `_renorm`,
       `_global` and `m_full` from the donor set, and `build()` re-normalises `E` by a MEDIAN
       row norm computed over all rows;
  (iii) rows no base donor carries switch from the `m_full` fallback to a real effect.

This splits `X` and `A` by row class and prints the per-class contribution before and after.
A global rescale control is included, because `W` is scale-invariant only in the linear limit.
"""
from __future__ import annotations

import sys

import numpy as np
from wkit import Ctx, uniform


def parts(c, dh, mask):
    w = c.w
    x = np.einsum("ij,ij->i", dh, c.d) * w
    a = np.einsum("ij,ij->i", dh, dh) * w
    return float(x[mask].sum()), float(a[mask].sum())


def main():
    name, extra = sys.argv[1], sys.argv[2]
    c2 = Ctx(name)
    b2 = c2.run(uniform)
    W2 = b2["W"]
    perts = list(c2.perts)
    base_cov = np.array([any(p in c2.model.pert_index[n] for n in c2.model.names)
                         for p in perts])
    print(f"########## {name}: base 2 donors -> +{extra}")
    print(f"  rows no BASE donor carries (m_full fallback): {int((~base_cov).sum())}")

    # global-scale control: W is only scale-invariant in the linear limit
    print("\n  global-scale control on the 2-donor base (pure rescale of E):")
    for s in (0.9, 0.95, 1.0, 1.05, 1.1):
        W = c2.XAW(c2.dh_exact(b2["E"] * s))[2]
        print(f"    E x {s:.2f} -> W {W:.6g}   mult {np.sqrt(W/W2):.4f}")
    E2 = b2["E"].copy()
    dh2 = c2.dh_exact(E2)
    del c2

    c3 = Ctx(name, extra_donors=(extra,))
    b3 = c3.run(uniform)
    assert list(c3.perts) == perts
    touch = np.array([extra in c3.model.pert_index and p in c3.model.pert_index[extra]
                      for p in perts])
    print(f"\n  rows the new donor carries: {int(touch.sum())}/{len(perts)}")
    print(f"  W 2-donor {W2:.6g}  ->  3-donor {b3['W']:.6g}   "
          f"mult {np.sqrt(b3['W']/W2):.4f}")

    # did the untouched rows' effects really move?
    dE = np.linalg.norm(b3["E"] - E2, axis=1) / np.maximum(np.linalg.norm(E2, axis=1), 1e-30)
    print(f"  relative ||dE|| on rows the new donor does NOT carry: "
          f"median {np.median(dE[~touch]):.4f}  max {dE[~touch].max():.4f}")
    print(f"  relative ||dE|| on rows it DOES carry:                "
          f"median {np.median(dE[touch]):.4f}  max {dE[touch].max():.4f}")

    dh3 = c3.dh_exact(b3["E"])
    print(f"\n  {'row class':<28}{'X2':>12}{'X3':>12}{'dX':>12}{'A2':>12}{'A3':>12}{'dA':>12}")
    for lbl, m in (("carried by new donor", touch), ("not carried", ~touch),
                   ("ALL", np.ones(len(perts), bool))):
        x2, a2 = parts(c3, dh2, m)
        x3, a3 = parts(c3, dh3, m)
        print(f"  {lbl:<28}{x2:>12.6g}{x3:>12.6g}{x3-x2:>+12.4g}"
              f"{a2:>12.6g}{a3:>12.6g}{a3-a2:>+12.4g}")

    # counterfactual: keep the 3-donor mixture ONLY on the rows the new donor carries
    Emix = E2.copy()
    Emix[touch] = b3["E"][touch]
    Wm = c3.XAW(c3.dh_exact(Emix))[2]
    print(f"\n  new mixture applied ONLY to carried rows: W {Wm:.6g}  "
          f"mult {np.sqrt(Wm/W2):.4f}")
    Eun = b3["E"].copy()
    Eun[touch] = E2[touch]
    Wu = c3.XAW(c3.dh_exact(Eun))[2]
    print(f"  new mixture applied ONLY to UNcarried rows: W {Wu:.6g}  "
          f"mult {np.sqrt(Wu/W2):.4f}")


if __name__ == "__main__":
    main()
