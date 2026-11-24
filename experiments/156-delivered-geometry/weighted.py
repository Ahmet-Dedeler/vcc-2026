"""156b — WHERE to put the mass correction. One scalar constraint, a family of directions.

The constraint is `sum_g clip(expm1(b_g + h_g), 0) = TS`, one scalar per row. Any `eps` with
`<eps, w> = -(S - TS)`, `w_g = exp(b_g + h_g)`, satisfies it. Minimising a weighted norm
`sum_g eps_g^2 / u_g` subject to that gives

    eps_g = lambda * u_g * w_g,      lambda = -(S - TS) / sum_g u_g * w_g^2

so the whole family is indexed by the weight `u`. Two members are already measured:

    u = 1        -> eps ∝ w        minimum L2      (`minnorm`, SHIPPED in a23)
    u = 1/w      -> eps ∝ 1        uniform shift   (`const`, measured -0.0056/-0.0031 pds)

`const` is far worse than `minnorm`, and the two differ only in how sharply the correction
concentrates on high-expression genes -- `minnorm` concentrates MORE. That gradient says to
test concentrating harder still (`u = w`, `u = w^2`), and it also raises the opposite question:
min-L2 dumps the correction on the highest-expressed genes, which are exactly the genes whose
reference effect is measured most reliably and which carry the most cosine. A confidence-aware
weight (`u ∝ 1/(|effect| + c)`) instead steers the correction toward genes we have no signal on.

Both directions are cheap to test and only one can be right, so this is a 1-D screen, not a
guess. Run AFTER a23 (which ships `minnorm`) so any winner is a clean single-variable successor.
"""
from __future__ import annotations

import json
import os

import numpy as np

from massfix import GATE_CPM, production_effect  # noqa: F401
from screen import TS, WORK, endpoints, render

TOL = 1e-12


def weighted_mass_fix(D, b, gate, weight="one", c=0.01, iters=64):
    """Row-wise weighted-minimum-norm mass correction. `weight` selects `u`."""
    h = D.astype(np.float64).copy()
    for _ in range(iters):
        live = gate[None, :] & (b[None, :] + h > 0.0)
        S = np.clip(np.expm1(b[None, :] + h), 0.0, None).sum(1)
        r = S - TS
        if np.abs(r).max() <= TOL * TS:
            break
        w = np.where(live, np.exp(b[None, :] + h), 0.0)
        if weight == "one":
            u = np.ones_like(w)
        elif weight == "inv_w":
            u = np.where(w > 0, 1.0 / np.maximum(w, 1e-300), 0.0)
        elif weight == "w":
            u = w
        elif weight == "w2":
            u = w * w
        elif weight == "inv_effect":
            u = 1.0 / (np.abs(D) + c)
        elif weight == "effect":
            u = np.abs(D) + c
        else:
            raise ValueError(weight)
        u = np.where(live, u, 0.0)
        denom = (u * w * w).sum(1)
        lam = np.where(denom > 0, -r / np.maximum(denom, 1e-300), 0.0)
        h = h + lam[:, None] * u * w
    return h


ARMS = [("minnorm  (u=1, SHIPPED)", "one", 0.0),
        ("const    (u=1/w)", "inv_w", 0.0),
        ("sharper  (u=w)", "w", 0.0),
        ("sharpest (u=w^2)", "w2", 0.0),
        ("low-effect (u=1/(|e|+0.001))", "inv_effect", 0.001),
        ("low-effect (u=1/(|e|+0.01))", "inv_effect", 0.01),
        ("low-effect (u=1/(|e|+0.1))", "inv_effect", 0.1),
        ("high-effect (u=|e|+0.01)", "effect", 0.01)]


def main() -> None:
    out = {}
    for surface in ("hct116", "hek293t"):
        D, T, b, b_ref, gate, mask, ptar = production_effect(surface)
        base = endpoints(render(D, b)[:, mask], T[:, mask])
        print(f"\n{'=' * 96}\n{surface}")
        print(f"  {'arm':<30} {'pds':>9} {'d pds':>9} {'cos':>9} {'d cos %':>9} "
              f"{'sse*':>10} {'|eps|/|D|':>10}")
        print(f"  {'delivered (production)':<30} {base['pds']:9.5f} {0.0:+9.5f} "
              f"{base['cos']:9.5f} {0.0:+8.1f}% {base['sse']:10.6f} {'-':>10}")
        rows = {}
        for name, wkey, c in ARMS:
            h = weighted_mass_fix(D, b, gate, weight=wkey, c=c)
            e = endpoints(render(h, b)[:, mask], T[:, mask])
            rel = float(np.median(np.linalg.norm(h - D, axis=1)
                                 / np.linalg.norm(D, axis=1)))
            rows[name] = dict(**e, rel=rel)
            print(f"  {name:<30} {e['pds']:9.5f} {e['pds'] - base['pds']:+9.5f} "
                  f"{e['cos']:9.5f} {100 * (e['cos'] / base['cos'] - 1):+8.1f}% "
                  f"{e['sse']:10.6f} {rel * 100:9.3f}%")
        rows["delivered (production)"] = dict(**base, rel=0.0)
        out[surface] = rows

    with open(os.path.join(WORK, "weighted.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\n{'=' * 96}\nWORST-SURFACE, vs production")
    b0 = {s: out[s]["delivered (production)"] for s in out}
    names = [n for n, _, _ in ARMS]
    tab = sorted(((min(out[s][n]["pds"] - b0[s]["pds"] for s in out), n) for n in names),
                 reverse=True)
    print(f"  {'arm':<30} {'worst d pds':>12} {'d cos% hct':>11} {'d cos% hek':>11} "
          f"{'sse* hct':>10} {'sse* hek':>10}")
    for worst, n in tab:
        print(f"  {n:<30} {worst:+12.5f} "
              f"{100 * (out['hct116'][n]['cos'] / b0['hct116']['cos'] - 1):+10.1f}% "
              f"{100 * (out['hek293t'][n]['cos'] / b0['hek293t']['cos'] - 1):+10.1f}% "
              f"{out['hct116'][n]['sse']:10.6f} {out['hek293t'][n]['sse']:10.6f}")


if __name__ == "__main__":
    main()
