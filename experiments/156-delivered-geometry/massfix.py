"""156 — the mass constraint, and the cheapest way to satisfy it.

THE MECHANISM
-------------
The builder renders a log-space effect `h` as a profile
    p = expm1(b + h) / sum(expm1(b + h))
and the scorer reads back `log1p(TS*p) - b`. Since `sum(expm1(b)) == TS` exactly
(the control profile is itself a normalised log1p(TS*p) object), the read-back is
the identity **if and only if**
    S(h) := sum_g expm1(b_g + h_g) == TS.
That is one scalar constraint per row. The truth satisfies it by construction.
Our `h` does not, and the renderer therefore silently applies its own correction:
a division by S(h)/TS in the LINEAR space, which in log space is a shift
concentrated on the highest-expressed genes and -- critically -- nearly the same
direction for every row. Measured: 22% of the row norm at row-row cosine 0.303,
against an intended row-row cosine of 0.0040. Whitening's whole job is to make
those rows mutually dissimilar.

expm1 is convex and S(0) = TS, so S has no second root along the ray {gamma*h}:
scaling cannot fix this, and the correction must move off the ray. But it only
has to satisfy ONE scalar, so we may choose the cheapest direction. dS/dh_g =
exp(b_g + h_g), so the minimum-L2 correction is
    eps = -(S(h) - TS) * w / ||w||^2,     w_g = exp(b_g + h_g)
supported on the genes production already lets move (the CPM gate), iterated to
convergence. That is Newton's method with a minimum-norm step.

Everything is priced on both held Orion lines, in the scorer's own profile space.
"""
from __future__ import annotations

import json
import os

import numpy as np

from screen import (FROZEN, REPO, TS, WORK, _unit, endpoints, render,
                    restore_norms, t_centroid, t_rowgram)

GATE_CPM = 10.0          # a15's shipped config.min_cpm


def production_effect(surface: str):
    """The effect matrix exactly as production hands it to the renderer."""
    z = np.load(FROZEN % surface, allow_pickle=True)
    t = np.load(os.path.join(WORK, f"truth-{surface}.npz"), allow_pickle=True)
    D = z["baseline"].astype(np.float64)
    b = z["bctrl"].astype(np.float64)
    cpm = z["ctrl_cpm"].astype(np.float64)
    gate = (cpm > GATE_CPM) & (b > 0)                      # builder line 1207
    D = np.where(gate[None, :], D, 0.0)

    ptar = np.asarray([str(x) for x in z["targets"]])
    ttar = np.asarray([str(x) for x in t["targets"]])
    sel = np.asarray([{s: i for i, s in enumerate(ttar)}[s] for s in ptar])
    T = t["effect"].astype(np.float64)[sel]
    b_ref = t["b_ctrl"].astype(np.float64)

    genes = np.asarray([str(x) for x in z["genes"]])
    gpos = {g: i for i, g in enumerate(genes)}
    drop = {gpos[s] for s in ptar if s in gpos}
    mask = np.asarray([i for i in range(len(genes)) if i not in drop])
    return D, T, b, b_ref, gate, mask, ptar


def _builder():
    """Load the PRODUCTION builder module, so the screen cannot drift from the build.

    Experiment 156's first pass had its own copy of the correction, which used the UNCLIPPED
    mass while the renderer clips the profile at zero. That is precisely the class of
    screen/production mismatch that cost a19 -0.01352, so the reimplementation is deleted and
    the shipped function is called directly.
    """
    import importlib.util
    import sys
    path = os.path.join(REPO, "scripts/make_submission_pm_cd4.py")
    spec = importlib.util.spec_from_file_location("prod_builder_156", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_156"] = mod
    spec.loader.exec_module(mod)
    assert mod.TS == TS and mod.MASS_FIX is None
    return mod


_PROD = None


def mass_fix(D, b, gate, mode="minnorm"):
    """Row-wise wrapper around the production `mass_consistent_effect`."""
    global _PROD
    if _PROD is None:
        _PROD = _builder()
    out = np.empty_like(D)
    for i in range(D.shape[0]):
        out[i] = _PROD.mass_consistent_effect(D[i], b, gate, mode=mode)
    return out


def main() -> None:
    out = {}
    for surface in ("hct116", "hek293t"):
        D, T, b, b_ref, gate, mask, ptar = production_effect(surface)
        assert np.abs(b - b_ref).max() == 0.0, "builder and reference control differ"
        n = D.shape[0]
        print(f"\n{'=' * 104}\n{surface}: {n} targets, {gate.sum()} gated genes "
              f"(CPM > {GATE_CPM}), {len(mask)} scored genes")

        Pd = render(D, b)
        rel = np.median(np.linalg.norm(Pd - D, axis=1) / np.linalg.norm(D, axis=1))
        mr = np.expm1(b[None, :] + D).sum(1) / TS
        print(f"  production as shipped: mass ratio mean {mr.mean():.6f} "
              f"(sd {mr.std():.6f}) | render residual {rel * 100:.2f}% of row norm")

        arms = {}
        arms["delivered (production)"] = Pd
        arms["intended (INFEASIBLE)"] = D
        for mode in ("minnorm", "const"):
            h = mass_fix(D, b, gate, mode=mode)
            err = np.abs(render(h, b) - h).max()
            dn = np.median(np.linalg.norm(h - D, axis=1) / np.linalg.norm(D, axis=1))
            print(f"  massfix-{mode:8s}: render identity to {err:.2e} | "
                  f"correction {dn * 100:.3f}% of row norm | "
                  f"row norm x{np.median(np.linalg.norm(h, axis=1) / np.linalg.norm(D, axis=1)):.5f}")
            arms[f"massfix-{mode}"] = render(h, b)

        # geometry, on the delivered object, with the mass constraint re-imposed after
        hmf = mass_fix(D, b, gate, mode="minnorm")
        for label, fn, kw in (("centroid d0.75", t_centroid, {"dose": 0.75}),
                              ("centroid d1.0", t_centroid, {"dose": 1.0}),
                              ("rowgram lam0.01", t_rowgram, {"lam": 0.01})):
            g = np.where(gate[None, :], fn(hmf, **kw), 0.0)
            arms[f"massfix + {label}"] = render(mass_fix(g, b, gate), b)
            g2 = np.where(gate[None, :], fn(D, **kw), 0.0)
            arms[f"{label} only (pre-render)"] = render(g2, b)

        print(f"\n  {'arm':<30} {'pds':>9} {'d pds':>9} {'cos':>9} {'d cos %':>9} "
              f"{'row-row':>9} {'sse*':>9}")
        base = endpoints(Pd[:, mask], T[:, mask])
        rows = {}
        for name, P in arms.items():
            e = endpoints(P[:, mask], T[:, mask])
            rows[name] = e
            print(f"  {name:<30} {e['pds']:9.5f} {e['pds'] - base['pds']:+9.5f} "
                  f"{e['cos']:9.5f} {100 * (e['cos'] / base['cos'] - 1):+8.1f}% "
                  f"{e['rr']:9.5f} {e['sse']:9.6f}")
        out[surface] = rows

    with open(os.path.join(WORK, "massfix.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"\n{'=' * 104}\nWORST-SURFACE (selection rule) — delta pds vs production")
    names = list(out["hct116"])
    base = {s: out[s]["delivered (production)"]["pds"] for s in out}
    tab = sorted(((min(out[s][k]["pds"] - base[s] for s in out), k) for k in names),
                 reverse=True)
    print(f"  {'arm':<30} {'worst':>10} {'hct116':>10} {'hek293t':>10} "
          f"{'cos hct':>9} {'cos hek':>9}")
    for worst, k in tab:
        print(f"  {k:<30} {worst:+10.5f} "
              f"{out['hct116'][k]['pds'] - base['hct116']:+10.5f} "
              f"{out['hek293t'][k]['pds'] - base['hek293t']:+10.5f} "
              f"{out['hct116'][k]['cos']:9.5f} {out['hek293t'][k]['cos']:9.5f}")


if __name__ == "__main__":
    main()
