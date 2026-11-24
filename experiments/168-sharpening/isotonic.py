"""168c — replace the fitted exponent with the MEASURED conditional mean, fitted cross-surface.

## Why

168b established that `E[true | predicted]` is strongly convex in log-log: the top 5% of `|D|`
carries ~70% of the aggregate alignment `X = sum(D*T)`, the middle bins carry 0-3% each and are
noise, and a global power-law fit reads q ~ 1 while the top ventile's local slope is ~3. Global
sharpening (`|D|^q`) and head multiplication (`top f% x g`) are both *parametric guesses* at that
curve, and their parameters are fitted -- which is exactly what will not transfer to the unseen
D/E/F panel.

But the curve itself is measurable. For an inner-product/cosine objective the optimal transform
of a noisy predictor is its conditional mean, so the right transform is

    |D| -> m(|D|) := E[ t * sign(D) | |D| ]

estimated non-parametrically on equal-count bins of `|D|`, made monotone (isotonic), and applied
as a per-coordinate multiplier `w = m(|D|)/|D|`. One interpretable temperature `tau` then spans
the family:

    w^tau     tau = 0 -> identity (what we ship)
              tau = 1 -> the conditional mean (optimal for squared error / cosine)
              tau > 1 -> sharper than the conditional mean, which is what `pds`, being a RANK
                         metric rather than an estimator, empirically prefers

## The protocol that makes this honest

The curve must NEVER be estimated on the surface it is scored on -- that would be fitting to the
held truth, and it would look spectacular. So the curve is fitted on ONE Orion line and applied
to the OTHER, both directions reported. That is also exactly how it would deploy: fit on the
donor lines (leave-one-donor-out), apply to A/B/C, whose truth we never see.

`tau` is the only free number, and a cross-surface sweep is a legitimate way to choose it
because the surface it is chosen on is not the surface it is read on.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, os.path.join(REPO, "experiments/156-delivered-geometry"))
sys.path.insert(0, HERE)
from massfix import mass_fix, production_effect          # noqa: E402
from screen import endpoints, render                     # noqa: E402
from sharpen import de_endpoints                         # noqa: E402

WORK = os.path.join(HERE, "work")
NBINS = 24
TAUS = (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0)


def fit_curve(D, T, gate, mask, nbins=NBINS):
    """Equal-count bins of |D| -> isotonic E[t*sign(D) | |D|]. Returns (edges, m, xbar).

    Both sides are mean-removed across targets first: the Orion harnesses carry 26.4% / 44.1%
    of their truth energy in a target-independent vector, and leaving it in is what produced
    168b's discarded first pass (sign accuracies of 0.97 in the smallest bins).
    """
    Dm, Tm, gm = D[:, mask], T[:, mask], gate[mask]
    Dc = Dm - Dm.mean(0, keepdims=True)
    Tc = Tm - Tm.mean(0, keepdims=True)
    sel = np.zeros_like(Dm, dtype=bool)
    sel[:, gm] = True
    sel &= np.abs(Dm) > 0
    p, t = Dc[sel], Tc[sel]
    ap, al = np.abs(p), t * np.sign(p)
    edges = np.quantile(ap, np.linspace(0, 1, nbins + 1))
    edges[0], edges[-1] = 0.0, np.inf
    m, xbar = np.zeros(nbins), np.zeros(nbins)
    for i in range(nbins):
        s = (ap >= edges[i]) & (ap < edges[i + 1])
        if s.sum() < 200:
            m[i], xbar[i] = np.nan, np.nan
            continue
        m[i], xbar[i] = np.mean(al[s]), np.mean(ap[s])
    # fill gaps, clip at zero (a negative conditional mean means "no signal", not "flip"),
    # then enforce monotonicity -- the mechanism says reliability rises with magnitude, and
    # isotonic regression is how you impose that without imposing a functional form.
    ok = np.isfinite(m)
    m = np.interp(np.arange(nbins), np.flatnonzero(ok), m[ok])
    xbar = np.interp(np.arange(nbins), np.flatnonzero(ok), xbar[ok])
    m = np.maximum(m, 0.0)
    m = np.maximum.accumulate(m)
    return edges, m, xbar


def apply_curve(D, gate, edges, m, xbar, tau):
    """Multiply each gated coordinate by (m/|D|)^tau, normalised so tau=0 is the identity."""
    if tau == 0.0:
        return D.copy()
    G = np.where(gate[None, :], D, 0.0)
    a = np.abs(G)
    idx = np.clip(np.searchsorted(edges, a, side="right") - 1, 0, len(m) - 1)
    # the multiplier implied by the curve, on the bin's own representative magnitude
    with np.errstate(divide="ignore", invalid="ignore"):
        mult = np.where(xbar > 0, m / np.maximum(xbar, 1e-300), 0.0)
    mult = mult / max(float(np.max(mult)), 1e-300)      # scale-free; renorm restores the level
    w = mult[idx] ** tau
    out = G * w
    before = float(np.median(np.linalg.norm(G, axis=1)))
    after = float(np.median(np.linalg.norm(out, axis=1))) or 1.0
    return out * (before / after)


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    S = {}
    for name in ("hct116", "hek293t"):
        D, T, b, bref, gate, mask, ptar = production_effect(name)
        S[name] = dict(D=D, T=T, b=b, gate=gate, mask=mask, ptar=ptar)

    curves = {}
    for name, d in S.items():
        edges, m, xbar = fit_curve(d["D"], d["T"], d["gate"], d["mask"])
        curves[name] = (edges, m, xbar)
        rel = m / np.maximum(xbar, 1e-300)
        print(f"[{name}] curve fitted: multiplier m/|D| spans "
              f"{rel[rel > 0].min():.3g} .. {rel.max():.3g} over {NBINS} bins "
              f"(ratio {rel.max() / max(rel[rel > 0].min(), 1e-300):.1f}x)")

    out = {}
    for scored in ("hct116", "hek293t"):
        fitted_on = [o for o in S if o != scored][0]
        d = S[scored]
        base_eff = mass_fix(d["D"], d["b"], d["gate"])
        base = endpoints(render(base_eff, d["b"])[:, d["mask"]], d["T"][:, d["mask"]])
        bde = de_endpoints(scored, base_eff, d["b"], d["gate"], d["mask"], d["ptar"])
        bn = bde.get("nmae_exact", float("nan"))
        print(f"\n{'=' * 104}\nscored on {scored}, curve fitted on {fitted_on} "
              f"(NO truth from {scored} used to build the transform)")
        print(f"  {'tau':>5} {'pds':>9} {'d pds':>9} {'cos':>9} {'d cos %':>9} {'sse*':>10} "
              f"{'nmae':>8} {'d nmae':>9} {'pds->ovr':>9} {'nmae->ovr':>10}")
        rows = []
        for tau in TAUS:
            eff = mass_fix(apply_curve(d["D"], d["gate"], *curves[fitted_on], tau),
                           d["b"], d["gate"])
            e = endpoints(render(eff, d["b"])[:, d["mask"]], d["T"][:, d["mask"]])
            de = de_endpoints(scored, eff, d["b"], d["gate"], d["mask"], d["ptar"])
            nm = de.get("nmae_exact", float("nan"))
            dp, dn = e["pds"] - base["pds"], nm - bn
            rows.append(dict(tau=tau, **e, **de, d_pds=dp))
            print(f"  {tau:>5.2f} {e['pds']:9.5f} {dp:+9.5f} {e['cos']:9.5f} "
                  f"{100 * (e['cos'] / base['cos'] - 1):+8.1f}% {e['sse']:10.6f} "
                  f"{nm:8.4f} {dn:+9.4f} {dp * 0.36867:+9.5f} "
                  f"{(-dn * 0.2778) if np.isfinite(dn) else float('nan'):+10.5f}")
        out[scored] = dict(fitted_on=fitted_on, rows=rows)

    with open(os.path.join(WORK, "isotonic.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"\nwrote {os.path.join(WORK, 'isotonic.json')}")
    print("\nCOMPARISON, worst-surface d pds (from 168b, same instruments):")
    print("  GLOBAL q=1.15  +0.01145 | top 2% x2  +0.02602 | top 2% x3  +0.03452")


if __name__ == "__main__":
    sys.exit(main())
