"""159 SS6 -- put PRODUCTION on the ceiling curve, on one surface, in one gene axis.

Everything else in this experiment measures the ceiling.  This measures the distance to it,
and it has to be done on a single surface or the comparison is meaningless.

The surface is the ORION HCT116 DONOR TABLE (all cells, 165,931 controls, its own 38,606-gene
axis restricted to the challenge's 18,533) rather than the thinned harness, because
`reconcile.py` shows the harness puts 87% of its effect energy into a cell-pool blob that the
donor table does not have.  The prediction is `experiments/104-context-heldout-credit-adapter/
work/hct116/frozen.npz`, the same frozen production artifact experiment 156 scored -- built
from the LEGAL donor pair (k562_gwps + orion_hek293t), so scoring it against orion_hct116 is
leakage-free.

Compared, on identical targets and genes, all mean-removed (the across-target mean is the
metric's zero and no predictor is credited for it):

  production                the shipped effect matrix
  1 donor / 2 / 3 / 4       equal-weight consensus of unit-normalised donor rows
  k -> infinity             the fit 1/cos^2 = A + B/k extrapolated, i.e. the across-line mean
  oracle same-line          the reference's own truth (cos = 1 by construction; the pds it
                            gets is the surface's own retrieval ceiling)
"""
from __future__ import annotations

import itertools
import json
import os

import numpy as np

from ceiling import endpoints, unit

SEED = 20260902
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
FROZEN = os.path.join(REPO, "experiments/104-context-heldout-credit-adapter/work/%s/frozen.npz")

# Only donors that actually CARRY the 300-target challenge panel can appear here.  099 SS1.1:
# rpe1 / nadig_jurkat / k562_essential are common-essential libraries and carry 0/300, because
# the 2026 panel was built to exclude common-essential genes.  That is also why this surface
# has four independent lines and not six.
REF = {"hct116": "orion_hct116", "hek293t": "orion_hek293t"}
DONORS = {"orion_hct116": ["data/donors/orion_hct116"],
          "orion_hek293t": ["data/donors/orion_hek293t"],
          "k562_gwps": ["data/donors/k562_gwps"],
          "kolf2": ["data/external/kolf2"],
          # the three CD4 conditions are ONE cell type; averaged so the consensus counts lines
          "cd4t": ["data/external/cd4t/cloud-compact/cd4_compact/cd4_rest",
                   "data/external/cd4t/cloud-compact/cd4_compact/cd4_stim8hr",
                   "data/external/cd4t/cloud-compact/cd4_compact/cd4_stim48hr"]}
BANNED = {"hct116": {"orion_hct116"}, "hek293t": {"orion_hek293t"}}


def load(rel):
    p = os.path.join(REPO, rel)
    return (np.nan_to_num(np.load(os.path.join(p, "delta.npy")).astype(np.float32)),
            [l.strip() for l in open(os.path.join(p, "perts.txt")) if l.strip()],
            [l.strip() for l in open(os.path.join(p, "genes.txt")) if l.strip()])


def kinf(ks, coss):
    """1/cos^2 = A + B/k  (exact for 'shared signal + independent per-donor component');
    returns (cos_infinity, R^2)."""
    x = 1.0 / np.asarray(ks, float)
    y = 1.0 / np.asarray(coss, float) ** 2
    A = np.column_stack([np.ones(len(x)), x])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    r2 = 1 - ((y - pred) ** 2).sum() / max(((y - y.mean()) ** 2).sum(), 1e-30)
    return float(1.0 / np.sqrt(max(coef[0], 1e-9))), float(r2)


def main() -> None:
    rng = np.random.default_rng(SEED)
    out = {}
    for surface in ("hct116", "hek293t"):
        z = np.load(FROZEN % surface, allow_pickle=True)
        ptar = [str(x) for x in z["targets"]]
        pgen = [str(x) for x in z["genes"]]
        raw = {k: [load(r) for r in v] for k, v in DONORS.items()}
        usable = [k for k in DONORS if k not in BANNED[surface]]

        tsets = [set.intersection(*[set(r[1]) for r in v]) for v in raw.values()]
        gsets = [set.intersection(*[set(r[2]) for r in v]) for v in raw.values()]
        ts = sorted(set(ptar).intersection(*tsets))
        gs = sorted(set(pgen).intersection(*gsets) - set(ts))
        print(f"\n=== {surface}: {len(ts)} targets x {len(gs)} genes "
              f"(production panel {len(ptar)}, donor lines {usable})")

        ip = {t: i for i, t in enumerate(ptar)}
        ig = {g: i for i, g in enumerate(pgen)}
        P = z["baseline"].astype(np.float64)[np.asarray([ip[t] for t in ts], dtype=np.int64)][
            :, np.asarray([ig[g] for g in gs], dtype=np.int64)]
        P = P - P.mean(0, keepdims=True)

        E = {}
        for k, mats in raw.items():
            acc = []
            for d, pl, gl in mats:
                a = {p: i for i, p in enumerate(pl)}
                b = {g: i for i, g in enumerate(gl)}
                v = d[np.asarray([a[t] for t in ts], dtype=np.int64)][
                    :, np.asarray([b[g] for g in gs], dtype=np.int64)].astype(np.float64)
                acc.append(unit(v - v.mean(0, keepdims=True)))
            v = np.mean(acc, axis=0)
            E[k] = v - v.mean(0, keepdims=True)
        T = E[REF[surface]]

        rows = [("production (shipped effect matrix)", endpoints(P, T))]
        cons = []
        for k in range(1, len(usable) + 1):
            combos = list(itertools.combinations(usable, k))
            if len(combos) > 10:
                combos = [combos[i] for i in rng.choice(len(combos), 10, replace=False)]
            vals = [endpoints(np.mean([unit(E[o]) for o in c], axis=0), T) for c in combos]
            e = {m: float(np.mean([v[m] for v in vals])) for m in vals[0]}
            cons.append((k, e))
            rows.append((f"consensus of {k} donor(s)  [{len(combos)} subsets]", e))
        ci, r2 = kinf([k for k, _ in cons], [e["cos"] for _, e in cons])
        rows.append(("oracle: the reference's own truth", endpoints(T, T)))

        print(f"  {'arm':<44}{'pds':>9}{'cos':>10}{'sse*':>10}{'row-row':>10}")
        for lbl, e in rows:
            print(f"  {lbl:<44}{e['pds']:9.4f}{e['cos']:10.5f}{e['sse']:10.6f}{e['rr']:10.5f}")
        print(f"  {'k -> infinity (fit 1/cos^2 = A + B/k)':<44}{'':>9}{ci:10.5f}"
              f"{1 - ci ** 2:10.6f}   R2 {r2:.5f}")
        prod = rows[0][1]["cos"]
        print(f"\n  PRODUCTION is at {prod / ci:.1%} of the k->infinity "
              f"context-independent ceiling on this surface "
              f"({prod:.5f} of {ci:.5f}),")
        print(f"  and at {prod / cons[-1][1]['cos']:.1%} of the best FINITE consensus we can "
              f"actually build here ({cons[-1][1]['cos']:.5f}, {len(usable)} donors).")
        # the pds-vs-transfer-level curve ON THIS SURFACE, so the placement and the curve
        # share a geometry.  P_t = tau*s_t + sqrt(1-tau^2)*z_t with z_t a random combination
        # of the OTHER targets' truths, orthogonalised against s_t.
        S = unit(T)
        W = rng.normal(size=(len(ts), len(ts)))
        np.fill_diagonal(W, 0.0)
        Z = W @ S
        Z = unit(Z - (Z * S).sum(1, keepdims=True) * S)
        curve = []
        print(f"\n  pds-vs-cosine curve on this surface ({len(ts)} targets, {len(gs)} genes)")
        print(f"    {'cos':>7}{'pds':>9}{'sse*':>10}   marker")
        marks = {round(rows[0][1]["cos"], 3): "<- PRODUCTION",
                 round(cons[-1][1]["cos"], 3): f"<- {len(usable)}-donor consensus",
                 round(ci, 3): "<- k->inf CEILING"}
        for tau in (0.0, 0.02, 0.03, 0.038, 0.043, 0.05, 0.072, 0.10, 0.123, 0.15,
                    0.20, 0.25, 0.30, 0.40, 0.50):
            e = endpoints(tau * S + np.sqrt(max(1 - tau ** 2, 0.0)) * Z, T)
            curve.append(dict(tau=tau, **e))
            print(f"    {tau:7.3f}{e['pds']:9.4f}{e['sse']:10.6f}   "
                  f"{marks.get(round(tau, 3), '')}")

        out[surface] = {"n_targets": len(ts), "n_genes": len(gs), "donors": usable,
                        "curve": curve,
                        "rows": [{"arm": l, **e} for l, e in rows],
                        "cos_kinf": ci, "kinf_r2": r2,
                        "production_frac_of_kinf": prod / ci,
                        "production_frac_of_best_finite": prod / cons[-1][1]["cos"]}

    with open(os.path.join(WORK, "place.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {os.path.join(WORK, 'place.json')}")


if __name__ == "__main__":
    main()
