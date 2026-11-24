"""168 — re-price magnitude sharpening with the instrument experiment 030 did not have.

## Why this is being reopened

Experiment 030 (2026-08-21) found that a one-parameter power transform of the predicted
per-gene magnitudes,

    P -> sign(P) * |P|^q ,   renormalised to preserve the median row norm

is a **large** `pds` lever: macro PDS 0.8017 -> 0.8415 at q = 2.5 on the LODO harness, +0.038 at
q = 2, and roughly additive with whitening (whitening alone +0.080, sharpening alone +0.072,
both +0.120 -- 79% additive, taking the incumbent 0.7378 -> 0.8415).

**It was never shipped, and its own write-up says why:**

> "Not shipped, and not resolvable with the tools that produced it. ... The `pds` gain is exact
> and large; the `nmae` cost is estimated by an instrument known to be miscalibrated on the axis
> that decides it."

The deciding axis was `nmae`, estimated through experiment 018's proxy, which put the
`nmae`-optimal alpha at ~0.75 while the count-level run measured 1.08 -- so the whole
decoupling argument rested on a number the proxy got wrong.

**Three things have changed since, and each closes part of that gap:**

1. `experiments/153-de-bench/debench.py:Surface.nmae()` computes `nmae` **exactly** and
   deterministically from the reference DE table -- no proxy, no emitter, no seed. 030 could not
   do this.
2. `experiments/156-delivered-geometry` established that the effect the scorer reads is not the
   effect we build, and supplies the renderer + `endpoints`. 030 scored the pre-render effect.
3. `mse` was not known to be a live member until experiment 096 (2026-08-28), a week after 030
   was parked. It is now our second-largest gap, and a13 proved the magnitude *distribution* is
   its strongest lever (raw `mse` 0.9013 -> 0.8493, `score_mse` +0.0548, the largest single-member
   gain in this project's history).

## Why sharpening should be the RIGHT version of a13's lever

a13 shrank by **donor agreement** and `nmae` collapsed to raw 1.0185 -- past its own baseline,
worse than predicting nothing. Its post-mortem named the mechanism: `nmae` scores fold-change
error over **the reference's** DE genes, and donor agreement is uncorrelated with which genes
those are, so the shrinkage landed squarely on the scored set. Measured here on HCT116, AUC for
predicting reference-DE membership:

    |predicted effect|            0.5966   (positive on 97.6% of targets)
    donor agreement (a13's key)   0.4980   (no information; below 0.5 on 57.5% of targets)

So sharpening concentrates magnitude on genes that **are** more likely to be scored, which is
the opposite of what killed a13. And experiment 154 measured top-decile-by-|effect| sign purity
at 0.675 against an overall 0.531, so the magnitude also predicts *correctness* -- which is the
quantity both `mse` and `nmae` are paid in.

A monotone power transform **preserves the |effect| ordering exactly**, so `keep_top` and every
|effect|-keyed policy are untouched by construction.

## What is measured here

Per arm, on every surface with truth, in the DELIVERED space (gate -> production
`mass_consistent_effect` -> render):
  pds, macro cosine, row-row cosine, sse* (the `mse` proxy at its own best amplitude),
  and -- the point of the experiment -- **exact `nmae`**, plus sign accuracy and head purity.
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
from debench import Surface                                    # noqa: E402
from massfix import GATE_CPM, mass_fix, production_effect       # noqa: E402
from screen import TS, endpoints, render                       # noqa: E402

WORK = os.path.join(HERE, "work")
QS = (1.0, 1.15, 1.25, 1.5, 1.75, 2.0, 2.5)
EPS = 1e-9

DE_TABLE = {"hct116": os.path.join(REPO, "experiments/095-clean-harness/de_thinned.npz")}
FROZEN = os.path.join(REPO, "experiments/104-context-heldout-credit-adapter/work/%s/frozen.npz")


def sharpen(D: np.ndarray, q: float, gate: np.ndarray) -> np.ndarray:
    """sign(D)*|D|^q on the gated coordinates, renormalised to the ORIGINAL median row norm.

    The renorm is what keeps `alpha` meaning one fixed thing: |D|^q changes the overall scale by
    orders of magnitude, exactly the reason `whiten.py:apply` renormalises too. Preserving the
    MEDIAN row norm (not each row's own) leaves the across-target magnitude spread free to move,
    which is the part of the transform being tested.
    """
    if q == 1.0:
        return D.copy()
    G = np.where(gate[None, :], D, 0.0)
    out = np.sign(G) * np.abs(G) ** q
    before = float(np.median(np.linalg.norm(G, axis=1)))
    after = float(np.median(np.linalg.norm(out, axis=1))) or 1.0
    return out * (before / after)


def de_endpoints(surface_name, eff, b, gate, mask, ptar):
    """Exact nmae + sign accuracy + head purity, via the cached reference DE table."""
    path = DE_TABLE.get(surface_name)
    if path is None or not os.path.exists(path):
        return {}
    s = Surface(path)
    z = np.load(FROZEN % surface_name, allow_pickle=True)
    ctrl_counts = np.load(os.path.join(REPO, "experiments/156-delivered-geometry/work",
                                       f"truth-{surface_name}.npz"),
                          allow_pickle=True)["ctrl_counts"].astype(np.float64)
    ctrl_cpm = 1e6 * ctrl_counts / ctrl_counts.sum()
    lin = np.clip(np.expm1(b[None, :] + eff), 0.0, None)
    p = lin / lin.sum(1, keepdims=True)
    L = np.log2((1e6 * p + EPS) / (ctrl_cpm[None, :] + EPS))

    gpos = {g: i for i, g in enumerate(s.genes)}
    cols = np.asarray([gpos.get(str(g), -1) for g in z["genes"]])
    prow = {pp: i for i, pp in enumerate(s.perts)}
    rows = np.asarray([prow.get(pp, -1) for pp in ptar])
    ok, keep = cols >= 0, rows >= 0
    M = np.full((len(s.perts), len(s.genes)), np.nan)
    M[np.ix_(rows[keep], cols[ok])] = L[np.ix_(np.flatnonzero(keep), np.flatnonzero(ok))]
    M = np.where(np.isfinite(M), M, 0.0)
    s.pred_lfc = M
    num = den = 0
    for pp in range(len(s.perts)):
        i = np.flatnonzero(s.real_sig[pp] & s.adj[pp])
        if i.size < 10:
            continue
        m = i[np.abs(M[pp, i]) >= np.percentile(np.abs(M[pp, i]), 90.0)]
        num += int((np.sign(M[pp, m]) == np.sign(s.real_lfc[pp, m])).sum())
        den += m.size
    return {"nmae_exact": s.nmae(), "sign_acc": s.sign_accuracy(),
            "head_purity": num / max(den, 1)}


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    out = {}
    for surface in ("hct116", "hek293t"):
        D, T, b, b_ref, gate, mask, ptar = production_effect(surface)
        print(f"\n{'=' * 118}\n{surface}: {D.shape[0]} targets, {int(gate.sum())} gated genes, "
              f"{len(mask)} scored genes")
        print(f"  {'q':>5} {'pds':>9} {'d pds':>9} {'cos':>9} {'d cos %':>9} {'row-row':>9} "
              f"{'sse*':>10} {'nmae':>9} {'d nmae':>9} {'sign':>7} {'head':>7}")
        base, rows = None, []
        for q in QS:
            eff = mass_fix(sharpen(D, q, gate), b, gate)
            scored = render(eff, b)
            e = endpoints(scored[:, mask], T[:, mask])
            de = de_endpoints(surface, eff, b, gate, mask, ptar)
            rec = dict(q=q, **e, **de)
            if base is None:
                base = rec
            rows.append(rec)
            nm = de.get("nmae_exact", float("nan"))
            print(f"  {q:>5.2f} {e['pds']:9.5f} {e['pds'] - base['pds']:+9.5f} "
                  f"{e['cos']:9.5f} {100 * (e['cos'] / base['cos'] - 1):+8.1f}% "
                  f"{e['rr']:9.5f} {e['sse']:10.6f} "
                  f"{nm:9.4f} {nm - base.get('nmae_exact', nm):+9.4f} "
                  f"{de.get('sign_acc', float('nan')):7.4f} "
                  f"{de.get('head_purity', float('nan')):7.4f}")
        out[surface] = rows

    with open(os.path.join(WORK, "sharpen.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"\n{'=' * 118}\nWORST-SURFACE d pds (the repo's selection rule)")
    print(f"  {'q':>5} {'worst d pds':>12} {'hct116':>10} {'hek293t':>10} "
          f"{'sse* hct':>10} {'sse* hek':>10} {'d nmae hct':>11}")
    for i, q in enumerate(QS):
        h, k = out["hct116"][i], out["hek293t"][i]
        dh = h["pds"] - out["hct116"][0]["pds"]
        dk = k["pds"] - out["hek293t"][0]["pds"]
        print(f"  {q:>5.2f} {min(dh, dk):+12.5f} {dh:+10.5f} {dk:+10.5f} "
              f"{h['sse']:10.6f} {k['sse']:10.6f} "
              f"{h.get('nmae_exact', float('nan')) - out['hct116'][0].get('nmae_exact', 0):+11.4f}")
    print("\nNOTE `nmae` here is EXACT (debench, reference DE table), not experiment 018's proxy.")
    print("     That proxy is the reason 030 was parked: it put the nmae-optimal alpha at ~0.75")
    print("     against a count-level measurement of 1.08.")


if __name__ == "__main__":
    sys.exit(main())
