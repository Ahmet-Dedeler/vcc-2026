"""161 SS5 — is the fitted gain CONTEXT TRANSFER, or just a re-shaped CONSENSUS?

THE PROBLEM THIS SOLVES
-----------------------
`residual.py`'s ablation puts the whole gain in the `T` family -- features computable
WITHOUT the recipient's basal profile.  That is not a small detail: it means the winning map
is, by construction, a function of donor-side quantities only, and therefore cannot be
carrying information about the held line.  Whatever it does, it does identically to every
recipient.

So the honest reading of a positive `T` arm is NOT "we predicted the line-specific residual".
It is "the equal-weight consensus is the wrong per-entry shape, and a learned monotone
re-calibration of it scores better" -- which is a CONSENSUS improvement, the 41st of those,
not a residual model.

This script tests that reading directly, with three closed-form re-calibrations that use only
`C` and the cross-donor spread `dsd`, no recipient information and no learning:

    power      C' = sign(C) * |C|^p                      pure monotone per-entry reshaping
    soft       C' = sign(C) * max(|C| - lam*dsd, 0)      soft-threshold on donor disagreement
    eb         C' = C * t2 / (t2 + kap * dsd^2 / k)      per-entry empirical-Bayes shrinkage

If one of these recovers most of the GBM's `T` gain, the GBM rediscovered a per-entry
consensus recalibration and the residual hypothesis is unsupported by it.

NOTE (AGENTS.md §12): EB shrinkage and inverse-variance pooling are on the closed list.  They
are used HERE as a diagnostic to explain a measured gain, not as a shipping proposal -- and if
the GBM's gain equals theirs, that is evidence the GBM is re-deriving a closed lever.

Scored in the delivered space, at restored row norms, exactly as `residual.py`.
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

from residual import (Surface, _pds_contrib, boot_pds, clusters, deliver, endpoints,
                      restore_norms, split_par_orth)
from common import WORK


def arms(C, DSD, k):
    """Every arm is a function of the DONORS only -- no recipient information at all."""
    out = [("consensus (baseline)", C)]
    for p in (0.5, 0.75, 1.25, 1.5, 2.0):
        out.append((f"power p={p}", np.sign(C) * np.abs(C) ** p))
    md = float(np.median(DSD[DSD > 0])) if (DSD > 0).any() else 1.0
    for lam in (0.1, 0.25, 0.5, 1.0):
        out.append((f"soft lam={lam}", np.sign(C) * np.maximum(np.abs(C) - lam * DSD, 0.0)))
    t2 = float(np.mean(C ** 2))
    for kap in (0.25, 1.0, 4.0, 16.0):
        out.append((f"eb kappa={kap}", C * (t2 / (t2 + kap * DSD ** 2 / k))))
    for q in (0.5, 0.75, 0.9):
        thr = np.quantile(np.abs(C), q, axis=1, keepdims=True)
        out.append((f"hard keep-top {1 - q:.0%}", np.where(np.abs(C) >= thr, C, 0.0)))
    return out, md


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--surface", default="core5")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out_path = a.out or os.path.join(WORK, f"recal-{a.surface}.json")

    S = Surface(a.surface, ctrl=False)
    print(f"surface {a.surface}: {S.nT} x {S.nG}, lines {', '.join(S.lines)}")
    rows = []
    for L in S.lines:
        pool = [x for x in S.lines if x != L]
        C, DSD, _ = S.consensus(pool)
        b, gate, Etruth = S.bL[L], S.gateL[L], S.EL[L]
        dn = np.median([np.median(np.linalg.norm(S.R[d], axis=1)) for d in pool])
        alpha = dn / np.median(np.linalg.norm(C, axis=1))
        Pbase = deliver(alpha * C, b, gate)
        base = endpoints(Pbase, Etruth)
        base_contrib = _pds_contrib(Pbase, Etruth)
        lab = clusters(C)
        print(f"\n{'=' * 100}\nHELD LINE {L}   donors {pool}")
        print(f"  CONSENSUS BASELINE  pds {base['pds']:.5f}  cos {base['cos']:.5f}")
        A, md = arms(C, DSD, len(pool))
        for nm, Q in A:
            if nm.startswith("consensus"):
                rows.append(dict(line=L, arm=nm, **base, d_pds=0.0, d_cos=0.0))
                continue
            # same parameterisation as residual.py: the change is a DIRECTION change only,
            # and only its component orthogonal to C can move a row-normalised metric.
            eff = restore_norms(Q, C) * alpha
            e = endpoints(deliver(eff, b, gate), Etruth)
            _, u = split_par_orth(Q - C, C)
            rel = float(np.median(np.linalg.norm(u, axis=1) / np.linalg.norm(C, axis=1)))
            mu, lo, hi = boot_pds(deliver(eff, b, gate), base_contrib, Etruth, lab, n=200)
            r = dict(line=L, arm=nm, **e, d_pds=e["pds"] - base["pds"],
                     d_cos=e["cos"] - base["cos"], relnorm=rel,
                     boot_lo=lo, boot_hi=hi)
            rows.append(r)
            print(f"  {nm:<24} pds {e['pds']:.5f} ({r['d_pds']:+.5f})  "
                  f"cos {e['cos']:.5f} ({r['d_cos']:+.5f})  |orth|/|C| {rel:.3f}  "
                  f"boot[{lo:+.5f},{hi:+.5f}]")
    with open(out_path, "w") as fh:
        json.dump(dict(surface=a.surface, rows=rows), fh, indent=1, default=float)
    print(f"\nwrote {out_path}")

    print(f"\n{'=' * 100}\nWORST-LINE RANKING (the repo's selection rule)")
    names = sorted({r["arm"] for r in rows})
    tab = []
    for nm in names:
        v = [r["d_pds"] for r in rows if r["arm"] == nm]
        if len(v) == len(S.lines):
            tab.append((min(v), float(np.mean(v)), nm))
    tab.sort(reverse=True)
    print(f"  {'arm':<24} {'worst':>10} {'mean':>10}")
    for worst, mean, nm in tab:
        print(f"  {nm:<24} {worst:+10.5f} {mean:+10.5f}")


if __name__ == "__main__":
    main()
