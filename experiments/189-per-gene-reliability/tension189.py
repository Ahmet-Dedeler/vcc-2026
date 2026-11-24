"""189 — resolving the 171-vs-pilot tension about expression.

THE TENSION
-----------
* The pilot: control expression is the strongest per-gene cross-line reliability signal on the
  board of candidates, spearman +0.5215 against `r_HE` and a **+0.0505** held-out decile spread.
* 171 §12c: our emitted p-order is expression-SET (Spearman(-log10 p_adj, expression) = +0.485
  inside the call set, against -0.391 for `|H|`) and treated that as the defect to remove.

If expression genuinely predicts which genes transfer, part of what 171 called a defect is
signal.  This script decides it with three measurements, on both surfaces:

  A. h1-test, INSIDE the reference DE pool (what `reach`, `fid`, `nmae` all read): does
     expression order our per-gene sign accuracy?   [known-answer: INBOX-4 §1's flat row]
  B. h1-test, OUTSIDE the pool (adjudicable but not reference-significant): same question.
     If expression works here and not there, CONDITIONING ON SIGNIFICANCE is the mechanism.
  C. donor surface: re-run the pilot's per-gene agreement with the gate changed from "above
     this line's own median |effect|" to "in this target's top-K by the Wald proxy
     |delta|*sqrt(cpm)" (172 §1.1's validated DE-pool proxy, median recall 0.828).  That is
     the donor-side analogue of conditioning on significance.  If the expression signal
     collapses under it, the pilot's +0.5215 is a measurement-precision artefact of the
     UNCONDITIONED gate and does not survive into any member.
"""
from __future__ import annotations

import json
import os

import numpy as np

from common189 import (WORK, bvec, centre, decile_spread, gate_mask, load_table, pair_agreement,
                       spearman, submatrix)
from harness189 import Surface189, check

LINES = ["orion_hct116", "orion_hek293t", "k562_gwps"]


def decile_curve(key: np.ndarray, ok: np.ndarray, mask: np.ndarray, genes_sorted=None):
    """Pooled sign accuracy by decile of a per-(target, gene) key over `mask`."""
    num = np.zeros(10)
    den = np.zeros(10)
    for p in range(mask.shape[0]):
        idx = np.flatnonzero(mask[p])
        if idx.size < 10:
            continue
        k = np.asarray(key[p], dtype=np.float64)[idx]
        o = ok[p, idx].astype(np.float64)
        o = o[np.argsort(-k, kind="stable")]
        e = (np.arange(11) * idx.size) // 10
        for d in range(10):
            a, b = e[d], e[d + 1]
            if b > a:
                num[d] += o[a:b].sum()
                den[d] += b - a
    return num / np.maximum(den, 1)


def main() -> None:
    T = Surface189()
    check(T)
    out = {}

    expr = np.broadcast_to(T.b_ctrl, T.H.shape)
    absH = np.abs(T.H)
    inpool = T.pool
    outpool = T.S.adj & ~T.S.real_sig
    print(f"\n{'=' * 96}\nA/B — h1-test: does expression order our per-gene sign accuracy?")
    print(f"    pool (reference-significant AND adjudicable): "
          f"{int(inpool.sum()):,} (target, gene) pairs")
    print(f"    outside the pool (adjudicable, not significant): {int(outpool.sum()):,}")
    for tag, m in (("INSIDE the DE pool", inpool), ("OUTSIDE the pool", outpool)):
        for kn, kk in (("expression", expr), ("|H|", absH)):
            d = decile_curve(kk, T.correct, m)
            print(f"  {tag:<22} {kn:<11} d1..d10  " + " ".join(f"{x:.3f}" for x in d)
                  + f"   d1-d10 {d[0] - d[9]:+.4f}")
            out[f"{tag}|{kn}"] = list(map(float, d))

    # ---- C: the donor surface under a significance-like gate --------------------------
    print(f"\n{'=' * 96}\nC — donor surface: the pilot's gate vs a significance-like gate")
    tabs = {t: load_table(t) for t in LINES}
    genes = sorted(set.intersection(*[set(t["genes"]) for t in tabs.values()]))
    targets = sorted(set.intersection(*[set(t["perts"]) for t in tabs.values()]))
    D, B, CPM = {}, {}, {}
    for t in LINES:
        D[t] = centre(submatrix(tabs[t], targets, genes).astype(np.float32))
        B[t] = bvec(tabs[t], genes)
        CPM[t] = np.expm1(B[t])

    def topk_gate(A: np.ndarray, cpm: np.ndarray, frac: float) -> np.ndarray:
        """Per TARGET, the top `frac` of genes by the Poisson-Wald proxy |delta|*sqrt(cpm)."""
        z = np.abs(A) * np.sqrt(np.maximum(cpm, 0))[None, :]
        k = max(1, int(round(frac * A.shape[1])))
        g = np.zeros(A.shape, bool)
        part = np.argpartition(-z, k - 1, axis=1)[:, :k]
        np.put_along_axis(g, part, True, axis=1)
        return g

    rows = []
    hct, hek, k562 = LINES
    gates = {"pilot: > own per-gene median |effect|":
             {t: gate_mask(D[t]) for t in LINES}}
    for frac in (0.05, 0.01):
        gates[f"significance proxy: per-target top {frac:.0%} by |d|*sqrt(cpm)"] = \
            {t: topk_gate(D[t], CPM[t], frac) for t in LINES}

    print(f"{'gate':<52}{'genes':>7}{'rho(r_HE,expr)':>16}{'d10 expr':>10}"
          f"{'d10 r_HE':>10}{'mean r':>9}")
    for gname, G in gates.items():
        kHE, nHE = pair_agreement(D[hct], D[hek], G[hct], G[hek])
        kHK, nHK = pair_agreement(D[hct], D[k562], G[hct], G[k562])
        sup = (nHE > 50) & (nHK > 50)
        rHE, rHK = kHE[sup] / nHE[sup], kHK[sup] / nHK[sup]
        b = B[hct][sup]
        s = spearman(rHE, b)
        _, _, d_expr = decile_spread(b, rHK)
        _, _, d_rHE = decile_spread(rHE, rHK)
        rows.append(dict(gate=gname, genes=int(sup.sum()), rho_expr=s, d_expr=d_expr,
                         d_rHE=d_rHE, mean_r=float(rHE.mean())))
        print(f"{gname:<52}{int(sup.sum()):>7}{s:>+16.4f}{d_expr:>+10.4f}{d_rHE:>+10.4f}"
              f"{rHE.mean():>9.4f}")

    out["donor_gates"] = rows
    json.dump(out, open(os.path.join(WORK, "tension189.json"), "w"), indent=2)
    print("\nwrote work/tension189.json")


if __name__ == "__main__":
    main()
