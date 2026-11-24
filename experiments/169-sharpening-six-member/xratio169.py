"""169 — how sharpening moves the three `mse` terms `A`, `B` and `X`, on surfaces with truth.

`mse - 1 = ( s^2*A - 2*s*X + Kc - min(s^2*B, claim) ) / sumD`.

`A`, `B` and `claim` are all computable on the REAL challenge panel with no reference
(`claim169.py` does that).  `X = Sum_p w_p <dh_p, d_p>` is the one term that needs the
withheld truth, so it can only be transferred -- and 156/162 warn explicitly against pricing
a profile arm by scaling a harness cosine.  What is measured here is therefore the RATIO
`rX(q) = X_q / X_1` on the two held Orion lines, reported per surface and never averaged, so
the claim being transferred is the one thing a ratio can carry: the direction of the move.

`claim169.py` prices `mse` at THREE settings -- `rX = 1` (the null: sharpening does not touch
alignment), and the two measured ratios -- and the GO/NO-GO uses the null.

Also reported: `rA(q)` and `rB(q)`, which are the same quantities the challenge panel gives
directly, so the two routes can be cross-checked against each other.

Deterministic.  No emitter, no cell-eval2, no build, no submission.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(REPO, "experiments/156-delivered-geometry"))
sys.path.insert(0, HERE)
from massfix import mass_fix, production_effect     # noqa: E402
from screen import render                          # noqa: E402
from sharp import sharpen, sharpen_lfc             # noqa: E402

QS = (1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30)


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    out = {}
    for surface in ("hct116", "hek293t"):
        D, T, b, b_ref, gate, mask, ptar = production_effect(surface)
        print(f"\n{'=' * 96}\n{surface}: {D.shape[0]} targets, {len(mask)} scored genes")
        for mode in ("row", "median", "lfc", "lfcL"):
            print(f"\n  --- renormalisation: {mode} ---")
            print(f"  {'q':>5} {'A':>12} {'rA':>9} {'X':>12} {'rX':>9} "
                  f"{'cos=X/sqrt(A.R)':>17} {'r cos':>9}")
            base = None
            rows = []
            for q in QS:
                if mode in ('lfc', 'lfcL'):
                    t = sharpen_lfc(D, q, b, gate, mode='row',
                                    restore='h' if mode == 'lfc' else 'lfc')
                else:
                    t = sharpen(D, q, gate, mode=mode)
                eff = mass_fix(t, b, gate)
                P = render(eff, b)[:, mask]
                Tm = T[:, mask]
                A = float((P * P).sum())
                X = float((P * Tm).sum())
                R = float((Tm * Tm).sum())
                rec = dict(q=q, A=A, X=X, R=R, cos=X / np.sqrt(A * R))
                if base is None:
                    base = rec
                rows.append(rec)
                print(f"  {q:>5.2f} {A:12.5f} {A / base['A']:9.5f} {X:12.5f} "
                      f"{X / base['X']:9.5f} {rec['cos']:17.6f} "
                      f"{rec['cos'] / base['cos']:9.5f}")
            out[f"{surface}|{mode}"] = rows

    with open(os.path.join(WORK, "xratio169.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)

    print(f"\n{'=' * 96}\nTRANSFERRED RATIOS (rule 1: report both, never average, "
          f"take the smaller for a GO)")
    for mode in ("row", "median", "lfc", "lfcL"):
        print(f"\n  {mode}:  {'q':>5} {'rX hct116':>11} {'rX hek293t':>12} {'rX min':>9} "
              f"{'rA hct116':>11} {'rA hek293t':>12}")
        h, k = out[f"hct116|{mode}"], out[f"hek293t|{mode}"]
        for i, q in enumerate(QS):
            rxh, rxk = h[i]["X"] / h[0]["X"], k[i]["X"] / k[0]["X"]
            print(f"         {q:>5.2f} {rxh:>11.5f} {rxk:>12.5f} {min(rxh, rxk):>9.5f} "
                  f"{h[i]['A'] / h[0]['A']:>11.5f} {k[i]['A'] / k[0]['A']:>12.5f}")


if __name__ == "__main__":
    sys.exit(main())
