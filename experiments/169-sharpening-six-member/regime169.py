"""169 — WHICH surface is the `nmae` analogue of the challenge?  The two disagree, so decide it.

`cliff169.py` measured `nmae = 1 + (Sum|h| - 2*Sum m)/Sum|r|` exactly, and then measured the
marginal efficiency `e = P(sign correct AND |h| < |r|)` by decile of `|h|`.  The two surfaces
disagree about its SIGN in the head:

    HCT116   e rises  0.504 -> 0.626 across the deciles, every decile above 1/2
    h1-test  e falls  0.513 -> 0.362 across the deciles, deciles 7-10 below 1/2

`e` factorises as `P(correct) * P(|h| < |r| | correct)`, and only the second factor is a
REGIME property -- how big our predicted fold changes are against the reference's.  The
challenge's regime is recorded, independently of this experiment, in `docs/CURRENT_STATE.md`:

    "our |lfc_hat| is 0.14-0.17x the reference's on its own gate with |h| <= |r| on 99.3%
     of gated genes"

and the same file records the consequence, `NMAE = 1 - <sign(lfc_ref), lfc_hat>` to a residual
below 0.004 -- which is EXACTLY the unsaturated limit of the identity above, since `m = |h|`
whenever `|h| <= |r|` and the correct signs then give `Sum|h| - 2*Sum m = -<sign(r), h>`.

So the challenge regime is fully unsaturated, and the discriminating measurement is the
saturated fraction and the `|h|/|r|` level on each surface.  That is what this script reports,
next to the board's own raw `nmae` for a15 (0.944251, `A15-BOARD-RESULT.md`).

Deterministic.  No emitter, no cell-eval2, no build, no submission.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cliff169 import NDEC, surface_h1test, surface_hct116  # noqa: E402

BOARD_NMAE_A15 = 0.944251        # 096/A15-BOARD-RESULT.md, raw
CHALLENGE_UNSAT = 0.993          # docs/CURRENT_STATE.md, "|h| <= |r| on 99.3% of gated genes"
CHALLENGE_RATIO = (0.14, 0.17)   # docs/CURRENT_STATE.md, "|lfc_hat| is 0.14-0.17x"


def main() -> None:
    print(f"challenge regime, recorded independently of this experiment "
          f"(docs/CURRENT_STATE.md):")
    print(f"  |h| <= |r| on {CHALLENGE_UNSAT * 100:.1f}% of gated genes; "
          f"|h|/|r| = {CHALLENGE_RATIO[0]}-{CHALLENGE_RATIO[1]}x; "
          f"board raw nmae (a15) = {BOARD_NMAE_A15:.6f}")

    for maker in (surface_hct116, surface_h1test):
        got = maker()
        if got is None:
            print("\nh1-test not built -- skipping")
            continue
        name, s, lfc = got
        M = lfc(1.0)
        allh, allr, corr = [], [], []
        for p in range(len(s.perts)):
            i = np.flatnonzero(s.real_sig[p] & np.isfinite(s.real_lfc[p]))
            if i.size < 10:
                continue
            r = s.real_lfc[p, i]
            h = np.where(np.isfinite(M[p, i]), M[p, i], 0.0)
            allh.append(np.abs(h))
            allr.append(np.abs(r))
            corr.append(np.sign(h) == np.sign(r))
        H = np.concatenate(allh)
        R = np.concatenate(allr)
        C = np.concatenate(corr)
        unsat = H <= R
        ratio = H / np.maximum(R, 1e-12)
        print(f"\n{name}: {H.size} (perturbation, gene) pairs on the reference's DE set")
        print(f"  nmae at q=1                 {s.nmae(pred_lfc=M):.6f}   "
              f"(board a15 {BOARD_NMAE_A15:.6f})")
        print(f"  |h| <= |r|                  {100 * unsat.mean():.2f}%   "
              f"(challenge {CHALLENGE_UNSAT * 100:.1f}%)")
        print(f"  L1 share on saturated genes {100 * H[~unsat].sum() / H.sum():.2f}%")
        print(f"  median |h|/|r|              {np.median(ratio):.4f}   "
              f"(challenge {CHALLENGE_RATIO[0]}-{CHALLENGE_RATIO[1]})")
        print(f"  Sum|h| / Sum|r|             {H.sum() / R.sum():.4f}")
        print(f"  L1-weighted sign accuracy   {(H * C).sum() / H.sum():.4f}   "
              f"(unweighted {C.mean():.4f})")
        print(f"  2*overlap / budget          "
              f"{2 * np.where(C, np.minimum(H, R), 0.0).sum() / H.sum():.4f}")

        # the factorisation: which of the two factors makes `e` fall on h1-test?
        k = np.clip((np.searchsorted(np.sort(H), H, side="right") * NDEC) // H.size,
                    0, NDEC - 1)
        print(f"  {'decile':>7} {'L1 share':>9} {'P(correct)':>11} {'P(unsat)':>10} "
              f"{'P(unsat|corr)':>14} {'e':>8}")
        for d in range(NDEC):
            m = k == d
            w = H[m]
            print(f"  {d + 1:>7} {w.sum() / H.sum():9.4f} "
                  f"{(w * C[m]).sum() / w.sum():11.4f} "
                  f"{(w * unsat[m]).sum() / w.sum():10.4f} "
                  f"{(w * (unsat[m] & C[m])).sum() / max((w * C[m]).sum(), 1e-30):14.4f} "
                  f"{(w * (unsat[m] & C[m])).sum() / w.sum():8.4f}")


if __name__ == "__main__":
    sys.exit(main())
