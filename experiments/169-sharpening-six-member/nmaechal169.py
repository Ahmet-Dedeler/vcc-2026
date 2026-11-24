"""169 — `nmae` in the CHALLENGE's regime: the saturation-free limit, on both surfaces.

`regime169.py` established that the two surfaces disagree on `nmae` because they sit at
different AMPLITUDES against their own reference, not because sharpening does different
things on them:

    |h| <= |r|                challenge 99.3%   HCT116 98.88%   h1-test 91.01%
    median |h|/|r|            challenge 0.14-0.17   HCT116 0.1266   h1-test 0.2825
    L1 on saturated genes                       HCT116  8.80%   h1-test 22.62%

`P(sign correct)` RISES with `|h|` on both (0.51 -> 0.83 and 0.52 -> 0.62), exactly as 154
measured.  What falls on h1-test is `P(|h| < |r|)`, i.e. the head runs out of headroom.  So
h1-test's `nmae` explosion is a property of ITS amplitude, and the challenge does not have it.

THE SATURATION-FREE LIMIT IS EXACTLY THE LAW THE BOARD ALREADY FOLLOWS
---------------------------------------------------------------------
With `m_g = min(|h_g|,|r_g|)*1[sign match]`, `nmae_p = 1 + (Sum|h| - 2*Sum m)/Sum|r|`.  When
`|h| <= |r|` everywhere, `m = |h|*1[correct]` and the whole thing collapses to

        nmae_p  =  1  -  <sign(r_p), h_p> / ||r_p||_1                                  (U)

which is precisely what `docs/CURRENT_STATE.md` records for the challenge, "residual < 0.004".
(U) uses only the reference's SIGNS, never its magnitudes, so it is the one form whose transfer
does not depend on the harness's amplitude being right.

Reported per surface: the true `nmae` (with saturation), the saturation-free `nmae` (U), and
their deltas.  If the two surfaces agree in (U) and disagree only in the saturated form, the
disagreement is the amplitude and (U) is the number to carry to the board.

Also reported: an amplitude-matched sweep -- rescale `h` by `c` until the surface's saturated
L1 share matches the challenge's, and read the saturated metric there.  Two independent routes
to the same quantity; both are printed.

Deterministic.  No emitter, no cell-eval2, no build, no submission.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, HERE)
from cliff169 import QS, surface_h1test, surface_hct116  # noqa: E402

BOARD_EXCESS = 0.944251 - 1.0     # 096/A15-BOARD-RESULT.md raw nmae for a15
CS = (1.3, 1.15, 1.07, 1.0, 0.85, 0.7, 0.5, 0.35, 0.25)


def both(s, M, c=1.0):
    """(true nmae, saturation-free nmae, saturated-L1 share) at amplitude `c`."""
    t, u, satL1, allL1 = [], [], 0.0, 0.0
    for p in range(len(s.perts)):
        i = np.flatnonzero(s.real_sig[p] & np.isfinite(s.real_lfc[p]))
        if i.size < 10:
            continue
        r = s.real_lfc[p, i]
        h = c * np.where(np.isfinite(M[p, i]), M[p, i], 0.0)
        den = np.abs(r).sum()
        if not np.isfinite(den) or den == 0:
            continue
        same = np.sign(h) == np.sign(r)
        m = np.where(same, np.minimum(np.abs(h), np.abs(r)), 0.0)
        t.append((np.abs(h).sum() - 2 * m.sum()) / den)
        u.append(-(h * np.sign(r)).sum() / den)
        sat = np.abs(h) > np.abs(r)
        satL1 += np.abs(h)[sat].sum()
        allL1 += np.abs(h).sum()
    return 1 + float(np.mean(t)), 1 + float(np.mean(u)), satL1 / allL1


def main() -> None:
    out = {}
    for maker in (surface_hct116, surface_h1test):
        got = maker()
        if got is None:
            print("\nh1-test not built -- skipping")
            continue
        name, s, lfc = got
        Ms = {q: lfc(q) for q in QS}
        print(f"\n{'=' * 96}\n{name}")
        print(f"  {'q':>5} {'nmae true':>11} {'d true':>10} {'nmae unsat':>12} "
              f"{'d unsat':>10} {'sat L1 share':>13}")
        rows = []
        base = None
        for q in QS:
            t, u, sh = both(s, Ms[q])
            if base is None:
                base = (t, u)
            rows.append(dict(q=q, true=t, unsat=u, sat=sh))
            print(f"  {q:>5.2f} {t:11.5f} {t - base[0]:+10.5f} {u:12.5f} "
                  f"{u - base[1]:+10.5f} {sh:13.4f}")
        out[name] = rows

        print(f"\n  amplitude-matched: rescale h by c and re-read the SATURATED metric.")
        print(f"  the board's own excess is {BOARD_EXCESS:+.5f}; the row whose "
              f"`nmae true - 1` matches it\n  is the amplitude at which this surface is the "
              f"challenge's analogue.")
        print(f"  {'c':>5} {'nmae true':>11} {'excess':>9} {'sat L1':>8} "
              + " ".join(f"{'d@q=' + f'{q:.2f}':>11}" for q in (1.10, 1.15, 1.20, 1.25)))
        for c in CS:
            t0, _, sh = both(s, Ms[1.00], c)
            cells = []
            for q in (1.10, 1.15, 1.20, 1.25):
                tq, _, _ = both(s, Ms[q], c)
                cells.append(f"{tq - t0:+11.5f}")
            print(f"  {c:>5.2f} {t0:11.5f} {t0 - 1:+9.5f} {sh:8.4f} " + " ".join(cells))

    if len(out) == 2:
        print(f"\n{'=' * 96}\nCROSS-SURFACE AGREEMENT in the saturation-free limit (U)")
        print(f"  {'q':>5} {'d unsat hct116':>16} {'d unsat h1-test':>17} "
              f"{'worst':>9} | {'d TRUE hct116':>15} {'d TRUE h1-test':>16}")
        for i, q in enumerate(QS):
            a = out["hct116"][i]["unsat"] - out["hct116"][0]["unsat"]
            b = out["h1-test"][i]["unsat"] - out["h1-test"][0]["unsat"]
            at = out["hct116"][i]["true"] - out["hct116"][0]["true"]
            bt = out["h1-test"][i]["true"] - out["h1-test"][0]["true"]
            print(f"  {q:>5.2f} {a:+16.5f} {b:+17.5f} {max(a, b):+9.5f} | "
                  f"{at:+15.5f} {bt:+16.5f}")

    with open(os.path.join(WORK, "nmaechal169.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"\nwrote {os.path.join(WORK, 'nmaechal169.json')}")


if __name__ == "__main__":
    sys.exit(main())
