"""241 — an UNPLANNED known-answer row: 214 already measured two of these rows on the
shipped renderer, and this experiment must reproduce them.

214 §1.1/§1.2 swept PerturbMatch and the CD4 dose on the SHIPPED renderer (mass fix unset),
4 PLOLO folds × challenge frames A/B/C, γ 0.65, α 0.982.  Those are exactly this experiment's
`NOMF` and `SHIPPED` bases, reached by independently written code.  If they disagree, one of
the two is wrong and neither table should be read.

214's rows used here (RESULTS.md §1.1 table and §1.2's two tables):
  PerturbMatch `pds` value, shipped renderer   f = 0.10 -> +0.002856 raw / +0.001065 overall
                                               f = 0.34 -> +0.003797 raw / +0.001415 overall
  PerturbMatch `pds` value, a23's renderer     f = 0.10 -> +0.005058 raw / +0.001885 overall
  CD4 dose value at pm 0.50                    f = 0.10 -> +0.007302 raw
                                               f = 0.34 -> +0.007260 raw
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
FOLDS4 = ("orion_hct116", "orion_hek293t", "k562_gwps", "kolf2")

EXPECT = [
    # (label, base tag, arm, 214's raw |d pds|, 214's overall, tol_raw)
    ("PerturbMatch, a23 renderer, f=0.10", "MF", "noPM", 0.005058, 0.001885, 5e-6),
    ("PerturbMatch, SHIPPED renderer, f=0.10", "NOMF", "noPM", 0.002856, 0.001065, 5e-6),
    ("PerturbMatch, SHIPPED renderer, f=0.34", "SHIPPED", "noPM", 0.003797, 0.001415, 5e-6),
    ("CD4 residual at pm 0.50, f=0.10", "NOMF", "noCD4", 0.007302, None, 5e-6),
    ("CD4 residual at pm 0.50, f=0.34", "SHIPPED", "noCD4", 0.007260, None, 5e-6),
]


def main() -> None:
    rows = json.loads(open(os.path.join(WORK, "pds241.json")).read())
    R = {(r["fold"], r["arm"], r["space"]): r for r in rows}

    def chal(fold, arm):
        vs = [R[(fold, arm, f"chal{c}/delivered")]["pds"] for c in "ABC"
              if (fold, arm, f"chal{c}/delivered") in R]
        return float(np.mean(vs)) if len(vs) == 3 else None

    print("=" * 108)
    print("CROSS-EXPERIMENT KNOWN-ANSWER ROW — 214's shipped-renderer PLOLO numbers")
    print("=" * 108)
    print(f"  {'row':<44}{'241 raw':>11}{'214 raw':>11}{'|d|':>10}"
          f"{'241 ovr':>11}{'214 ovr':>11}   ")
    ok_all, out = True, {}
    for label, tag, arm, want_raw, want_ovr, tol in EXPECT:
        per = [chal(f, f"{tag}:{arm}") - chal(f, f"{tag}:BASE") for f in FOLDS4]
        raw = -float(np.mean(per))                    # value of HAVING the component
        ovr = raw / (0.95228 - 0.50020) * 1.011 / 6
        d = abs(raw - want_raw)
        ok = d < tol
        ok_all &= ok
        out[label] = dict(raw=raw, overall=ovr, expect_raw=want_raw, expect_overall=want_ovr,
                          passed=bool(ok))
        print(f"  {label:<44}{raw:>+11.6f}{want_raw:>+11.6f}{d:>10.2e}"
              f"{ovr:>+11.6f}"
              + (f"{want_ovr:>+11.6f}" if want_ovr is not None else f"{'—':>11}")
              + ("   PASS" if ok else "   FAIL"))
    print(f"\n  VERDICT: {'PASS' if ok_all else 'FAIL'} — 241's PLOLO reproduces 214's "
          f"independently written shipped-renderer sweep")
    print("  (and 219's '214's PerturbMatch is 1.77x smaller on the shipped renderer' is "
          f"exactly {out['PerturbMatch, a23 renderer, f=0.10']['raw'] / out['PerturbMatch, SHIPPED renderer, f=0.10']['raw']:.2f}x)")
    json.dump(out, open(os.path.join(WORK, "ka214_241.json"), "w"), indent=1, default=float)
    print("\nwrote work/ka214_241.json")


if __name__ == "__main__":
    sys.exit(main())
