"""192 — the COMPOSED v14 retrodiction and the final six-member ledger.

`retro192` and `close192` each retrodicted ONE of v14's three DE members and missed another,
because each modelled ONE of the two variables v14 moved.  v14 = v11 + `--tail-spread 1.0
--keep-top 850`, which changes the call set's SIZE (5,961 -> 1,563) AND its COMPOSITION
relative to our own `|predicted effect|` ordering.  The composed model uses:

  `fid`   from the SIZE channel  -- key-free, through the coverage term (`cover192`), because
          185 §4 measured precision flat at 0.507-0.517 for every key at every K, so coverage
          is the whole shape of `fid(n_pred)`;
  `jac`   from the COMPOSITION channel (`close192` row B), because `jac`'s numerator is
          membership, which coverage cannot move;
  `reach` from the composition channel too, and it is the member this instrument is already
          known to get wrong on emission-law arms (audit 181; rule 13's toll).
"""
from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
INV = dict(reach=1.12927, fid=3.32482, jac=2.72081, nmae=1.64787)
SHRINK = 0.4277
V14 = dict(fid=-0.035470, jac=+0.004271, reach=+0.013731)
TOLL = -0.0146            # rule 13, emission-law class: mean(a21, a22) scaled `reach`


def main() -> None:
    cov = json.load(open(os.path.join(WORK, "cover192.json")))
    clo = json.load(open(os.path.join(WORK, "close192.json")))
    pri = json.load(open(os.path.join(WORK, "price192.json")))
    f = str(cov["best_f"])
    P = cov["prop"][f]["P"]
    pred = cov["prop"][f]["pred"]            # v14, v15, a8, a15, v11 order
    print("=" * 100)
    print("COMPOSED RETRODICTION OF v14 (falsifier 2).  One fitted parameter (f), three members.")
    d_fid_raw = pred[0] - pred[4]            # v14 minus v11, coverage channel, key-free
    rows = [("fid", d_fid_raw * INV["fid"], V14["fid"], "SIZE / coverage, key-free")]
    B = clo["rows"]["B tail-only  (top-1,563 by flat p-order OFF the |lfc_pred| top-850)"]
    v11 = clo["v11"]
    for m, lab in (("jac", "COMPOSITION, off our |lfc_pred| head"),
                   ("reach", "COMPOSITION (known-bad channel)")):
        d = (B["raw"][m] - v11[m]) * INV[m] * (SHRINK if m == "reach" else 1.0)
        rows.append((m, d, V14[m], lab))
    print(f"  {'member':<7}{'predicted scaled':>18}{'board scaled':>15}{'ratio':>8}{'sign':>7}"
          f"   channel")
    for m, p, b, lab in rows:
        print(f"  {m:<7}{p:>+18.6f}{b:>+15.6f}{p / b:>8.2f}"
              f"{'OK' if np.sign(p) == np.sign(b) else 'WRONG':>7}   {lab}")
    print(f"\n  `fid` and `jac` are BOTH retrodicted at 0.7x in size with the right sign, from "
          f"one\n  fitted parameter.  `reach` inverts -- the same member, on the same arm class, "
          f"that audit\n  181 found this instrument inverts on a21, and that rule 13's "
          f"{TOLL:+.4f} toll exists to charge.")

    print(f"\n{'=' * 100}")
    print("THE SIX-MEMBER LEDGER, both ways (rule 12).  Arm: call set = top-K by the fitted")
    print("gene-level membership key, realised through 191's emitter at bit-identical totals.")
    band = ["0.4", "0.5"]                     # the retrodicted band, from two board members
    for K in (1600, 2200, 2800, 3600, 4474):
        vals = [pri["rows"][f"{b}|{K}"] for b in band]
        dfid = np.mean([v["overall"]["fid"] for v in vals])
        djac = np.mean([v["overall"]["jac"] for v in vals])
        drch = np.mean([v["overall"]["reach"] for v in vals])
        tot = dfid + djac + drch
        print(f"\n  K = {K}   (mean over the retrodicted band f = {band})")
        print(f"    {'member':<8}{'d overall':>12}   provenance")
        print(f"    {'pds':<8}{0.0:>+12.6f}   PINNED — 191 falsifier 1, 900/900 delivered "
              f"per-gene panel totals digest-identical")
        print(f"    {'mse':<8}{0.0:>+12.6f}   PINNED — 185 §6, `claim` at the reference cap "
              f"0.24653 for every rc arm")
        print(f"    {'nmae':<8}{0.0:>+12.6f}   PINNED — 191 falsifier 2, +0.00011 raw h1-test, "
              f"0.564x challenge driver")
        print(f"    {'fid':<8}{dfid:>+12.6f}   this experiment, h1-test frontier at the "
              f"retrodicted n_real level")
        print(f"    {'jac':<8}{djac:>+12.6f}   this experiment")
        print(f"    {'reach':<8}{drch:>+12.6f}   this experiment, x0.4277 (164 §6 shrink); a "
              f"CHECK, not a target")
        print(f"    {'SUM':<8}{tot:>+12.6f}   WITHOUT rule 13's toll")
        print(f"    {'SUM+toll':<8}{tot + TOLL / 6:>+12.6f}   WITH rule 13's emission-law toll "
              f"({TOLL:+.4f} scaled reach)")
        print(f"    gate: shape-only +0.0019 -> "
              f"{'PASS' if tot + TOLL / 6 > 0.0019 else 'FAIL'} with toll, "
              f"{'PASS' if tot > 0.0019 else 'FAIL'} without")


if __name__ == "__main__":
    main()
