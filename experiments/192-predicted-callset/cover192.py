"""192 — the KEY-FREE retrodiction: fitting the board's own `n_real` distribution from the
coverage term alone.

`close192.py` found that the size-only retrodiction of v14's `fid` flips sign depending on
which p-order is used to define the truncation (a15's tiered order: -0.0366 scaled; a8's flat
order: +0.0062).  That is an unmodelled variable and it disqualifies a key-dependent fit.

But `fid` factorises:

    fid = mean_p  k_p / max(n_pred, n_real_p)  =  mean_p prec_p * min(1, n_pred / n_real_p)

and 185 §4 measured that PRECISION IS 0.507-0.517 FOR EVERY KEY AT EVERY K on h1-test -- no
key selects a high-precision subset.  So the whole SHAPE of `fid(n_pred)` is the coverage term,
which depends on the `n_real` distribution and on NOTHING ELSE.  Fitting

    fid_hat(n) = P * mean_p min(1, n / n_real_p)          [ P = the saturated board precision ]

to the board's five same-profile points is therefore key-free, and it MEASURES the quantity
rule 9 says our harnesses cannot read.

Two shape families are fitted, because the proportional family is an assumption:
  * PROPORTIONAL   n_real_p = round(f * n_real_p^{h1})            (the mandate's instrument)
  * TAIL-CAPPED    n_real_p = min(cap, n_real_p^{h1})             (a21's light-tail evidence)
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E185 = os.path.join(REPO, "experiments/185-rank-controlled-emission")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, HERE)
from debench import Surface  # noqa: E402

BOARD = [("v14", 1563, 0.49126), ("v15", 1911, 0.49654), ("a8", 3836, 0.5021),
         ("a15", 4424, 0.5022), ("v11", 5961, 0.50178)]


def cover(nr, n):
    return float(np.mean(np.minimum(1.0, n / np.maximum(nr, 1))))


def main() -> None:
    s = Surface(os.path.join(E185, "work/de185.a15.npz"))
    z = np.load(os.path.join(E185, "work/de185.a15.npz"), allow_pickle=True)
    rp = np.where(np.isfinite(z["real_p"]), z["real_p"], np.inf)
    base = ((z["real_padj"] < 0.05) & np.isfinite(z["real_padj"])).sum(1)
    G = rp.shape[1]
    nr0 = np.maximum(base, 1)
    obs = np.array([v for _, _, v in BOARD])
    ns = np.array([n for _, n, _ in BOARD], float)

    def fit(nr):
        c = np.array([cover(nr, n) for n in ns])
        P = float((obs * c).sum() / (c * c).sum())          # least-squares precision
        pred = P * c
        return P, pred, float(np.sqrt(np.mean((pred - obs) ** 2)))

    print("=" * 100)
    print("KEY-FREE FIT of the board's n_real distribution from the coverage term of `fid`.")
    print("  Board:  " + "  ".join(f"{a}(n={n}) {v:.5f}" for a, n, v in BOARD))
    print(f"\n  PROPORTIONAL family  n_real_p = round(f * h1-test n_real_p)")
    print(f"  {'f':>6}{'median':>8}{'mean':>8}{'P>1563':>8}{'P (fitted prec)':>17}"
          f"{'RMS':>10}   predicted fid")
    outp = {}
    for f in (0.15, 0.25, 0.35, 0.5, 0.7, 0.85, 1.0, 1.3, 1.7, 2.2, 3.0):
        nr = np.maximum(np.round(f * nr0).astype(int), 1)
        P, pred, rms = fit(nr)
        outp[f] = dict(P=P, rms=rms, pred=pred.tolist(), med=float(np.median(nr)),
                       mean=float(nr.mean()))
        print(f"  {f:>6}{np.median(nr):>8.0f}{nr.mean():>8.0f}{(nr > 1563).mean():>8.3f}"
              f"{P:>17.5f}{rms:>10.6f}   " + " ".join(f"{x:.5f}" for x in pred))
    bp = min(outp, key=lambda k: outp[k]["rms"])
    print(f"  BEST proportional f = {bp}  (RMS {outp[bp]['rms']:.6f}, fitted precision "
          f"{outp[bp]['P']:.5f}; median n_real {outp[bp]['med']:.0f}, mean "
          f"{outp[bp]['mean']:.0f})")

    print(f"\n  TAIL-CAPPED family  n_real_p = min(cap, h1-test n_real_p)")
    print(f"  {'cap':>6}{'median':>8}{'mean':>8}{'P>1563':>8}{'P (fitted prec)':>17}"
          f"{'RMS':>10}   predicted fid")
    outc = {}
    for cap in (200, 400, 700, 1200, 2000, 3500, 6000, 10779):
        nr = np.minimum(nr0, cap)
        P, pred, rms = fit(nr)
        outc[cap] = dict(P=P, rms=rms, pred=pred.tolist(), med=float(np.median(nr)),
                         mean=float(nr.mean()))
        print(f"  {cap:>6}{np.median(nr):>8.0f}{nr.mean():>8.0f}{(nr > 1563).mean():>8.3f}"
              f"{P:>17.5f}{rms:>10.6f}   " + " ".join(f"{x:.5f}" for x in pred))
    bc = min(outc, key=lambda k: outc[k]["rms"])
    print(f"  BEST cap = {bc}  (RMS {outc[bc]['rms']:.6f}, fitted precision "
          f"{outc[bc]['P']:.5f}; median n_real {outc[bc]['med']:.0f}, mean "
          f"{outc[bc]['mean']:.0f})")

    # a21's independent constraint, stated but NOT fitted (its emission law also moved)
    print(f"\n  a21's independent check (n_pred 1,141; board raw fid 0.494779, and its emission")
    print(f"  law separately moved precision +0.0072 on h1-test, 185's a21_like):")
    for lab, nr in (("proportional f=" + str(bp), np.maximum(np.round(bp * nr0).astype(int), 1)),
                    ("cap=" + str(bc), np.minimum(nr0, bc))):
        P = outp[bp]["P"] if lab.startswith("prop") else outc[bc]["P"]
        pr = (P + 0.0072) * cover(nr, 1141)
        print(f"    {lab:<22} predicted a21 fid {pr:.5f}   observed 0.49478   "
              f"d {pr - 0.49478:+.5f}")

    with open(os.path.join(WORK, "cover192.json"), "w") as fh:
        json.dump(dict(prop={str(k): v for k, v in outp.items()},
                       cap={str(k): v for k, v in outc.items()},
                       best_f=bp, best_cap=bc), fh, indent=1, default=float)
    print(f"\nwrote {WORK}/cover192.json")


if __name__ == "__main__":
    main()
