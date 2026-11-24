"""169 — the deliverable: both arms, all six members, rule 13's gate, from the saved JSONs.

`gate169.py` did this for the `h`-space arm alone and read `nmae` from the DETERMINISTIC
profile-based evaluator.  Two corrections are folded in here:

  1. the `lfc`-space arm (`sharpen_lfc`, builder-norm restored) is carried alongside;
  2. `nmae` on h1-test is taken from the EMITTED tables -- the scorer's own path, five matched
     seeds, sem 0.00006 -- rather than from the profile-based evaluator.  The two differ by ~2x
     on this surface because the profile-based form uses the UNFLOORED profile while
     `emit_counts` floors `n*L*profile`, and the DE gate here is dominated by low-count genes
     where that bias is largest.  098's `validate.py` checked the analytic form against the
     realised one to +0.0003 on 40 rows of a12; it does not hold at this gate depth on this
     artifact, so the emitted number is the one used and the analytic ladder
     (`nmaegate169.py`) is kept only for the gate-DEPTH comparison, where the bias is common
     to both arms of every ratio.

Nothing is re-derived by hand.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, HERE)
from debench import ANCHOR, Surface                                    # noqa: E402
from six169 import overlap_ratio, precision                            # noqa: E402
from reach169 import K_REACH, per_pert_reach                           # noqa: E402

SPAN = {k: ANCHOR[k][1] - ANCHOR[k][0] for k in ANCHOR}
SPAN_NMAE = 0.60                      # conservative; our own board points solve to 0.6674
TOLL = -0.014569 / 6
K_LIT, K_FIT = 1.011, 0.899
SCORE = lambda r: min(max((0.987862879 - r) / 0.950162566, 0.0), 1.0)   # noqa: E731
NR, NP = 240.0, 4424.0                # the board's own counts (164 §3.3), for the `jac` rebuild

ARMS = {
    "h-space   (sharpen)":      dict(pspace="chal%s", claim="claim169.json", tag="_mf",
                                     nmaekey="hct116|h"),
    "lfc-space (sharpen_lfc)":  dict(pspace="chal%sL", claim="claim169-lfc.json", tag="L_mf",
                                     nmaekey="hct116|lfc"),
}
BASE = ["q1.00_mf", "q1.00_mf_s2", "q1.00_mf_s3", "q1.00_mf_s4", "q1.00_mf_s5"]
QS = (1.05, 1.10, 1.15, 1.20, 1.25, 1.30)


def de_members(names):
    """Seed-averaged h1-test members for a list of matched-seed arms."""
    reach, prec, rec, nmae = [], [], [], []
    for n in names:
        p = os.path.join(HERE, f"de169.{n}.npz")
        reach.append(per_pert_reach(p, 1.0).mean())
        s = Surface(p)
        prec.append(precision(s))
        rec.append(overlap_ratio(s))
        nmae.append(s.nmae())
    return dict(reach=np.array(reach), precision=np.array(prec),
                recall=np.array(rec), nmae=np.array(nmae))


def main() -> None:
    P = json.loads(open(os.path.join(WORK, "plolo169.json")).read())
    R = {(r["fold"], r["arm"], r["q"], r["space"]): r for r in P}
    folds = sorted({r["fold"] for r in P})
    NG = json.loads(open(os.path.join(WORK, "nmaegate169.json")).read())
    bs = np.load(f"{REPO}/experiments/096-mse-floor/work/board_solve.npz")
    sD, Kc, X0 = float(bs["sumD"]), float(bs["v11_Kc"]), float(bs["X"])
    s15 = 0.982 / 1.0796

    b_de = de_members([n for n in BASE
                       if os.path.exists(os.path.join(HERE, f"de169.{n}.npz"))])
    nb = len(b_de["reach"])

    for label, cfg in ARMS.items():
        C = json.loads(open(os.path.join(WORK, cfg["claim"])).read())

        def dpds(q):
            v = [np.mean([R[(f, "a23", q, cfg["pspace"] % c + "/delivered")]["pds"]
                          for c in "ABC"])
                 - np.mean([R[(f, "a23", 1.00, cfg["pspace"] % c + "/delivered")]["pds"]
                            for c in "ABC"]) for f in folds]
            return float(np.mean(v)) / SPAN["pds"], float(np.std(v)) / SPAN["pds"]

        def dmse(q):
            def raw(qq):
                c = C[repr(float(qq))]
                return 1 + (s15**2 * c["A"] - 2 * s15 * X0 + Kc
                            - min(s15**2 * c["B"], c["claim"])) / sD
            return (SCORE(raw(q)) - SCORE(raw(1.00))) / 6

        print(f"\n{'=' * 132}\n{label}   (all values OVERALL board units vs a15 + a23)")
        print(f"  {'q':>5} {'pds PLOLO':>11} {'scaled':>10} {'pds/6':>10} {'mse':>10} "
              f"{'reach':>10} {'fid':>10} {'jac':>10} {'nmae h1':>10} {'nmae hct':>10} "
              f"{'SUM(h1)':>10} {'SUM(hct)':>10} {'GATE':>6}")
        for q in QS:
            names = [f"q{q:.2f}{cfg['tag']}"] + [f"q{q:.2f}{cfg['tag']}_s{i}"
                                                 for i in (2, 3, 4, 5)]
            names = [n for n in names if os.path.exists(os.path.join(HERE, f"de169.{n}.npz"))]
            # legacy seed-2 name for the h-space q=1.15 twin
            if not names and q == 1.15 and cfg["tag"] == "_mf":
                names = ["q1.15_mf", "q1.15_mf_b"]
            if q == 1.15 and cfg["tag"] == "_mf":
                names = ["q1.15_mf", "q1.15_mf_b", "q1.15_mf_s3", "q1.15_mf_s4", "q1.15_mf_s5"]
            if q == 1.25 and cfg["tag"] == "_mf":
                names = ["q1.25_mf", "q1.25_mf_s2", "q1.25_mf_s3"]
            m, sd = dpds(q)
            pds6 = m * K_LIT / 6
            mse6 = dmse(q)
            hct = NG.get(f"hct116|{cfg['nmaekey'].split('|')[1]}|0.35|{q}")
            nmae_hct = -hct / SPAN_NMAE / 6 if hct is not None else float("nan")
            if names:
                a = de_members(names)
                k = min(len(names), nb)
                dr = float(np.mean(a["reach"][:k] - b_de["reach"][:k]))
                dp = float(np.mean(a["precision"][:k] - b_de["precision"][:k]))
                dq = float(np.mean(a["recall"][:k] - b_de["recall"][:k]))
                dn = float(np.mean(a["nmae"][:k] - b_de["nmae"][:k]))
                reach6 = dr * K_REACH / SPAN["reach"] / 6
                fid6 = dp / SPAN["fid"] / 6
                r0 = float(np.mean(b_de["recall"][:k]))
                j0 = r0 * NR / (NP + NR - r0 * NR)
                j1 = (r0 + dq) * NR / (NP + NR - (r0 + dq) * NR)
                jac6 = (j1 - j0) / SPAN["jac"] / 6
                nmae_h1 = -dn / SPAN_NMAE / 6
            else:
                reach6 = fid6 = jac6 = nmae_h1 = float("nan")
            z = lambda x: 0.0 if np.isnan(x) else x                      # noqa: E731
            s_h1 = pds6 + mse6 + z(reach6) + z(fid6) + z(jac6) + z(nmae_h1)
            s_hc = pds6 + mse6 + z(reach6) + z(fid6) + z(jac6) + z(nmae_hct)
            print(f"  {q:>5.2f} {m:>11.6f} {sd:>10.6f} {pds6:>+10.6f} {mse6:>+10.6f} "
                  f"{reach6:>+10.6f} {fid6:>+10.6f} {jac6:>+10.6f} {nmae_h1:>+10.6f} "
                  f"{nmae_hct:>+10.6f} {s_h1:>+10.6f} {s_hc:>+10.6f} "
                  f"{'PASS' if s_h1 > 0.0019 else 'FAIL':>6}")

        m15, _ = dpds(1.15)
        print(f"\n  rule 13 at q = 1.15:")
        print(f"    Gate 1  |0.899 * d pds| = {abs(K_FIT * m15):.5f}  vs 0.00246  "
              f"{'PASS' if abs(K_FIT * m15) >= 0.00246 else 'FAIL'}")
        print(f"    Gate 4  d pds = {m15:+.6f} scaled  vs the +0.0126 shape-only budget  "
              f"({m15 / 0.0126:.2f}x)  {'PASS' if m15 > 0.0126 else 'FAIL'}")
        print(f"    literal (blanket reach toll, k=1.011)         "
              f"{m15 * K_LIT / 6 + TOLL:+.6f}  "
              f"{'PASS' if m15 * K_LIT / 6 + TOLL > 0.0019 else 'FAIL'}")


if __name__ == "__main__":
    sys.exit(main())
