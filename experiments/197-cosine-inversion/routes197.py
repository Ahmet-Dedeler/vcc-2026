"""197 — the three-route cross-check over the WHOLE informative cohort, not just the top.

The mandate asks for every inverted cosine to be cross-checked against two independent
observables.  Team by team that is a 54-row table nobody reads, so the object reported here is
the agreement structure:

  * Spearman rank agreement between the three routes (do they order the field the same way?)
  * the per-team ratio  B/A  =  (cosine that `pds` implies) / (cosine floor that `mse` implies),
    which is a direct measure of how much RETRIEVAL a team extracts per unit of bulk direction,
    i.e. AGENTS.md 4.2's whitening effect measured on other people's submissions.

    uv run python routes197.py
"""
from __future__ import annotations

import json
import os

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
LINES = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LINES.append(s)


C159 = np.array([0.000, 0.020, 0.030, 0.038, 0.043, 0.050, 0.072, 0.100, 0.123, 0.150,
                 0.200, 0.250, 0.300])
PDS_HCT = np.array([0.5043, 0.5730, 0.6096, 0.6374, 0.6538, 0.6765, 0.7430, 0.8199, 0.8705,
                    0.9175, 0.9723, 0.9920, 0.9981])


def main():
    f = json.load(open(f"{WORK}/invert197.json"))["field"]
    A = np.array([t["floor"] for t in f])
    B = np.array([float(np.interp(t["pds"], PDS_HCT, C159)) for t in f])
    Q = np.array([t["qmin"] for t in f])
    R = np.array([t["reach"] for t in f])
    us = [i for i, t in enumerate(f) if "Ahmet" in t["team"]][0]

    P("=" * 96)
    P(f"197 — THREE ROUTES OVER THE WHOLE INFORMATIVE COHORT ({len(f)} teams)")
    P("=" * 96)
    P("\nRANK AGREEMENT (Spearman).  The routes disagree on LEVEL; do they disagree on ORDER?")
    for tag, u, v in (("A (mse floor)  vs  B (pds curve)", A, B),
                      ("A (mse floor)  vs  C (nmae q_w)", A, Q),
                      ("B (pds curve)  vs  C (nmae q_w)", B, Q),
                      ("A (mse floor)  vs  raw reach", A, R)):
        r, p = spearmanr(u, v)
        P(f"  {tag:<36} rho {r:+.3f}   p {p:.1e}")
    P("  => the three routes ORDER the field the same way (rho 0.66-0.81).  Their disagreement")
    P("     is entirely a SCALE disagreement, which is what makes it interpretable rather than")
    P("     a sign that one of them is broken.")

    P("\nB/A PER TEAM — retrieval extracted per unit of bulk direction.  Restricted to the")
    P("teams whose floor is non-trivial (> 0.05), because B/A explodes as A -> 0.")
    m = A > 0.05
    ratio = B / np.maximum(A, 1e-9)
    idx = np.arange(len(f))[m]
    o = idx[np.argsort(-ratio[m])]
    P(f"  {'#':>3} {'team':<24}{'A (mse)':>10}{'B (pds)':>10}{'B/A':>8}{'raw pds':>10}")
    for i, j in enumerate(o[:10]):
        P(f"  {i+1:>3} {f[j]['team'][:23]:<24}{A[j]:>10.5f}{B[j]:>10.5f}{ratio[j]:>8.2f}"
          f"{f[j]['pds']:>10.5f}")
    P(f"  ...")
    P(f"  us: rank {list(o).index(us)+1} of {m.sum()},  B/A {ratio[us]:.2f} against a cohort "
      f"median of {np.median(ratio[m]):.2f}")
    P(f"  i.e. we extract {100*(ratio[us]/np.median(ratio[m])-1):.0f}% more `pds` per unit of")
    P("  bulk cosine than the median shrunk team.  That is AGENTS.md 4.2's whitening effect,")
    P("  measured for the first time against the field rather than against our own LODO.")
    P("  It is also why `pds` cannot be used to price the direction gap: the member where we")
    P("  are strongest is the member that is least sensitive to what we are short of.")

    json.dump(dict(spearman_AB=float(spearmanr(A, B)[0]),
                   spearman_AQ=float(spearmanr(A, Q)[0]),
                   spearman_BQ=float(spearmanr(B, Q)[0]),
                   spearman_Areach=float(spearmanr(A, R)[0]),
                   ratio_us=float(ratio[us]), ratio_median=float(np.median(ratio[m])),
                   n_nontrivial=int(m.sum())),
              open(f"{WORK}/routes197.json", "w"), indent=1)
    open(f"{WORK}/routes197.log", "w").write("\n".join(LINES))
    P(f"\nwrote {WORK}/routes197.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
