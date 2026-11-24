"""438: re-level factors from the FOLD-ROUTE energies, for a panel with no 3-donor operator (efr).

The fold route (each design's realised render energy averaged over the three cancer folds, anchored at 417a's measured
energy through A/B/C's own fold-route 417a) reproduces the full-panel prices on abc and real -- the two panels where
every fold carries all targets, as efr's do -- to <= 0.0011 levelled (surr438). ntclevel's rule is then applied as in
full438: s_c = sqrt(rho* claim_c / B_c), rho* = exp(2 sd Phi^-1(0.148)), sd 0.01.

    uv run python foldlevel438.py <panel> <claims_full_tag>      -> work/<panel>/level438.json
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402

import numpy as np  # noqa: E402

import price438 as P8  # noqa: E402

sys.path.insert(0, str(P8.HERE))


def main(panel, claims_tag, sd=0.01):
    from scipy.stats import norm
    claims = json.loads((P8.WORK / claims_tag / "full438.json").read_text())["claims"]
    E = P8.fold_energies(panel, "held", claims)
    rho_t = float(np.exp(2.0 * sd * norm.ppf(0.31 / (0.31 + 1.78))))
    fac = {n: {c: float(np.sqrt(rho_t * claims[c] / e["AB"][c]["B"])) for c in P8.SLOTS} for n, e in E.items()
           if not n.startswith("CTRL:")}
    out = dict(panel=panel, route="fold", claims_from=claims_tag, sd=sd, rho_target=rho_t, claims=claims, factors=fac)
    (P8.WORK / panel / "level438.json").write_text(json.dumps(out, indent=1))
    print(f"{panel}: {len(fac)} designs; 417a analogue factors " + "/".join(f"{fac[P8.REF][c]:.4f}" for c in P8.SLOTS)
          + "; stack " + "/".join(f"{fac['h2c1.5|K3+sh1.25'][c]:.4f}" for c in P8.SLOTS))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
