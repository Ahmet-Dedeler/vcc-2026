#!/usr/bin/env python3
"""604: unit tests of def_identity_prior's two reference guards on synthetic profiles (no data files read beyond CLEAN4).
  G1 same-call guard: a reference with the context's own identity call is never used.
  G2 duplicate-call guard: only one reference per identity call (bundle contexts first).
  G3 sister guard (626): a reference whose deviation-Spearman similarity to the context exceeds 0.3 is dropped.
  G4 the sister guard never leaves fewer than 2 references.
  G5 lam 0 still reproduces CLEAN4 exactly with the guards active."""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import def_identity_prior as D

with np.load(D.CLEAN4_PATH) as z:
    genes = list(map(str, z["genes"]))
P, ps = D.load_clean4(genes)
rng = np.random.default_rng(604)
G = len(genes)
base = rng.gamma(0.6, 30.0, G)                       # a shared expression backbone
def line(scale=0.6, rel=None, w=0.0):
    own = base * np.exp(rng.normal(0, scale, G))
    return own if rel is None else np.exp((1 - w) * np.log(own) + w * np.log(rel))
X = line(); Y = line(); Z = line(); H = line()
SIS = line(rel=X, w=0.85)                            # a sister of X (shares most of X's identity deviations)
ok = []
def rec(n, c):
    ok.append(c); print(("PASS " if c else "FAIL ") + n)

# G1 + G2: a bundle context X and a pinned reference with X's call; a duplicate-call pair Y / Y2
r = D.identity_prior(genes, {"X": X, "Y": Y}, {"X": "LX", "Y": "LY"}, {"pinX": X, "pinY": Y, "h1": H, "Z": Z},
                     {"pinX": "LX", "pinY": "LY", "h1": "H1_ESC", "Z": "LZ"}, 0.5, P, ps)
rx = r["X"]["prov"]["references"]
rec("G1 same-call reference never used (X's refs exclude pinX)", "pinX" not in rx and "X" not in rx)
rec("G2 one reference per call (Y from the bundle, not pinY)", "Y" in rx and "pinY" not in rx)
# G3: sister of X in the pool is dropped (pool large enough)
r = D.identity_prior(genes, {"X": X}, {"X": "LX"}, {"SIS": SIS, "Y": Y, "Z": Z, "h1": H},
                     {"SIS": "LSIS", "Y": "LY", "Z": "LZ", "h1": "H1_ESC"}, 0.5, P, ps)
pv = r["X"]["prov"]
rec(f"G3 sister dropped (sim {pv['sister_similarity'].get('SIS')})", "SIS" in pv["sister_dropped"] and "SIS" not in pv["references"])
# G4: only two references, one of them a sister -> keep both (never below 2)
r = D.identity_prior(genes, {"X": X}, {"X": "LX"}, {"SIS": SIS, "Y": Y}, {"SIS": "LSIS", "Y": "LY"}, 0.5, P, ps)
pv = r["X"]["prov"]
rec("G4 never fewer than 2 references", len(pv["references"]) == 2 and pv["sister_dropped"] == [])
# G5: lam 0 KA with guards active
r = D.identity_prior(genes, {"X": X}, {"X": "LX"}, {"SIS": SIS, "Y": Y, "Z": Z, "h1": H},
                     {"SIS": "LSIS", "Y": "LY", "Z": "LZ", "h1": "H1_ESC"}, 0.0, P, ps)
cov = np.isfinite(P)
rec("G5 lam 0 == CLEAN4 with guards", np.array_equal(r["X"]["psign"], ps) and np.array_equal(r["X"]["pdown"][cov], P[cov]))
print(f"{sum(ok)}/{len(ok)} PASS"); sys.exit(0 if all(ok) else 1)
