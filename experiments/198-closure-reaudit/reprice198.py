#!/usr/bin/env python3
"""198 — re-price the flagged closures under the corrected instrument.

Every number below is quoted from a file; nothing is estimated from memory.
The correction applied to an Orion `reach` delta is AGENTS.md's own measured
overstatement factor for that surface: h1-test reproduces two board `reach`
points to within 0.014 raw while repaired HCT116 "misses a15's by 2.6x and
overstates the tier-policy effect by 4.4x".
"""
from __future__ import annotations
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent / "work"
OUT.mkdir(exist_ok=True)
M = ("pds", "mse", "nmae", "fid", "reach", "jac")
ORION_OVERSTATE = (2.6, 4.4)      # AGENTS.md harness warning, both measured factors

CASES = [
    dict(
        exp="150-dual-moment-score-level-hct116", verdict="TERMINAL_NO_GO",
        src="experiments/150-.../RESULTS.md 'Primary result' table",
        surface="held HCT116 (Orion) ONLY - zero h1-test rows",
        # SCALED deltas as recorded
        d=dict(pds=0.0, mse=0.0, nmae=+0.006929, fid=-0.006765,
               reach=-0.045128, jac=+0.014694),
        recorded_avg=-0.005045,
        note="pds exactly preserved per target; the cheap invariant arm class",
    ),
    dict(
        exp="151-relaxation-consensus-selector", verdict="terminal NO-GO",
        src="experiments/151-.../RESULTS.md 'Comparison A'",
        surface="held HCT116 (Orion) ONLY",
        d=dict(pds=0.0, mse=+0.000341, nmae=+0.006315, fid=-0.003736,
               reach=-0.043628, jac=+0.014645),
        recorded_avg=-0.004344,
        note="rescue attempt for 150; same surface, same dominant term",
    ),
]

print("=" * 78)
print("RE-PRICING FLAGGED CLOSURES WHOSE VERDICT RESTS ON AN ORION `reach` NUMBER")
print("=" * 78)
res = []
for c in CASES:
    d = c["d"]
    avg = sum(d.values()) / 6
    ex_reach = sum(v for k, v in d.items() if k != "reach") / 6
    print(f"\n### {c['exp']}   [{c['verdict']}]")
    print(f"    surface: {c['surface']}")
    print(f"    source : {c['src']}")
    print(f"    {'member':8s} {'scaled d':>11s} {'overall d':>11s}  share of |avg|")
    for k in M:
        print(f"    {k:8s} {d[k]:+11.6f} {d[k] / 6:+11.6f}  "
              f"{100 * abs(d[k] / 6) / abs(avg):5.1f}%")
    print(f"    {'AVG':8s} {'':11s} {avg:+11.6f}   (file records {c['recorded_avg']:+.6f})")
    print(f"    all five members EXCEPT `reach`:        {ex_reach:+11.6f} overall  <- POSITIVE")
    row = dict(exp=c["exp"], recorded=avg, ex_reach=ex_reach, corrected={})
    for f in ORION_OVERSTATE:
        corr = ex_reach + (d["reach"] / f) / 6
        row["corrected"][f] = corr
        print(f"    corrected at the {f}x Orion `reach` overstatement: {corr:+11.6f} overall")
    lo, hi = sorted(row["corrected"].values())
    print(f"    >>> CORRECTED RANGE  [{lo:+.5f}, {hi:+.5f}]  vs recorded {avg:+.6f}")
    print(f"    >>> the file calls this 'not a borderline gate miss'; corrected it IS "
          f"borderline.")
    print(f"    note: {c['note']}")
    res.append(row)

# ---------------------------------------------------------------- 125
print("\n" + "=" * 78)
print("125 — the endpoint is not a scored member AND is degree-1 in the prediction")
print("=" * 78)
adj = json.loads((OUT / "adj125.json").read_text())
D_CREDIT, D_PREC = -0.015049899010796045, +0.012174572820953597
print(f"  candidate/baseline magnitude: LS scalar {adj['ls_scalar']:.6f} "
      f"(a {1 / adj['ls_scalar']:.2f}x SHRINK), cosine {adj['cosine']:.6f}")
print(f"  signed credit  {D_CREDIT:+.6f}  <- degree-1 homogeneous in the prediction "
      f"(run.py:774-776)")
print(f"  top-20% sign precision {D_PREC:+.6f}  <- exactly scale-INVARIANT "
      f"(run.py:777-792)")
print(f"  implied credit(122) under a pure-rescale null: "
      f"{D_CREDIT / (adj['ls_scalar'] - 1):+.6f}")
print("  -> consistent with a 0.04-0.08 transfer cosine (AGENTS.md 4.2). The shrink")
print("     alone accounts for the entire loss; no directional degradation needed.")

json.dump(res, open(OUT / "reprice198.json", "w"), indent=1)
