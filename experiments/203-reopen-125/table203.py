"""203 — aggregate the PLOLO arms with 162's pre-committed instrument and price the gate.

Instrument (162, fixed in advance, retrodicts 12 board points):
  * frame  = CHALLENGE (unweighted mean of contexts A/B/C, delivered);
  * folds  = unweighted mean over all four, no fold selected;
  * scale  = board `pds` is reference-scaled; 162's span is PDS_SPAN and its
             fitted transfer slope for `pds` is 1.011;
  * noise  = 162's leave-one-out sigma for the `pds` predictor, 0.00149 raw.

`mse` is priced through 096's `claim` route the way 197 inverted it: raw `mse` is an
affine function of the transfer cosine, `raw = 1.023611 - 2.074495 * cos`, verified out of
sample on `rc_hybrid` to 4.9e-7.  So a measured `d cos` maps to `d raw mse` directly, and
197's exchange rate `+0.00365 overall per +0.01 cosine` prices the whole direction channel.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
PDS_SPAN = 0.95228 - 0.50020          # 162 / 188
PDS_SLOPE = 1.011                     # 162's fitted h1->board transfer for pds
PDS_LOO_SIGMA = 0.00149               # 162, raw
COS_RATE = 0.00365 / 0.01             # 197: overall per unit cosine
MSE_SLOPE = -2.074495                 # 197: d raw mse / d cos

rows = json.loads((OUT / "arms203.json").read_text())
d = {(r["fold"], r["arm"]): r for r in rows}
arms = sorted({r["arm"] for r in rows})


def chal(fold, arm, key):
    s = d[(fold, arm)]["spaces"]
    return float(np.mean([s[f"chal{c}/delivered"][key] for c in "ABC"]))


def table(key):
    return {a: {f: chal(f, a, key) for f in FOLDS} for a in arms}


P, C, Q = table("pds"), table("cos"), table("sgn_w")
base = "a15"

print("=" * 118)
print("PLOLO, challenge frame (mean of A/B/C delivered), all four folds, unweighted mean "
      "over folds -- 162's pre-committed instrument")
print("=" * 118)
hdr = (f"{'arm':<30}{'alpha':>7}{'shrink':>8}" + "".join(f"{f[:11]:>12}" for f in FOLDS)
       + f"{'MEAN pds':>11}{'d pds':>10}{'d scaled':>10}")
print(hdr)
print("-" * len(hdr))
summary = {}
for a in arms:
    m = float(np.mean([P[a][f] for f in FOLDS]))
    dm = m - float(np.mean([P[base][f] for f in FOLDS]))
    r0 = d[(FOLDS[0], a)]
    print(f"{a:<30}{r0['alpha']:>7.3f}{r0['shrink_vs_a15']:>8.2f}"
          + "".join(f"{P[a][f]:>12.6f}" for f in FOLDS)
          + f"{m:>11.6f}{dm:>+10.6f}{dm/PDS_SPAN:>+10.6f}")
    summary[a] = dict(pds=m, d_pds=dm, d_pds_scaled=dm / PDS_SPAN)

print()
hdr2 = (f"{'arm':<30}" + "".join(f"{f[:11]:>12}" for f in FOLDS)
        + f"{'MEAN cos':>11}{'d cos':>10}{'MEAN q_w':>11}{'d q_w':>10}")
print(hdr2)
print("-" * len(hdr2))
for a in arms:
    mc = float(np.mean([C[a][f] for f in FOLDS]))
    mq = float(np.mean([Q[a][f] for f in FOLDS]))
    dc = mc - float(np.mean([C[base][f] for f in FOLDS]))
    dq = mq - float(np.mean([Q[base][f] for f in FOLDS]))
    print(f"{a:<30}" + "".join(f"{C[a][f]:>12.5f}" for f in FOLDS)
          + f"{mc:>11.5f}{dc:>+10.5f}{mq:>11.5f}{dq:>+10.5f}")
    summary[a].update(cos=mc, d_cos=dc, qw=mq, d_qw=dq)

print("\n" + "=" * 118)
print("GATE (rule 13, as amended 2026-09-03).  This arm CHANGES THE EFFECT MATRIX, so it is")
print("amplitude-changing and profile-changing: floor +0.0089, and `mse` is NOT exempt.")
print("=" * 118)
print(f"{'arm':<30}{'d pds/6*1.011':>15}{'d raw mse':>12}{'d mse/6':>10}"
      f"{'composed':>11}{'vs +0.0089':>12}{'cos route':>12}")
print("-" * 102)
for a in arms:
    s = summary[a]
    pds_term = s["d_pds_scaled"] / 6.0 * PDS_SLOPE
    d_raw_mse = MSE_SLOPE * s["d_cos"]
    # board `mse` is reference-scaled the same way `claim` inverts it; 096's route reports
    # the scaled member delta as d_raw / the reference denominator used in 197 (1.0 by
    # construction of the affine map, so the raw delta IS the member delta on that scale).
    mse_term = -d_raw_mse / 6.0        # lower raw mse is better
    composed = pds_term + mse_term
    cos_route = s["d_cos"] * COS_RATE
    s.update(pds_term=pds_term, d_raw_mse=d_raw_mse, mse_term=mse_term,
             composed=composed, cos_route=cos_route)
    print(f"{a:<30}{pds_term:>+15.6f}{d_raw_mse:>+12.6f}{mse_term:>+10.6f}"
          f"{composed:>+11.6f}{'PASS' if composed > 0.0089 else 'VETO':>12}"
          f"{cos_route:>+12.6f}")

print(f"\n162 LOO sigma on the `pds` predictor: {PDS_LOO_SIGMA:.5f} raw "
      f"= {PDS_LOO_SIGMA/PDS_SPAN:.5f} scaled = {PDS_LOO_SIGMA/PDS_SPAN/6:.6f} overall")
print("per-fold spread of d pds (the four folds are the only replicate we have):")
for a in arms:
    if a == base:
        continue
    v = np.array([P[a][f] - P[base][f] for f in FOLDS])
    print(f"  {a:<30} " + " ".join(f"{x:+.5f}" for x in v)
          + f"   mean {v.mean():+.6f}  sd {v.std(ddof=1):.6f}  "
          f"signs {int((v>0).sum())}/4 positive")

(OUT / "table203.json").write_text(json.dumps(summary, indent=1))
print(f"\nwrote {OUT/'table203.json'}")
