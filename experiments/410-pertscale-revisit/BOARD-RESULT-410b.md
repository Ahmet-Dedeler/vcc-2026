# 410b BOARD RESULT — entry `4WgVi6FIE4iVw9PmRUP3` ("a spreadsheet achieved enlightenment"), published 2026-09-22 19:23Z

**0.2365722 = +0.0126358 vs 395 (0.2239363) — NEW BEST, rank 15 at publication.** Predicted ≈ +0.0125 (BOARD-RESULT-410 §Consequence).
Package sha 51703467…, prediction 70ad6d81…, builder 40772dd1 (commit ebc1bed6). One variable vs 410: `VCC_PERT_SCALE_BALLAST=1`.

| member (raw) | 395 | 410 (welded) | 410b (ballast) | 410b − 395 | 410b − 410 |
|---|---:|---:|---:|---:|---:|
| pds | 0.854738 | 0.850914 | 0.851156 | −0.003582 | +0.000242 |
| **mse** | 0.809898 | 0.730213 | **0.730994** | **−0.078904** | +0.000780 |
| nmae | 0.937101 | 0.991476 | 0.936458 | −0.000643 | −0.055018 |
| fid | 0.514228 | 0.513441 | 0.514273 | +0.000045 | +0.000833 |
| reach | 0.283909 | 0.273245 | 0.283822 | −0.000087 | +0.010577 |
| jac | 0.040891 | 0.040911 | 0.040847 | −0.000044 | −0.000064 |

## Readings
- **Ballast unwelds the DE channel on the board.** All four DE members sit within 0.0007 raw of 395's; the welded 410's
  nmae +0.054 and reach −0.011 are gone. The mechanism (the pseudobulk sums raw counts, the DE table weighs cells 1/400)
  holds on the server scorer exactly as on the local one.
- **The bulk kept the whole gain.** mse −0.0789 raw vs 410's −0.0797: the ratio-mode + pair-jackknife-boost construction
  preserved the credit on A/B/C (single-ballast blocks at ~0.79× of the carriers' jackknife were enough given the 1.58–1.86×
  margin measured in jkmargin410.py).
- **pds −0.0036 raw remains** (410: −0.0038) — a property of the per-target scaling of the bulk (sum-to-one render), not of
  ballast. Recoverable only by a render that keeps the scaled direction; the target-gene mass sink failed on H1, open.
- **Structural consequence:** the bulk (pds, mse) and the DE channel (fid, reach, jac, nmae) are now separately designable.
  Every bulk arm previously closed because "nmae pays" (158, 224/228, 169, 187) is reopened under ballast.
