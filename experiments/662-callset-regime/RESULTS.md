# 662 — Which truth-power level does the board match for CALL-SET changes? (after K20 failed 610's card)

[21dd4c], 2026-10-04, coordinator's ask after the K20 read (59543bfe). Analysis only: no emission, no cell-eval2, no
`vcc`, every process < 1 GB (peak RSS 0.75 GB).

**Answer: f ≈ 0.61 (0.60–0.62). It is neither f 0.35 nor f 1.**
- At f 0.61 ONE fid multiplier (×2.69) fits both board call-set pairs: cp (composition) and K20 (size).
- Their jac multipliers agree there (0.43 / 0.49).
- Every priced sign of all four board DE pairs is right there: cp fid / jac, K20 fid / jac, zs nmae, coarse-MDE nmae.
- f 1 matches K20's fid SIGN only:
  - it over-reads K20's fid 6× (H1 −0.0100 vs board −0.0017);
  - it mis-signs K20's jac (−0.0029 vs +0.0002);
  - it needs fid multipliers 17× apart for the two pairs (cp 2.83, K20 0.17).
- f ≤ 0.5 (611's f 0.35 included) mis-signs K20's fid.
- **Corrected instrument I-662:** price call-set changes on H1 at f 0.61 with fid × 2.7 and jac × 0.43–0.49.
  - Registered tests: K16 on 10-06 (`PRED-662-K16.md`: −0.0067 overall vs the f 1 reading's −0.0030) and CPH on 10-07
    (`PRED-662-CPH.md`).

## What was run (all read-only re-scoring of runs on disk, each with a known answer)
| script | runs | grid | known answer |
|---|---|---|---|
| `regime662.py` | 537's H1 real-emitter comparisons (cp, K20, K16, the CLEAN4 hook P1 / P8 / P9 on the cp base and on 446bj2z, CPH totals); 3 seeds; 528's rows_at on lib521's state | PowerSurface f 0.35 … 1.0 in 0.05 steps, 0.01 steps over 0.55–0.70, + cmp | 48 (comparison, level) cells == `537/work/price537.json` (max err 0) |
| `delivery662.py` | 611's two delivery pairs on h1-test (zs: 435 coarse vs coarse+zs, 3 seeds; coarse MDE: 426 vs 435, seed 0) | same | == `610/work/desk611b.json` (max err 0) |
| `tv662.py` | 619's three pairs on 618's h1-trainval surface (zs, coarse MDE, cp [worn]) | same | == `619/work/score619.json` (max err 0) |
| `t496_662.py` | 496a's own H1 pair (τ 0.05 hook on the plain coarse base vs 435) | same | members == each run's metrics json at f 1 (2e-16); Δnmae −0.001458 = 496's −0.00146 |
| `fit662.py` | the fit, the delivery check, the K16 / CPH prices | | — |

PowerSurface level f: each target's truth pool is the reference's top round(f × n_real) genes by p (lib521.pools). fid
reads f only through n_real in max(n_pred, n_real). That is why a call-set SIZE change flips sign with f: shrinking K
raises sign precision on targets with n_real < K, and loses yield on targets with n_real > K.
- On h1-test, 15 % of targets (13 / 85) have n_real ≥ 2,000 at f 1;
- 3.5 % at f 0.5;
- 5.9 % at f 0.6.

## The fit (board raw Δ vs 446bj2z: cp fid −0.00358, jac +0.00399; K20 fid −0.00169, jac +0.000185, nmae −0.00006)
| f | cp fid (H1; board/H1) | cp jac board/H1 | K20 fid (H1; board/H1) | K20 jac (H1; board/H1) | K20 nmae H1 | delivery nmae board/H1: h1-test zs / MDE; 618 zs / MDE | priced signs missed |
|---|---|---|---|---|---|---|---|
| 0.35 | −0.00183; 1.95 | 0.48 | +0.00348; **−0.48** | +0.00241; 0.08 | −0.00218 | 1.18 / 1.59; 1.17 / 1.04 | 1 (K20 fid) |
| 0.50 | −0.00176; 2.03 | 0.45 | +0.00088; **−1.92** | +0.00131; 0.14 | −0.00028 | 1.20 / 1.52; 1.98 / 0.98 | 1 |
| 0.55 | −0.00161; 2.22 | 0.44 | −0.00006; 26 | +0.00091; 0.20 | −0.00027 | 1.27 / 1.54; 1.48 / 0.99 | 0 |
| 0.60 | −0.00144; 2.49 | 0.43 | −0.00048; 3.51 | +0.00046; 0.41 | −0.00050 | 0.91 / 1.34; 1.29 / 1.03 | 0 |
| **0.61** | −0.00133; **2.70** | **0.43** | −0.00063; **2.69** | +0.00038; **0.49** | −0.00036 | 0.92 / 1.37; 1.30 / 1.05 | **0** |
| 0.62 | −0.00122; 2.93 | 0.43 | −0.00078; 2.17 | +0.00026; 0.70 | −0.00025 | 0.95 / 1.39; 1.32 / 1.05 | 0 |
| 0.65 | −0.00114; 3.15 | 0.43 | −0.00134; 1.27 | +0.00002; 12 | −0.00018 | 0.90 / 1.39; 1.40 / 1.08 | 0 |
| 0.70 | −0.00118; 3.03 | 0.44 | −0.00290; 0.58 | −0.00049; **−0.37** | +0.00018 | 0.89 / 1.47; 1.53 / 1.13 | 1 (K20 jac) |
| 1.0 | −0.00126; 2.83 | 0.44 | −0.00997; 0.17 | −0.00286; **−0.06** | +0.00034 | 1.19 / 1.74; 1.34 / 1.14 | 1 |

- **Common-fid-multiplier level:** board_cp / H1_cp(f) = board_K20 / H1_K20(f), linear in f between grid points.
  - The mean over seeds gives one root: **f\* = 0.610, m_fid = 2.69**.
  - Per seed: f 0.614 / 0.601 / 0.614, m 2.95 / 2.62 / 2.54.
- **Out-of-fit check: jac.**
  - Three parameters (f, m_fid, m_jac) were fitted to four numbers (two pairs × fid and jac); jac is the one left over.
  - At f\*, cp's jac ratio is 0.43 and K20's is 0.49. With m_jac 0.43, I-662 retrodicts K20's jac as +0.00016 against
    the board's +0.00019.
  - K20's nmae at f\* is −0.00036 vs the board's −0.00006: 0.0001 overall, tiny. K20 is not an nmae lever.
- **Delivery pairs at f\*: nmae ratios 0.92–1.37 on both panels.** At f 0.35 they are 1.04–1.59, at f 1 1.14–1.74. So f\* is
  no worse for delivery levers. Their reach stays mis-signed at f\*, as 619 found at f 0.35 on 618: reach stays unpriced.
- **618 (second panel, cp only, worn):** at f\*, cp fid 0.81×, jac 0.49×.
  - The jac multiplier holds across panels and across pair types (0.43–0.49).
  - The fid multiplier is panel-specific: 618 0.8× vs h1-test 2.7×. That was already 619's finding; only the level changes.
- **192 §4.2's earlier estimate** was f 0.35–0.50 (fid coverage best at f 0.5, RMS 0.00126 vs 0.00427 at f 0.7; jac level
  f ≈ 0.35). I-662's 0.61 sits just above that band. K20 is the first board pair that reads the n_real tail around 2,000
  directly.

## I-662 (the corrected instrument for call-set changes)
- **H1 level:** H1 real emitter at f 0.61 (band 0.60–0.62).
- **fid × 2.7** on h1-test: one value for composition (cp) and size (K20). On 618 × 0.8 (cp, worn).
- **jac × 0.43–0.49** on both panels.
- **nmae** for delivery levers × 0.92–1.37. For a hook, use T from 496a's own board read at f\* (−0.001979 on H1).
- **reach:** unpriced.
- **Status:** fitted on two pairs, with one out-of-fit check (K20's jac) passing. K16 (10-06) is its first out-of-sample
  test.
- **The dial:** K16's Δfid maps to a level through K20's anchor: −0.0130 → f 0.60; −0.0116 → 0.61; −0.0084 → 0.65;
  −0.0053 → 0.70; −0.0031 → f 1.

## Consequences
1. **611's rule and 619's call-set line are amended** (in their own files):
   - call-set changes are priced at f 0.61, not f 0.35;
   - fid × 0.8–2.7 by panel, replacing × 0.8–2 at f 0.35;
   - 611's ×2 / ×0.5 survive roughly (×2.7 / ×0.43–0.49), but at the new level.
2. **Every f 0.35 call-set price is void.** That includes the K16 card, K16 / K20 in 537's tree and CPH's f 0.35 band.
   CPH's price barely moves (its members keep their signs at every level). K-size prices flip.
3. **October K stays 2,200.** I-662 says K 1,600 / 2,000 lose on the board, as the K20 read did.
   - For D/E/F, a K change could only be priced from D/E/F's own n_real tail, which is unknown until D/E/F reads exist.
   - Class: instrument calibration (not a recipe change). Tries: one fitted level on two board pairs; no board slot spent
     by 662.
4. **Nmae level for magnitude levers: registered for the reads, not settled here** (12a6809b, before each read).
   - The delivery pairs cannot tell f 0.35 from f\* for nmae, so 661's f 0.35 lines stand.
   - An I-662 line was added beside each (`661/alt662.py`; known answers vs stacks618b / score618 / score661oct exact):
     - 630: +0.0035 vs +0.0030;
     - **DNMAG: +0.00065 vs +0.0015** (DN's increment is −0.0049 at f 0.35 but −0.0020 at f\*: the discriminating read);
     - idall3: +0.0002, a borderline TIE / WIN (vs TIE at f 0.35);
     - 625's rule on 618 passes {4} (vs {} at f 0.35).
   - `PRED-662-496a.md` prices 10-05's plain-hook read at f\*: h1-test −0.0022 (−0.0027 … −0.0016), 618 −0.0014
     (−0.0017 … −0.0011). Both say TRANSFERS. That read also sets CPH's T.

## Files
- **Scripts:** `regime662.py`, `delivery662.py`, `tv662.py`, `t496_662.py`, `fit662.py`.
- **Outputs:** `work/*.json`, `work/*.out`.
- **Predictions:** `PRED-662-K16.md`, `PRED-662-CPH.md`, `PRED-662-496a.md`; the I-662 lines in `661/PRED-661-*.md`.
