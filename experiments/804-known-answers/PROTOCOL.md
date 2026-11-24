# 804 — the checker's known-answer retrodiction: protocol (frozen before any checker number exists)

Written 2026-10-15 ~22:50Z (committed fe4dba9d, 22:53:39Z) by the lab lead ([c61d4c]). AGENTS 11.1 rule 3: the checker is trusted only after it
retrodicts known board results. This file fixes HOW it is graded before the engines (801 bulk, 803 DE) return anything.
It is amended only by a dated section at the bottom, written before the numbers it affects.

## The rows (board side: `work/board_pairs.json`, from docs/BOARD-LEDGER.md by `board_pairs804.py`)
| row | board Δoverall | class | checker arms (cand vs base) | surfaces |
|---|---|---|---|---|
| 495 | +0.0016341 | bulk-key | cap6 vs cap3 (both with the uniform q 1.75 of 446bj2) | BULK-ABC, BULK-W |
| 479p2g | +0.0010365 | bulk-shape | cap3+P2g vs cap3 | BULK-ABC, BULK-W; DE-H1 if expressible |
| 498p2g | +0.0010083 | bulk-shape | pct+P2g vs pct | BULK-ABC, BULK-W; DE-H1 if expressible |
| 479sp2 | −0.0007141 | bulk-shape | cap3+S_P2 vs cap3 | BULK-ABC, BULK-W; DE-H1 if expressible |
| 702 | −0.0025464 | bulk-external | pct+P2g + PIE blend (w 0.25, PIE-seen targets, norm-preserving, on the fold line's PIE rows) vs pct+P2g | the K562 folds only (PIE trained on Orion: Orion folds are seen surfaces and may only bound a cost) |
| 446bj2z (zs) | +0.0005494 | de-mag | zs vs none | DE-H1 (h1test, h1tv) |
| 496a (hook) | +0.0005552 | de-sign | c4h vs none | DE-H1 |
| DN | +0.0003736 | de-mag | 630+DN vs 630 | DE-H1 |
| extras | 630arm +0.0009296 (de-mag), idall3 vs hook −0.000501 (de-sign), K20 −0.0007163 / K16 −0.0021068 (de-callset), 471c −0.000379 (bulk-union), 460dr +0.0028693 (bulk-key) | | where the engines express them | |

## How the checker is graded
1. **Units.** The checker's price is `check800.py aggregate` with default settings: board INV (DESIGN §5); bulk members
   pds (literal) + mse (518's identity, k 1); DE members fid + jac + nmae at regime f 0.35 (619: board-scale for nmae;
   reach reported, unpriced); paired target bootstrap 4,000 draws, seed 800; best-of-N with N = 1 (known answers are
   single configurations). Bulk and DE blocks are priced separately; a bulk arm's DE side effects are reported from DE-H1
   when the arm is expressible there.
2. **Decision.** SHIP iff the 2.5 % bound > 0; HARM iff the 97.5 % bound < 0 or point ≤ −0.0005; else WEAK / NULL by sign.
3. **Hit.** A row is a hit when the checker's decision matches the board's: board Δ ≥ +0.0003 → the checker must SHIP;
   board Δ ≤ −0.0003 → the checker must NOT ship (HARM, NULL or WEAK all count; HARM is reported separately as a strong
   hit). |board Δ| < 0.0003 → not graded (inside the class noise).
4. **Over-read ratio.** ρ = checker point / board Δ for graded rows with matching sign, per class and per member (pds and
   mse separately for bulk rows: the board's scaled member deltas are in board_pairs.json). Reported as median, range and
   leave-one-out range; never fitted to make a row pass.
5. **Surfaces are read separately and pooled.** BULK-ABC (the board's panel, worn), BULK-W (fresh lab panels, 18 folds),
   DE-H1 (h1test, h1tv). The headline hit rate uses the surface set the checker will use for new candidates: BULK-W for
   bulk rows and DE-H1 for DE rows. BULK-ABC is reported beside it.
6. **Pass criterion (DESIGN §9):** on the headline surfaces the checker does not ship 702 or 479sp2, and ships at least 4
   of the 5 positive core rows (495, 498p2g, zs, hook, DN; 479p2g counted as an extra). Otherwise the report names the
   class it cannot price, and that class is priced only with its calibration range spanning zero.
7. **The old instruments' record is the baseline.** Beside each row: what the instrument used at the time predicted
   before the board read (from the row's PREREGISTRATION / RESULTS), so the lab's hit rate is compared with the record's.

## What would invalidate a row
- An engine whose own known answers (801: 522/551 rows |d| < 1e-12; 803: 435 / 618 stats) failed.
- An arm whose recorded config differs from the board arm's in anything but the intended switch (rule 10).

## Amendment 1 (2026-10-15 22:55Z, before any checker number exists)
- **Bulk rows are decided on the COMBINED price when the arm is expressible on DE-H1:** combined = bulk block (pds + mse on
  the BULK folds) + DE block (fid + jac + nmae on DE-H1, the same arm vs the same base), their bootstrap draws summed
  (independent surfaces). Reason: 479sp2's board loss was nmae through ballast (+0.0018 raw), which no bulk fold can see;
  deciding on the bulk block alone would grade the checker on a member it cannot measure. Where the arm is not expressible
  on DE-H1, the bulk block decides and the row is marked "DE side effects unmeasured".
- The old instruments' record (`old_record804.py`, written before this amendment) is 8/13 hits; all 5 misses are false
  SHIPs (479sp2, 702, K20, 458d, 471c). The lab's checker is compared with that.

## Amendment 2 (2026-10-15 23:12Z, before any checker number exists): prospective rows
Three bulk arms are queued on production's A/B/C board after this protocol was written: **705** (`lodo_strength` θ 0.25 rider,
10-17 s1, vs its base), **707** (strong-tier swap β 0.25, 10-17 s2 / 10-18 s1, vs ALLSW), **710** (H1 0.70 = ×2.5, 10-18 s1, vs
ALLSW). The checker prices each with the same frozen settings as the retrodiction rows and **commits its prediction (point,
95 % interval, decision) to `experiments/804-known-answers/work/prospective.json` before the row's entry is submitted**
(git commit time is the proof). These are graded exactly like the others but reported separately: a prospective hit is
worth more than a retrodicted one because the checker could not have been shaped by the answer. If the engine cannot express
an arm in time, the row is recorded as "not priced", never filled in afterwards. Arms: 705 = pct × exp(0.25 z_src) above
the floor (522's `def_strength_key522` form); 707 = top z_src tertile rows blended toward the donors' unit mean U at β 0.25
(`experiments/707-strong-tier-direction/swap707.py`'s construction); 710 = composition H1 ×2.5 vs ×2 at CD4 ×1.5 (A/B/C rules).
