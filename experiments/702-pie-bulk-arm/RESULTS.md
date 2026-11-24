# 702: PIE bulk arm (ALLSW + PIE-blended bulk dump, seen targets, w 0.25). RESULT: HARM, not adopted

Written 2026-10-15T04:43Z by the science session [bf003b]. Class attempted: **transfer** (a better bulk direction from Arc's PIE). It did
not transfer. The PIE bulk family closes with 1 board try (this entry); held-out tries: 3 w × 2 surfaces (amendment 1),
plus 702c's post hoc read.

## Config (record, not description)
- = ALLSW (4WNmsNfdQhYx6bfib0dR) + `VCC_BULK_LOAD = work/dump_pie702_w025_seen.npz`
  (sha256 2bc63f58…, `make_pie_dump700.py`, 5 CPM rule, seen rule 154/300, norm-preserving blend w 0.25; amendments 1–3).
- Re-levelled r1 → ntclevel → r2: rho A 0.9788, B 0.9789, C 0.9815. Config diff vs ALLSW: exactly the 4 expected keys
  (amendment 3 added `bulk_sharpen_pt`'s float-noise edges; 0/300 targets re-tiered). CLIFF OK, zero-gene 0, novelty rc 0.
- Package `submissions/702pie702w25-…-pie.vcc`, content sha256 8b5202f841cb6784… (SUBMITTED_LEDGER), slot 2 of 10-15.
- Board entry **U7BlIZ53wTlLbaZTTET2**, captured in `data/board-entries/`.

## Board read (`read702.py`, the preregistered reader)
| | 702 | ALLSW | Δ |
|---|---|---|---|
| overall | 0.2735387 | 0.2760851 | **−0.0025464** |
| pds raw | 0.8632739 | 0.8674656 | −0.0041918 (scaled −0.00906) |
| mse raw | 0.6166028 | 0.6117205 | +0.0048823 worse (scaled −0.00512) |
| fid / reach / jac / nmae raw | | | +0.00001 / −0.00089 / +0.0000004 / +0.00006 |

- **Primary Δ(pds + mse) = −0.0023631 → HARM** (bin ≤ −0.0005). DE members −0.00018 overall: the re-roll class
  (reach rms 0.0013 raw, AGENTS rule 13), so not a falsifier.
- Prediction (amendment 1): primary +0.0008 … +0.0017, Δmse +0.0013 … +0.0031, Δpds −0.0005 … −0.0014.
  Realised: Δmse −0.00085, Δpds −0.00151 (overall units). **The mse sign was wrong**; the pds cost sat just past the
  bottom of the predicted range (it matches h1-test's mini-pds cost, not K562's).
- Follow-up rule (PREREGISTRATION-FOLLOWUPS.md): NULL/HARM → 702c only if pds and mse moved in opposite directions. Both
  fell, so **no 702c and no 702x**. `dump_hook` is decided **off** for October (coordinator, 99df0055).

## Why the held-out surfaces over-read PIE (`overread702.py` → `work/overread702.json`)
Target bootstrap (4,000 draws) of the ew-cosine gain of the same blend over each surface's base, split by the surface's own
control CPM:

| surface | targets | genes | gain (ew) | 95 % CI | base ew |
|---|---|---|---|---|---|
| K562 (unseen context) | all 269 | all 7,448 | +0.00139 (+3.7 %) | +0.00017 … +0.00256 | 0.0377 |
| K562 | all 269 | CPM ≥ 20 | +0.00224 | +0.00076 … +0.00363 | 0.0417 |
| K562 | all 269 | CPM < 20 | −0.00014 | −0.00158 … +0.00128 | 0.0204 |
| K562 | seen 135 | all | +0.00277 | +0.00044 … +0.00522 | 0.0443 |
| h1-test (SEEN context) | seen 44 | all 10,485 | +0.01825 | +0.01273 … +0.02344 | 0.0971 |
| h1-test | seen 44 | CPM ≥ 20 / < 20 | +0.01998 / +0.01047 | both CIs > 0 | |

Where 702's change sits on A/B/C (dump vs the cap-6 base, non-panel genes 18,233):
- share of ||new − base||² **outside K562's gene axis 0.407**, outside H1's kept axis 0.167, outside both 0.124;
- by A/B/C control CPM: [0,5) 0.081 of the change (0.170 of base energy); [5,20) **0.277** (0.141); [20,100) 0.477 (0.407);
  [100,∞) 0.166 (0.282). PIE moves the low- and mid-expressed genes twice as hard as the base does.

Mechanism, in order of weight:
1. **The unseen-context read was weak and was priced at its point estimate.** K562's all-target gain is +3.7 %
   [+0.45 %, +6.8 %]. Through amendment 1's own chain (518's +1 % X ≈ +0.0007 overall, k 0.5; K562's mini-pds cost
   −0.0005), its lower bound gives primary ≈ −0.0003 (NULL), not +0.0008. The range should have been built from the CI.
2. **h1-test was the optimistic end because H1 is a context PIE trained on.** Its gain is 6.6× K562's on the same
   construction, and every one of its CPM bins is positive while K562's low-CPM bin is null. Seen-context reads are not
   transfer evidence (AGENTS 11.1c); h1-test should have bounded the pds cost (it did: −0.0014 vs board −0.0015), never the
   gain.
3. **The instruments priced only part of the change.** About 41 % of the change energy on A/B/C lies on genes K562's truth
   does not cover, and inside K562's axis the < 20 CPM part reads null. The board scores all 18,233 genes, so the price
   rested on the 59 % the instrument could see.
4. **The translation was never validated for this family.** 518's identity was fitted on changes of our own donor
   direction. PIE had no earlier board point to retrodict (AGENTS 6.5), so the step from ew gain to board units had no
   known-answer row.

Point 3 explains why the mse sign could flip although every in-axis subset of K562 was ≥ null. Points 1–2 explain why we
expected a gain at all.

## What this changes for D/E/F pricing
- It does **not** undermine the decided October transfer switches. Each one has multi-fold held-out reads (3/3 folds,
  9/9 October-like re-draws, bootstrap CIs > 0) and the sealed R6 (`docs/october-decisions.json`). 702 had one unseen line
  with a CI barely above zero, and one seen line.
- New pricing rules for any **direction change** (proposed for FINDINGS via the coordinator):
  1. Price from the target-bootstrap CI, not the point estimate. A gain whose lower bound translates below the WEAK bin
     is a NULL prediction.
  2. A seen-context surface may bound a cost; it never supplies the gain end of a range.
  3. Report the share of the change outside the instrument's gene axis and by control-CPM bin. Above ~25 % outside the
     axis, call the price partial.
  4. If the family has no board point, the translation is unvalidated: widen the range to include zero.

## State of the PIE line
- PIE as bulk direction: **closed** (this read). PIE as DE sign prior: withdrawn on held-out (703, a0 0.5115 vs CLEAN4
  0.5171). PIE's DE-gene ranking (AUPRC 0.175 vs 0.148 on h1-test) is the one property never tested in a builder slot. It
  needs a call-set use, not a bulk use; see 700 RESULTS.
- Artifacts stay on disk (dumps, ctx files, cloud predictions) for the write-up; nothing of PIE enters the October recipe.
