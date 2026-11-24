# 458 — Arc's own programme as every context's block prior (preregistration, 2026-09-26 04:50Z)
Written before 458 is built and before any 09-27/28/29/30 entry exists.

## The arm
458 = 446b + `VCC_SAMELINE=experiments/458-arc-lab-programme/work/arch1_programme_abc.npz:32b82850…:0`.
- The artifact (`artifact458.py`) gives A, B and C the same block-prior default: P_down over Arc's 2025 H1 screen
  (`h1_2025`, 235 non-panel perturbations, 10x Flex, 400 cells each; 396's rule). It has no rows and w 0.
- The builder path is 451's KA-2 path (programme default via VCC_SAMELINE), proved equal to `clean4q80ja` for A.
- q80 protection, the coarse MDE, the zero repair and the bulk are unchanged.
- Expected config diff vs 446b: exactly `sameline`.
- The prior's sign differs from CLEAN4's on 39.7 % of all genes.

## Why
- 450h/450k: among cancer lines a bank's programme is mostly its lab's and protocol's signature. Weissman lines agree at
  r 0.43–0.71, while K562 across labs agrees at 0.09. So the right programme source for Arc's contexts may be Arc's
  own lab, not another lab's screen of the same line (446bja's null).
- 450l: Arc's H1 programme shares r 0.13–0.20 with both Weissman and Xaira programmes. Those two labs agree with each
  other at ≈ 0.
- If this transfers, it is a lever for **every** October context: D/E/F come from the same lab and protocol. No same-
  line bank would be needed.

## Prediction (458 − 446b, raw)
Uncertain by design; it tests a hypothesis.
- fid −0.001 … +0.004. The first prior (CLEAN4 through 392) bought fid +0.0082 raw on the board.
- nmae −0.004 … +0.002. reach ±0.010 (the small-gate coin flips; unpriced).
- pds and mse identical (bulk untouched).
- Overall −0.002 … +0.003.

## Bins (vs 446b)
- **Falsifiers:** pds and mse must equal 446b's to ≤ 1e-6 raw, and jac must stay within ±0.0010.
- **fid (primary):**
  - ≥ +0.0015: the same-lab programme carries across lines. Carry it to October as every context's default where no
    same-line bank matches, after a D/E/F-safe artifact (it must not use D/E/F rows).
  - +0.0005 … +0.0015: partial. Try a CLEAN4 + Arc-H1 blend on the board before October.
  - |Δ| < 0.0005: no. Arc's H1 programme is ESC-dominated or the lab component is too small.
  - ≤ −0.0010: harm (ESC-specific signs). Close it.
- **nmae (secondary):** ≤ −0.0015 confirms; ≥ +0.0015 contradicts fid.
- **Overall:** ≥ +0.0010 means adopt for A/B/C riders; ≤ −0.0010 means close.

## Slot
09-30 or later (the rule is fixed after the 09-27 reads, before the entry exists). Name: a fresh joke. Description: empty.

## Addendum 05:20Z (before the build and before any entry): the dose
On each context's CPM block (2,200 top-CPM genes of its own controls), Arc's H1 programme differs from CLEAN4's sign on
**39.4 / 41.1 / 41.7 %** of genes (A / B / C). It is much less down-skewed: 51–52 % down against CLEAN4's 66–68 %.
- q80 protection keeps the strongest 20 % of each block at its delivered sign, so ~32 % of block genes change their prior.
- That is a larger intervention than 446bja (A only, ~25 %), so the read should be decisive either way.
- The harm bin (fid ≤ −0.0010) is the likely outcome if the ESC-specific content dominates.

## Addendum 09:35Z (before any 09-27…09-30 entry): the motivating numbers, corrected
450h v3 re-measured the lab-vs-line split on matched perturbation sets. Different lines within a lab agree at r 0.54; the
same line across labs at 0.21 (not 0.09 as written above). The lab remains the larger factor (~+0.37 vs the line's ~+0.2),
so the hypothesis and the bins are unchanged.
- 450l on shared sets: Arc-H1's programme agrees with each other lab at r 0.10–0.25, the generic cross-lab level. Other
  labs' within-lab cross-line agreement is 0.54 (Weissman, Xaira). If Arc's is similar, the same-lab prior should beat
  CLEAN4 (cross-lab ~0.18) on A/B/C. The ESC-specific content is the stated risk. Bins unchanged.

## Addendum 10:23Z 09-26 (before any 09-27…09-30 entry): 458b is dropped, sized analytically
`resid458b.py` builds the handoff's follow-up: Arc-H1's z-scored P_down with Pan-H1's regressed out (β 0.261 on
15,981 genes), keeping Arc's level.
- On each context's 2,200-gene CPM block it flips only **8.0 / 8.4 / 8.6 %** of 458's signs (A / B / C).
- 458 itself flips about 40 % relative to CLEAN4.
- After q80 protection, 458b − 458 differs on about 6 % of block genes, so it cannot read distinguishably from 458 on
  the board. It is not worth a slot in either branch.
- A 458 null means the lab component is too small, not that the ESC part diluted it.
- Context (462, same day): on the Broad-16, one lab's cancer lines share about a third of their honest programme
  (cross-line r 0.043 vs within-line 0.121, after removing shared-control noise). That supports a modest 458 effect,
  inside the preregistered range. The bins are unchanged.

## Addendum 10:30Z 09-26 (before 458d is built, before any 09-27…09-30 entry): 458d replaces 458 in the 09-30 slot
**Why.** 462c: the `h1_2025` bank's programme is its harness's control noise to a large degree.
- Setup: 400 cells per target against 18,399 control rows (15,571 unique cells), thinned to ~16k UMI.
- Expressed genes: split-half r 0.79 with a shared control against 0.39 with disjoint controls. On the A/B/C blocks,
  0.87 against 0.48.
- About a third of 458's programme variance is one control sample's noise, which carries no information about A/B/C.
- A null from 458 would be ambiguous, and a positive would still need the clean programme tested. So the clean one is
  the right single test of the hypothesis.

**The arm.** 458d = 446b + `VCC_SAMELINE=experiments/458-arc-lab-programme/work/arch1d_programme_abc.npz:6951dbc7…:0`.
It is 458's construction with the programme from ALL of Arc's public H1 cells (`denoise458d.py`):
- 300 perturbations × ~1,071 cells (275 off the A/B/C panel) against all 38,176 control cells, at raw depth. KA: all 300
  cell counts equal the official `pert_counts_*.csv`.
- The bank is `data/donors/h1_2025_full` (delta `26662c52…`).
- Disjoint-control reliability 0.69 (458's bank: 0.39).
- Against 458's programme: r 0.68 on expressed genes (0.77 above 100 CPM). It flips 21.8 / 20.8 / 21.6 % of block signs.
- It agrees with CLEAN4 on 0.650 / 0.611 / 0.591 of block genes (458: 0.604 / 0.589 / 0.586).
- Config diff vs 446b: exactly `sameline`.

**Prediction and bins: 458's, unchanged** (fid primary: ≥ +0.0015 carries; +0.0005 … +0.0015 partial; |Δ| < 0.0005 no;
≤ −0.0010 harm; falsifiers: pds and mse equal to ≤ 1e-6 raw, jac within ±0.0010). Expected: at least 458's effect on fid
if the lab hypothesis holds. It is a slightly smaller intervention vs CLEAN4 (35–41 % of block signs differ).

**Slot.** 09-30 s1 as "the pipette learned to meditate". Fallback: 458 if 458d is not built and gated. s2 stays 459dr.
October: if it carries, `DEF_DAY_PROGRAMME_SOURCE=h1_2025_full` (the denoised bank), not `h1_2025`.

## Addendum 10:49Z 09-26 (before 458e is built, before any 09-27…09-30 entry): 458e, the dose-curve partner
**The arm.** 458e = 446b + `VCC_SAMELINE=experiments/458-arc-lab-programme/work/arch1e_blend_abc.npz:b67af1f5…:0`.
- The block prior is sign(0.5·z(Arc-full P_down − 0.5) + 0.5·z(CLEAN4 P_down − 0.5)), z-scaled per context on its own
  CPM block (`artifact458e.py`).
- It flips 15.2 / 12.8 / 13.8 % of 458d's block signs and 19.8 / 26.1 / 27.1 % of CLEAN4's. It sits midway on the
  w = 0 (446b), 0.5, 1 (458d) curve.
- Config diff vs 446b: exactly `sameline`. On H1 (Arc's own line), w 0.5 reads fid-like 0.5456 (w 1: 0.5545, CLEAN4:
  0.5244; 462g).

**When it goes.** Only if 458d reads positive or partial (fid ≥ +0.0005 vs 446b). Then the two points set the weight
October ships. It is never submitted after a null or harmful 458d.

**Prediction and bins (458e − 446b, raw).** pds and mse identical to ≤ 1e-6. fid between 0.4× and 1.1× of 458d's Δfid.
- Δfid(458e) > Δfid(458d): the consensus adds information. October uses w 0.5 (a blend source).
- 0.4–1.0× of 458d's: the curve is monotone. October uses w 1.
- < 0.4×: non-linear. Keep w 1, and note it.
- Name: "the lab coat filed a complaint". Description: empty.

## Addendum 11:31Z 09-26 (before any 09-27…09-30 entry): the reach contingency
A programme prior can buy fid and pay reach. 446bja's board read had reach −0.0079 raw, and Pan's same-line programme on
H1 had fid +0.0097 with reach −0.041. The prior flips signs on block genes that rank high in our significance order,
and reach needs 90 % sign purity along that order.
- If 458d reads fid ≥ +0.0015 **and** reach ≤ −0.005 raw with overall ≤ 0: fid carries, but the October switch is NOT
  set as is. The follow-up is a reach-protecting variant (a larger q protection on the programme's block, 395's
  mechanism), board-tested first.
- 462l (queued after 458e) prices Arc's own programme on H1's harness, all six members, as the upper-reference
  profile.

## Addendum 17:04Z 09-26 (before any 09-27…09-30 entry): 458g joins 458d on 09-30
458g = 458d + the strength gate `VCC_DE_PRIOR_GATE=q33` (467 on H1). It takes 09-30 s2, read against 458d the same day.
Its bins and the combined October table (458d's P × the gate's read) are in `468-prior-gate/PREREGISTRATION.md`. 458d's
bins here are unchanged.

## Addendum 2026-09-28 19:40Z (coordinator; before any 09-29 or 09-30 entry exists): how to diagnose the jac falsifier
Source: worker 517 §0.1, verified here on the build metas. The bins and the falsifier thresholds are unchanged.
- Arc's programme as the block prior makes the K-class emitter leak calls it cannot suppress. The build metas record
  `output.rc_callset.call_set_unsuppressible`: **446b median 0, max 86; 458d median 1, max 1,655; 458g median 0, max 1,655**.
  On H1 the same leak cost jac −0.0027 (ungated) / −0.0022 (gated) raw and fid ≈ −0.0016.
- **Pre-specified diagnosis:** if the falsifier's jac leg fires (|Δjac| > 0.0010) with Δjac in **−0.001 … −0.004 raw**,
  and pds/mse are identical to ≤ 1e-6, attribute it to this leak (not a bug) and read fid, nmae and reach as
  preregistered, noting the leak's fid cost (≈ −0.001 … −0.002 raw) beside the fid bin. Any other falsifier pattern
  (pds/mse moving, jac positive beyond +0.0010, or jac below −0.004) remains a bug: stop and diagnose before reading.
