# 399 — Preregistration: per-target neighbourhood block tilt, six-bank family, shrink 400 (`VCC_DE_SIGN_PRIOR=knn25o400q80`) on 395. Frozen 2026-09-16 ~05:40Z, before the build and before any entry exists (supersedes the 05:50Z draft that used the four-bank family with shrink 50; changed on 398c/398d's evidence, all pre-entry).

**The arm.** 395's exact configuration (`VCC_NORMREL_EXP=-1.5`, `VCC_NORMREL_MODE=cap`, `VCC_CD4_PARTIAL_G=0.5`, q80-protected
block flips) with the sign prior replaced. 395's prior is the marginal CLEAN4 tilt (one P_down per gene over every perturbation
of four banks). 399's prior is **per target**: for each panel target present in K562-GWPS (269/300), its **25 nearest GWPS
perturbations** (Pearson correlation of GWPS response rows, the target itself excluded) vote on every gene through their rows
in the **six-bank family** (K562-GWPS, K562-essential, Nadig Jurkat, Nadig HepG2, raw X-Atlas/Orion HCT116 and HEK293T; non-zero
entries only), shrunk toward the family marginal with prior weight 400; sign = down if P > 0.5. Targets absent from GWPS use
the family marginal sign; genes without any vote the CLEAN4 sign. Everything downstream (call set, q80 protection, C5b mass
balance, dose bisection, emitter) is byte-identical to 395 — the prior vector handed to the emitter is simply indexed by
target. Artifact `experiments/399-knn-tilt-arm/work/knn25o400_prior_challenge.npz`, sha256 `1711f71c…` (JSON sidecar has the
six bank hashes; rebuilt by `build_prior399.py 25 400 orion`; median 2,685 of 18,533 signs differ from CLEAN4 per target;
family marginal covers 98.4 % of the axis). Expected config diff vs 395's meta: exactly `de_sign_prior` (clean4q80 →
knn25o400q80), `de_sign_prior_sha256`, `de_sign_prior_targets_with_row` (absent → 269). Package
`submissions/399-cap15-partial50-knn25o400q80.vcc`, model name "the fax machine joined a book club".
**Reference for Δ: 395's entry, whatever its score** (second slot of 2026-09-16, chained after 395 is terminal). If 395's
falsifiers fail, the reference is 392 and the reading is voided as a one-variable test.

## Evidence before the entry (block-sign accuracy on the held bank's top-2,200 block, its own targets excluded)
| instrument | marginal tilt | this construction | Δ [95 % CI] |
|---|---|---|---|
| 398d LOLO, Nadig Jurkat (2,252 targets), six-bank marginal 0.5417 | 0.5417 | **0.5753** | **+0.0335 [+0.0317, +0.0355]** |
| 398d LOLO, Nadig HepG2 (2,143) | 0.5641 | 0.5858 | +0.0217 [+0.0202, +0.0233] |
| 398d LOLO, RPE1 (2,187; epithelial) | 0.5970 | 0.5988 | +0.0018 [+0.0003, +0.0033] |
| vs the CLEAN4 marginal that 395 ships (398/398c): Jurkat 0.5441 → 0.5753, HepG2 0.5625 → 0.5858, RPE1 0.5666 → 0.5988 | | | +0.031 / +0.023 / +0.032 |
| 397 H1 analytic (four-bank family, shrink 50): called-block precision 0.5225 → 0.530; fid +0.0174 vs delivered (CLEAN4 +0.0116); nmae −0.0008 (CLEAN4 +0.0017) | | | |
Random neighbour sets of the same size score below the marginal on every bank (the gain is the neighbourhood). The Orion rows
add an epithelial marginal (+0.030 on RPE1) that CLEAN4 lacks and cost 0.002 on Jurkat; contexts B (HeLa) and C (CAL33) are
epithelial carcinomas, A is Jurkat.

## Predictions (board units; slopes per raw unit fid +0.561, nmae −0.255, reach +0.188)
Block precision gain over 395's prior ≈ +0.03 on A and +0.02…+0.03 on B/C → panel fid ≈ +0.02 raw (+0.011 overall) if the
public-bank gains transfer 1:1; 392 realised 64 % of its H1 forecast, so **central Δ ≈ +0.007 overall**, plausible range
+0.002 … +0.014; nmae ≈ −0.003 raw (+0.001); reach unpriced (H1 says ≈ 0, board noise ±0.0045).

## Falsifiers (checked first)
1. `score_pds` and `score_mse` equal 395's to ≤ 1e-9. 2. jac within ±0.0005 raw of 395's.

## Primary reading: Δ score_fid vs 395
| bin | reading |
|---|---|
| ≥ +0.030 (≥ +0.009 raw) | the neighbourhood tilt transfers at ≥ half its public-bank size; 399 is the incumbent; next: the K / shrink dose and neighbourhood signs on uncalled genes (397c) |
| +0.010 … +0.030 | partial transfer; incumbent if overall ≥ +0.0005 |
| −0.005 … +0.010 | no transfer at this construction; A/B/C do not behave like the public lines on the block |
| ≤ −0.005 | anti-informative; close the per-target construction |

## Overall bins (Δ score_avg vs 395)
| bin | reading |
|---|---|
| **≥ +0.004** | transfers; incumbent |
| +0.0005 … +0.004 | partial; incumbent if fid is in one of its two top bins |
| −0.0005 … +0.0005 | flat |
| ≤ −0.0005 | harm; 395 (or 392) stays |

**Stated in advance.** The displayed row becomes this entry's score and is not corrected. Upload ~12 min.

## Addendum ~05:40Z (before the build) — H1 known-answer row for the shipped construction (`397/knn397d.py`)
Six-bank family, K=25, shrink 400, family-marginal default, signs on the called block, |E0| magnitudes: called-block precision
**0.5338** (CLEAN4 0.5223, six-bank marginal 0.5262); vs delivered E0: fid **+0.0207 [+0.0141, +0.0277]** (CLEAN4's construction
+0.0116), nmae −0.0011 [−0.0061, +0.0048], reach −0.0025 [−0.034, +0.029]; board-slope estimate +0.0114 vs E0, i.e. ≈ +0.0055
over 392's construction on H1. Bins above unchanged.

## Addendum 05:55Z (before the build's artifact exists) — builder patch applied and tested
Builder sha `76431184…` (commit ed4ee6bf): unset/clean4/clean4q80/down paths unchanged; the per-target lookup is inside the target
loop (verified by scope check). Builder unit tests with the stale CD4-residual test skipped: **56 passed, 0 failed, 1 skipped**
(`work/builder-tests2.log`; the plain runner aborts at that stale test since 2026-09-04). Build launched 05:49:53Z.

## Addendum 05:57Z (build running, no entry) — the shipped construction END TO END on H1 (400 harness, real emitter + scorer)
`experiments/400-uncalled-knn/work/metrics-k3knn.json` (cap15 known answer reproduced exactly first): pds 0.889636 / mse 0.974946
identical to cap15 (bulk untouched), **fid 0.516022 = +0.0185 raw vs cap15** (392's construction: +0.0122; 395's: +0.0115),
reach 0.271758 (+0.0174), nmae 0.950231 (−0.0006, better), jac −0.00004. The strongest DE-side H1 row of any construction so far.
Bins unchanged.

## Addendum 06:02Z (no entry exists) — first build crashed; emitter guard; rebuild
The first build (05:49:53Z) raised `UnboundLocalError` in `rc_emit_k.rc_hybrid_block_k` on context A's first target (ABCD1, one of
the 31 targets without a GWPS row, so it used the six-bank marginal default vector; `[de_sign] A: down_frac=0.8183` over the full
axis): the prior's demand was infeasible at every dose (λ = 0), and the pinned emitter's tolerance block dereferenced steering
diagnostics that do not exist on that path. Fix (commit ee9782bf): a target with no feasible dose is delivered unsteered with
`steer_unsteered=1`, λ = 0 recorded; every other path byte-identical (392's smallest λ was 0.087, 395's 0.265, so neither ever
reached it). Emitter re-pinned `dbbb9121…`, builder `3e1f0cc0…`; the config diff vs 395 therefore adds `de_sign_emitter_sha256`.
**New falsifier-adjacent diagnostic:** the meta's `prior_lambda` distribution — if a large share of targets sit at λ = 0 or
λ ≤ 0.5 (no flips delivered), the arm is under-dosed on A/B/C relative to H1 (400 harness: λ median 1.0, 2 % below 1) and a
middle-bin fid reading is dose, not transfer; the remedy is a less down-heavy default/shrinkage target, not a different prior.
Bins unchanged. Rebuild launched 06:02Z.

## Addendum 06:36Z — rebuild diagnostics, written after the build and BEFORE the entry exists (upload in progress); bins unchanged
Rebuild (guarded emitter) completed 06:33Z: config diff vs 395 exactly the four expected keys (PASS), novelty OK, package sha256
`5df81d42…`, artifact 11.79 GB. Real A/B/C dose: `prior_lambda` median 1.0, **min 0.0** (at least one target delivered unsteered —
the guard fired), `prior_shrink_k` median 0.70 (max 1.0: a target where the down-side shrink could not balance the block mass;
`prior_mass_residual` max 4,032 counts), flips per target median 902 (765–1,002), `sign_flip_frac` as 395. The summary carries no
per-target counts, so the λ ≤ 0.5 share on A/B/C is unknown; H1 (400 harness) had 2 % below 1. If fid lands in the middle bin, the
first question is dose (the down-heavy six-bank default competing for the high-library class), not transfer.
