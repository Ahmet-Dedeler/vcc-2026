# 617 — Which FINDINGS §2 "Closed" families rest only on thinned-Orion DE readings? (read-only)

[21dd4c], 2026-10-02 ~20:00Z, for the coordinator (after 616's finding that the thinned Orion harnesses' DE truth was ~98 %
selection artefact). Nothing was recomputed.

## Method
1. **Cited experiments.** `scan617.py` (`work/scan617.json`) lists every experiment cited in each §2 bullet. The August
   bullet names methods without numbers, so those were mapped from the history sections: SVD 004 / 017 / 029 / 030 /
   073 / 078; EB 010 / 022 / 032; MODZ 016; inverse-variance 015 / 022 / 078; per-donor whitening 024; log1p
   re-expression 026; noise-covariance whitening 021; CD4 blends 047 / 050 / 051 / 054 / 062; de-whitening 063 / 098.
2. **Code-level test** (`work/codehits617.json`, precise). An experiment counts as an Orion-DE reader only if its CODE
   loads one of:
   - a thinned / greedy Orion harness (`orion-*-thinned.h5ad`, `orion-hct116.h5ad`, `-depthmatched`, `refde_orion`);
   - 480's `gate-orion_*` tables, directly or through `fold480` / `inst636` / `ds612b`.
   Mentions of Orion in a RESULTS are not enough: the Orion DONOR tables and PLOLO folds are clean (616: pooled controls,
   blob share 0.064).
3. **Reading.** Each flagged experiment's deciding evidence was read by hand.

## Answer: no §2 closure rests ONLY on thinned-Orion DE readings
Among the cited experiments, the code-level readers are 054, 067, 192, 477 / 480, 610, 612 and 613. (The 040 / 100 hits
were number-parsing false positives.)

| §2 family | Orion-reading evidence | what actually decided it | verdict |
|---|---|---|---|
| CLEAN4 bank composition (613f / 613g) | 613's block instrument, whose 458d / 392 validation was 2 / 3 Orion folds. 613c's Orion-only "Jurkat alone" hint | the READ is on 12 independent Broad lines (Jurkat-alone wins 1 / 12; every subset and superset loses) | **stands** (status quo). Caveat: its instrument's validation is now K562-only (616), so the Broad-12 read is "consistent", not "validated" |
| Disagreement shrink (612 ds612 / ds612b) | the fold model's Orion rows (already discredited: 612 CORRECTION 2, "do not cite") | H1 real emitter: KGID vs no hook −0.0015 … −0.0046 in every regime | **stands** |
| 480's consensus rider (477) | 480's fold analysis motivated the rider | analytic on H1 through 490's exact reconstruction | **stands** |
| Per-cell CPM coupling (610 cpmcouple610) | none (decile610, a different reading, used 636's surfaces) | H1 truth + 435's emitter stats | **stands** |
| "Everyone over-calls" (192: K < n_real is catastrophic via fid's yield) | §4.4 only: a 3-feature membership logistic fitted on Orion HCT116-thinned "does not transfer" (2.01 vs plain cpm 2.05) | h1-test call-set curves + fid arithmetic + board constants | **stands.** The side result is Orion-only; see re-openable item 1 |
| Control-variance transport (067, August) | HCT116's context controls came from `orion-hct116.h5ad` (greedy) | LODO macro-PDS deltas on donor-table folds (−0.0002 / −0.0004): bulk, ≈ 0 | **stands** (bulk; the DE artefact does not apply) |
| Pre-whitening CD4 blending (054, August) | scored on the depth-confounded HCT harness (bulk cosine) | superseded by 062's post-whitening CD4 residual (production) | **moot** |
| All other §2 bullets (630 capacity, 651, 645, q80 484 / 490 / 475, LODO keys 541 / 542 / 544, 545, 494, 473 / 414, 483, 211 / 216 / 317, 318, 325–342, 343–382, 09-11 / 09-13 / 09-16 / 09-22 / 09-22-23 sets, fifth donor line) | no code reads of thinned Orion | H1, board, PLOLO (clean donor tables), kolf2, Broad | **unaffected** |

Notes outside §2, for completeness:
- **§1 "434's closure stands on two instruments":** the second instrument (480's cancer folds) is now K562 only.
  434 itself is H1-decided.
- **History, 098:** "every nmae calibration route is closed" was read on h1-test + the Orion harness. It was already
  overturned by later board evidence (coarse MDE −0.0168, zs −0.0022).
- **History, 153 / 156 / 157 / 159:** they read the thinned harnesses. H.4 (09-02) already restricted them, and 157
  itself preferred h1-test on reach.
- **History closures of the sign and ordering channels and operator geometry** (193 / 212 / 226 / 228–234 / 243): no
  Orion harness reads.

## Re-openable (Orion-only evidence), all low priority
1. **A fitted call-membership model transferred across lines (192 §4.4).**
   - Its "does not transfer" was read only on Orion HCT116-thinned, whose DE calls were artefact.
   - Re-test on K562 (614's reference) and the Broad lines before calling it closed.
   - Expected value is small: on H1 the fitted model beat plain cpm only same-line (2.77 vs 2.05).
2. **604 v3's exclusion of identity-only genes** (from 616's list, not a §2 bullet).
   - The stated evidence (0.27–0.34 sign accuracy) was HEK-only.
   - The exclusion is conservative. The identity-only extension is re-openable on H1 / Broad-16.
3. **622's "public 3′ references alone are harmful"** (HEK-only; 7772dc withdrew it at 30866a06; not in §2).

## Why so little of §2 is exposed
- Since 09-02 (H.4), §2 closures have been decided on h1-test, the board, PLOLO's donor-table folds or the Broad lines,
  not on Orion harness absolutes.
- The council-era Orion DE readings (480 → 636 → 613 / 639 / 612 / 516 …) mostly fed **instruments and open calls**
  rather than closures. 616's audit covers those, and the coordinator has re-ruled them.
