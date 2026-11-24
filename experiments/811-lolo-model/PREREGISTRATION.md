# 811 — preregistration of the BASELINE (step 2 of the mandate), written before any 811 number exists

Written 2026-10-16 ~07:45Z by the 811 worker (Claude Opus 5.5). Scope: only the baseline the mandate requires to be beaten
first ("the consensus of all non-held lines' unit responses, equal weight, blended at λ into the composition, priced on
BULK-ABC/W"), its control, and one non-shippable ceiling diagnostic. The learned model waits for PROPOSAL.md and the judge.

## Class
**Transfer** (AGENTS 11.1c.1): a different prediction of the held line's response direction from public lines. Nothing in it
is fitted to A/B/C or to any fold truth. The day-of path is identical: D/E/F's panel targets' rows from the same public banks.

## Construction (frozen)
Per fold (panel P, held line h) and per target t of the fold:
- Line groups and their banks (the inventory, `work/inventory811.json`): K562 = {k562_gwps, viperturb_k562, jiang_K562};
  HCT116 = {orion_hct116}; HEK293T = {orion_hek293t}; CD4 = {cd4_rest, cd4_stim8hr, cd4_stim48hr} (the panel's GWCD4i
  compact banks); H1 = {h1_pan2025 (115's GSE295214 compact donor), h1_vcc2025 (data/donors/h1_2025_full)}; KOLF2 =
  {kolf2 local rows}; TeloHAEC = {igvf-telohaec pilot}; DLD1 = {223's dld1_lowmoi}; A549 = {223's altair_a549, jiang_A549};
  RPE1 = {rpe1}; Jurkat = {nadig_jurkat}; BXPC3, HAP1, HT29, MCF7 = {jiang_<line>}.
  Excluded by construction (not by any number): perturbmatch (the same K562 cells as k562_gwps), LINCS L1000 tables (inferred
  genes, a different assay; 176 closed), IGVF cardiomyocyte (significant hits only, median 4 genes per row), k562_essential
  (0 panel targets).
- The held line's whole group is out (on K562 folds every K562 bank is out; on ABC/kolf2 the KOLF2 group is out).
- Bank row r_b(t): the bank's delta row mapped to the official 18,533-gene axis, panel columns zeroed (as 706/707); unit-
  normalised. Line row = normalise(mean of its banks' unit rows that carry t). U_all(t) = mean over lines carrying t of the
  line rows (equal weight per LINE). U_src(t) = the same over the fold's composition SOURCES only (707's U, all targets).
- Arm: on every target with U(t) ≠ 0, the composition row is blended norm-preserving on the non-panel genes,
  fin'_t = |fin_t| · normalise((1 − λ) fin_t/|fin_t| + λ U(t)/|U(t)|) (707's formula, all targets instead of the top tier),
  then the incumbent's key w, sharpening q, render and ballast unchanged (`engine801.price`, base = the incumbent).

## Arms and tries (family `811-consensus`, 6 priced configurations + 2 known answers + 1 diagnostic)
| arm | what | role |
|---|---|---|
| `inc` | the fold's incumbent (804's `inc`) | KA: arrays must equal 804's `inc.npz` exactly |
| `707` | 804's 707, rebuilt by my code (sources, top z_src tertile, β 0.25) | KA: arrays must equal 804's `707.npz` exactly |
| `all_l10`, `all_l25`, `all_l50` | U_all at λ 0.10 / 0.25 / 0.50 | **the baseline** (3 tries) |
| `src_l10`, `src_l25`, `src_l50` | U_src at λ 0.10 / 0.25 / 0.50 | control: separates "more lines" from "all targets" (3 tries) |
| `oracle_vip` (K562 folds only) | λ 1 swap to viperturb_k562's own K562 rows where it carries t | NOT shippable (same line as the truth): ceiling diagnostic for any cross-line model |

Abort rule: if either KA row differs from 804's arrays (any |d| > 1e-12) the run is void and no number is read.

## Predictions (bins, before the numbers)
Board-unit overall from `check800.py aggregate` (pds literal, mse via 518), pooled W (18 folds) and ABC (3 cancer folds; kolf2
reported apart).
1. **Baseline vs the incumbent.** 707 (sources, top tier, β 0.25) read W +0.0049, ABC +0.0004 (WEAK). 706: on A/B/C larger
   β and more targets are worse. Prediction: W positive at every λ (largest at λ 0.25–0.5: +0.003 … +0.008); ABC ≤ +0.0005 at
   λ 0.10 and negative at λ 0.50 (−0.001 … −0.004); ABC/kolf2 negative at every λ.
2. **More lines (all_lX − src_lX at the same λ).** Sources saturate (218/225), CD4 and H1 already sit in the composition,
   and the only full-coverage extra line on disk for L folds is a second K562 bank. Prediction: |Δ| < 0.0005 overall on W
   and on ABC (NULL), with a negative lean on cancer truths (706's mechanism: Arc-like / non-cancer directions do not help
   public cancer truths).
3. **Same-line ceiling.** viperturb (Flex, ~48 cells per target) as a predictor of K562-GWPS truth. 424 found it loses as a
   donor on PLOLO. Prediction: the λ 1 swap reads r between −20 % and +10 % against the incumbent, i.e. a second lab's
   same-line rows are NOT a better direction than the composition at this depth. If instead r > +20 %, line-specific
   biology is large and a cross-line model has real headroom.

Bins: SHIP-like = 2.5 % bound > 0 after the best-of-6 shift (e(6) ≈ 1.27 sd); NULL = interval spans 0; HARM = 97.5 % bound < 0.
**Falsifier for the lane's premise ("more public lines carry more transferable direction"):** prediction 2 reads NULL or HARM
on both W and ABC.

## What is unmeasured and why
fid, jac, reach, nmae: the bulk engine does not measure them. Mechanism: ballast K (re-roll class). The engine reports
dK_rel; any arm with |dK_rel| > 5 % is flagged for DE-H1 before any dossier.
