# 446bja — context A's block prior from its own line (Nadig Jurkat). Preregistration, 2026-09-25 03:50Z
Written before 446bja is built and before any 09-26 entry exists.

## The arm
446bja = 446b with `VCC_DE_SIGN_PRIOR=clean4q80ja`: context A's CPM-block prior is Nadig 2025 Jurkat's per-gene P_down
(396's artifact `d5dee640…`, 2,270 essential perturbations). B and C keep CLEAN4. q80 protection, the coarse MDE and the
zero repair are unchanged.
- Config diff vs 446b: exactly `de_sign_prior`, `de_sign_prior_sha256`, `de_sign_prior_context_A`.
- The bulk is ballast-set, so pds and mse are pinned by construction.
- On A's block, Nadig's sign differs from CLEAN4's on ~25 % of genes (CLEAN4 already gives Jurkat a quarter of its
  weight). The H1 analogue (Pan vs CLEAN4) differs on 27 %.

## Why this slot (the October question)
450 found that a same-line bank from another lab and chemistry carries the recipient's programme (H1: gate-block sign
0.625 vs CLEAN4 0.581). 451 builds an October same-line mode on that finding. This arm is the only board test of its
programme part. A is Jurkat, and Nadig's screen is essential-only, which 450 §5 shows carries the programme as well
(0.627).

## H1 price (450, the shipped harness, Pan's programme replacing CLEAN4 on all 85 targets, 3 seeds, paired)
| member | Δ raw | 95 % bootstrap over targets |
|---|---:|---|
| fid | **+0.0096** | [+0.0046, +0.0150] |
| nmae | **−0.0086** | [−0.0147, −0.0026] |
| reach | −0.0407 | [−0.0865, −0.0034] |
| jac, pds, mse | 0 | |
| overall (435's board slopes) | +0.0001 ± 0.0003 (seeds) | [−0.0098, +0.0091] |

- The reach loss is four small-gate targets. Most of it is CNN3 and COL6A1, whose gate is one gene each, flipping
  1.0 → 0. AASS (16 genes) and ZNF664 (7) make up the rest. 69 of 82 targets do not move. Without the two
  single-gene coin flips the overall reads +0.0048.
- fid improves on 49 / 85 targets and nmae on 38 / 58.
- Since the same-line prior is more accurate on gate genes, small-gate reach should rise in expectation. H1's four
  targets are too few to price it.

## Prediction (446bja − 446b, raw; A only, so about a third of H1's full-panel effect, then the observed transfer rates)
- 392's CLEAN4 step transferred at fid 0.67× H1, nmae 1.7× H1 (sign right), and reach with the WRONG sign.
- **fid +0.0015 … +0.0040**; **nmae −0.0015 … −0.0060**; reach −0.010 … +0.010 (unpriced, noisy); jac |Δ| ≤ 0.0005.
- **Overall +0.0005 … +0.0030.**

## Falsifiers (void the reading)
1. pds and mse equal 446b's to ≤ 1e-6 raw. The ballast pins them; a bigger move means the prior leaked into the bulk.
2. jac within ±0.0010 raw of 446b's: the call set is the |z| design, and the prior moves signs only.

## Bins (vs 446b, whatever 446b scores)
- **fid (primary):** Δ ≥ +0.0015 raw means the same-line programme transfers through the DE channel. The October
  same-line mode keeps its programme part.
  - +0.0005 … +0.0015: partial.
  - |Δ| < 0.0005: no transfer. Challenge A's block signs are not Nadig-like beyond CLEAN4's quarter.
  - Δ ≤ −0.0010: harm. The same-line mode drops its programme part and keeps only per-target rows.
- **nmae (secondary):** Δ ≤ −0.0015 confirms, |Δ| < 0.0015 is neutral, and Δ ≥ +0.0015 contradicts fid.
- **reach:** reported, not binned. H1 cannot price it (see above), and the board's ordering-channel sd is large.
- **Overall:** ≥ +0.0010 means 446bja becomes the incumbent and every later rider carries it. −0.0010 … +0.0010 means
  flat: keep 446b's prior for A/B/C, and decide the October programme on fid alone. ≤ −0.0010 means harm.

## Stated in advance
Displayed row not corrected. Upload ~12–19 min. Slot: 09-26 second slot, after 446b's capture, whatever 446b scores (the
LODO-key rider 446bj2 moves to 09-27). Model name "the vortex mixer wrote a sonnet", description empty.

## 446bja verified 05:11Z (still before any 09-26 entry)
Package `submissions/446bja-…-mde010coarse-zr5.vcc`, sha256 `512a7e41…`, novelty rc 0.
- Config diff vs 446b is exactly `de_sign_prior`, `de_sign_prior_context_A`, `de_sign_prior_sha256`.
- 0 zero-total gated genes.
- Render energy vs 410b 0.9603 / 1.0282 / 1.0785, identical to 446b: the bulk is pinned, as the falsifier requires.
- Extremes: 17 gated cells below −5 and 37 beyond −4.5 / +3, the same counts as 446b. The Jurkat prior creates no new
  extreme cells.
- The slot-2 waiter (`arm0926_s2b.sh`, pid 86417) finds its gates. An APFS clone of the artifact
  (`work/446bja-ka-copy.h5ad`) stays until 451's KA-2 passes.
