# Draft §3(i) — What an L1 error on significance-selected genes rewards (for Draft v11; merge after the 433a board read)

**F17 — An L1 fold-change error scored only on the reference's significant genes has a closed-form optimum under sign
uncertainty. It sits near the minimum detectable effect, and it can be bought without better biology. Provenance: the
statistics are textbook (the L1 minimiser is a median; selection on significance truncates |l| from below); the
application to a perturbation-prediction metric, the calibration on real data and the delivery are NEW. Confirmation:
H1 (5 seeds); the board read of 433a on 2026-09-24 is pending.**

*The metric.* `de_wilcoxon_lfc_nmae` scores, per perturbation, Σ|l̂ − l| / Σ|l| over the genes the reference calls
significant (Wilcoxon, BH 0.05, control mean > 5 CPM), with l the log2 ratio of per-cell-normalised arithmetic means.

*The optimum.* Take a gene whose predicted sign is right with probability a, and write F for the distribution of |l|
given significance.
- For a predicted magnitude m, the expected error changes at rate 1 − 2a·P(|l| > m).
- So the posterior median, and the L1-optimal m, is **m\* = F⁻¹(1 − 1/(2a))**.
- Below |l|, each unit of magnitude buys (2a − 1).
- Significance bounds |l| from below by the minimum detectable effect. For 400 perturbed against 18,400 control cells
  under Poisson counting, that is MDE ≈ 0.24/√(control count per cell).

*Calibration on H1.*
- On h1-test's gates (58 targets, 64,931 gene cells), the 10th percentile of |l| tracks 0.24/√count across the
  expression range.
- The pipeline's delivered magnitudes were ≈ 0.17× the truth, far below the MDE.
- Adding sign · min(c/√count, cap) to every steered gene wins an exponent scan at exactly ½ (the Poisson law).
- It is worth nmae **−0.0096 ± 0.0002** at c = cap = 0.1 (5 seeds), with the other five members flat.

*Delivery.* The add must not move the pseudobulk (the two direction members) or the Wilcoxon ranks (the three
call-set members).
- It is written into the per-cell class split of each gene's fixed total.
- The emitter's library classes make that a channel of about ±0.5 log2.
- Two constraints bind:
  - **class capacity.** A fine-step polish drained the middle library class, and in the line with the least slack (A)
    the plain delivery reached zero dose. Moving whole counts between the extreme classes only makes an A-balanced
    set of moves count-balanced, and restores the dose (0 → 0.96).
  - **call leakage.** Doses above c = 0.1 create hundreds to thousands of calls the design cannot suppress.

*What the optimum implies about everyone else.* Flip the true gate signs at rate 1 − a and give each gene its
L1-optimal magnitude. H1's nmae then maps to sign accuracy: 0.95 at a = 0.56 (this pipeline, measured), 0.88 at 0.62,
and 0.84 at 0.65.
- The public board's leaders sit at 0.82–0.84, so they run gate sign accuracy near 0.65, against 0.56 here.
- 71 % of the gate's |l| mass lies in genes this pipeline does not call, where its sign accuracy is 0.537.
- A fitted conditional-median model over our own statistics and the panel's donor tables adds nothing beyond the
  1/√count law (434).
- The remaining gap is direction information, not magnitude.

*Why this belongs in a paper about what references measure.*
- A benchmark member that reads fold changes only where the reference is significant rewards a calibrated magnitude
  floor, even from a predictor whose signs are barely better than chance.
- Its scaled value therefore mixes two quantities that should be reported separately:
  - sign information, the biology;
  - magnitude calibration against the selection threshold, an artefact of the selection.
- We recommend that benchmark designers report the sign-accuracy-equivalent alongside the scaled score, using the
  mapping above, which needs only the reference.

*Board confirmation (2026-09-24, entry `FCbBcvGxF1HOFPxPISUC`):* 433a − 417e: Δnmae = **−0.016786** raw, 1.76× H1's −0.00953 and inside the predicted −0.009 … −0.022. pds and mse are **bit-identical** (16 digits); overall +0.0050669.
