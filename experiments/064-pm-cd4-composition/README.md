# 064 — PerturbMatch × post-whitening CD4 composition

## Question

Do the two already-selected effect-level changes compose without retuning?

The frozen 2×2 comparison is:

- **A:** incumbent production union-whitened baseline;
- **B:** A with the K562 source slot replaced by the fixed 50:50
  incumbent/PerturbMatch processing ensemble selected in experiment 058;
- **C:** A plus the fixed post-whitening CD4 mean3 `gene_norm` residual at
  fraction `0.05` selected in experiment 062;
- **D:** both frozen changes.

Only Orion HCT116 and Orion HEK293T are valid target roles because
PerturbMatch and the incumbent K562 bank are reanalyses of the same Replogle
assay.

## Frozen protocol

- Full gallery: every challenge target shared by incumbent K562,
  PerturbMatch, both Orion banks, and all three CD4 conditions.
- All 300 challenge panel coordinates excluded from predictions, whitening,
  and retrieval.
- Gamma `0.65`.
- Candidate-independent source-only union whitener from incumbent K562 and
  the non-held Orion donor.
- Exact production pre/post-whitening median-norm preservation.
- PerturbMatch fraction `0.50`, applied only on genes shared with the
  incumbent K562 slot; incumbent values are retained outside that shared
  support.
- CD4 residual: equal mean of Rest, Stim8hr, and Stim48hr, `gene_norm`,
  fraction `0.05`, added after whitening with the experiment-062 global
  preserve scalar.
- Raw centered, untransformed held-Orion truth.
- No count emission, `cell-eval2`, challenge prediction, or submission.

The primary estimand is the macro PDS change `D - A`. The report also includes
both conditional main effects, the factorial interaction
`D - B - C + A`, and a paired target-cluster bootstrap for `D - A`.

