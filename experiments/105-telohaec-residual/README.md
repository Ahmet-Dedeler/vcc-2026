# 105 — TeloHAEC residual adapter

Status: **extraction repaired; direct-source pre-gate NO-GO**

Date: **2026-08-30**

## Scope

This experiment evaluates the released IGVF TeloHAEC CRISPRi Perturb-seq
screen as a source-only residual adapter for the existing VCC predictor.

Hard restrictions:

- writes are limited to `experiments/105-telohaec-residual/` and
  `data/external/igvf-telohaec/`;
- no AWS or other cloud compute;
- no KOLF data or outcomes;
- no `cell-eval2`, count-level prediction build, challenge submission, or git
  commit;
- no VCC-panel outcome may select preprocessing, representation, model,
  regularization, scaling, gating, or blend strength.

## Frozen source

- Analysis set: `IGVFDS7943IIWZ`
- Cell-by-gene-and-guide H5AD: `IGVFFI8645YMGX`
- Guide library: `IGVFFI9754AGFB`
- Assay: CRISPR interference Perturb-seq
- Context: telomerase-immortalized human aortic endothelial cells
- Declared library scope: 2,285 targeted genes

Expected source metadata and checksums are recorded in `PROVENANCE.json`.

## Preregistered extraction

The extractor will operate in backed/chunked mode and will never materialize
the complete cell-by-gene matrix.

1. Verify that the expression matrix is sparse, finite, nonnegative, and
   count-like.
2. Resolve guide-to-target labels from the released guide manifest. Gene
   symbols are taken from the manifest's explicit symbol/description fields,
   with guide-name parsing used only as a checked fallback.
3. Retain single-guide cells when the H5AD exposes guide multiplicity. If it
   exposes only one assigned guide per row, use that assignment directly and
   record the limitation.
4. Form one pseudobulk per eligible guide by summing raw counts in row blocks.
   A guide requires at least 5 cells, per the count-only amendment below.
5. Construct the control profile as the coordinate-wise median of
   guide-level log-pseudobulk profiles across released safe-targeting and
   negative-control guides. Each eligible control guide receives equal weight.
6. For each targeted gene, subtract the control profile from every eligible
   guide profile, then take the coordinate-wise median across guides. Retain a
   target only when guides that individually pass the five-cell floor
   contribute at least 20 cells in total.
7. Use the project-standard transform

   `b = log1p(50000 * summed_counts / total_counts)`

   and store the raw uncentered effect `delta = b_target - b_control`.
8. Collapse duplicate feature symbols by summing their counts before
   normalization. Preserve the source gene axis; alignment to the VCC axis is
   a later modeling step.
9. Set every coordinate whose gene symbol is one of the 300 VCC perturbation
   targets to exactly zero in the exported effect matrix. Record and test this
   invariant.
10. Emit only compact source artifacts: `delta.npy`, `b_ctrl.npy`,
    `genes.txt`, `perts.txt`, `guide_qc.csv`, and `extraction_report.json`.

Primary extraction choices above are frozen. Schema-specific column names may
be selected after inspection, but their semantic roles and every threshold
are fixed here.

## Count-only protocol amendment

Recorded on **2026-08-30**, after inspecting only guide labels and cell counts
and before reading any expression values.

The preregistered minimum of 10 cells per guide was reduced to **5 cells per
guide**. The released matrix is unusually sparse at the individual-guide
level: among targeting guides, the observed cell-count quantiles are
`1, 2, 3, 5, 8, 11, 35`. Under the repaired rule that the 20-cell target
minimum counts only guides meeting the guide floor, the exact retained-target
counts are **1,567 at 10 cells/guide** versus **2,167 at 5 cells/guide**. The
five-cell floor therefore retains 600 additional adequately represented
targets.

This amendment used no expression values, downstream effects, VCC outcomes, or
held evaluation truth. All other extraction choices and the 20-cell target
minimum remain frozen.

## Extraction eligibility repair

A read-only audit found that the first extraction counted cells from
sub-threshold guides toward the 20-cell target minimum even though those guides
were excluded from the exported pseudobulk. That run retained 46 targets with
fewer than 20 contributing cells, including ten supported by only one eligible
guide.

The repaired extractor applies the 20-cell minimum to the sum over guides that
individually pass the five-cell floor. The compact artifacts and hashes below
were regenerated after this correction.

## Completed extraction

The repaired bounded extractor and five synthetic tests completed successfully on
**2026-08-30**.

- Source shape: 214,449 cells by 17,472 feature rows.
- Sparse entries checked: 638,734,499.
- Matrix validation: finite, nonnegative, integer-valued; zero explicit zeros;
  maximum stored count 1,977.
- Feature axis: 17,464 unique symbols after collapsing eight duplicate rows.
- Working axis: 17,164 genes after removing all 300 VCC panel coordinates
  before library-size normalization or any effect calculation.
- Eligible guides: 19,564 total, comprising 19,015 targeting and 549 control
  guides.
- Retained clean targets: 2,167, each with at least 20 cells from guides that
  individually pass the five-cell floor.
- VCC target coverage: 65 of 300 as source perturbations; all 300 symbols are
  present on the source expression axis.
- Export: `delta.npy` has shape 2,167 by 17,464 and `b_ctrl.npy` has length
  17,464, both float32.
- Exact-zero verification passed for all 300 panel coordinates in both
  exported arrays.
- Full repaired extraction runtime: 18.27 seconds wall, 17.01 seconds user, 0.98
  seconds system.
- Peak resident set reported by `/usr/bin/time -lp`: 3,593,207,808 bytes
  (3.35 GiB); no swaps.
- The temporary disk-backed guide-sum matrix was removed after completion.

Artifact byte sizes and SHA-256 hashes are frozen in `PROVENANCE.json` and
`data/external/igvf-telohaec/extraction_report.json`.

## Preregistered residual-adapter evaluation

The more flexible residual adapter below remains design-only. It was not run
because a stricter direct-source gate failed on both held somatic contexts; see
`RESULTS.md`.

For each outer held donor `d`, the held donor is test-only: none of its
perturbation values may enter a fitted basis, scale, residual label, model,
hyperparameter choice, or selection decision. Training residuals come only
from the other contexts. For each training context `c`, construct its incumbent
while excluding both `c` and `d`, then learn a low-capacity ridge map from
TeloHAEC effects to those non-held-context residuals.

The deterministic non-panel split is assigned once by target symbol and reused
across every training context, so a target identity cannot cross fit, tune, and
inner-test partitions. Fit and tune use only non-held contexts. The outer held
donor remains untouched until the selected model and all preprocessing are
frozen.

For every validation gallery `G`, every coordinate whose gene symbol belongs
to a perturbation target in `G` is zeroed in TeloHAEC inputs, incumbent and
residual vectors, fitted bases and scales, predictions, truth, and PDS. For the
challenge evaluation, `G` is exactly the fixed 300-target panel. This prevents
the perturbed gene's own coordinate from revealing target identity and applies
the same masking policy used at deployment.

The frozen search space is:

- deterministic non-panel target split: 60% fit, 20% tune, 20% untouched test;
- input ranks: 16, 32, 64;
- residual-output ranks: 8, 16, 32;
- ridge strengths: `1e-3, 1e-2, 1e-1, 1, 10`;
- blend fractions: `0, 0.02, 0.05, 0.10, 0.20`;
- identity (`blend = 0`) is mandatory;
- every fitted basis and scale excludes the held donor and the complete target
  coordinate mask for its validation gallery.

Selection uses only non-held-context, non-panel fit/tune targets. The untouched
inner non-panel test and the completely untouched outer donor must pass before
any one-time panel evaluation is authorized.

Promotion gate:

- nonnegative PDS delta in every outer held-donor fold;
- macro raw-PDS delta at least `+0.003`;
- target-cluster bootstrap 95% lower bound above zero;
- observed gain above the 95th percentile of target-label permutation nulls
  stratified by source coverage and effect norm;
- unsupported targets remain exact incumbent fallbacks.

No held evaluation truth is opened by the implementation/acquisition stage.
