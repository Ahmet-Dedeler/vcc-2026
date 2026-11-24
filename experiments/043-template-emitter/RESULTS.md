# 043 — Fixed-margin control-template emitter

**Date:** 2026-08-22. Design/prototype only. No `cell-eval2` run, real-data build, or external
submission was performed.

## Finding

The current `even` and `quasi` emitters fix each gene's panel total, so they protect the
`bulk_lognorm` profile used by PDS and MSE's plug-in distance. But both replace the selected
controls' library-size distribution with one common depth, `median(lib_sizes)`. They also
synthesize each gene largely independently of the real control cells' zero pattern. This is why
the resulting per-cell CPM distribution has only a few distinct values and why its arithmetic
mean CPM can disagree with the same profile's count-weighted pseudobulk.

## Implementation audit

The local environment has `cell-eval2 0.16.0`. The relevant implementation paths confirm:

- `src/vcc2026/donor_transfer.py::emit_counts`: `even` and `quasi` set
  `L = median(lib_sizes)` and derive fixed per-gene totals from `n * L * profile`; the individual
  sampled control depths are not used.
- `cell_eval2/prep.py::bulk_lognorm_means`: PDS/MSE bulk profiles depend only on raw group column
  sums and the grand total.
- `cell_eval2/moments.py::jackknife_correction`: MSE's prediction-side correction is dense in
  cells × genes and depends on the complete cell layout, not just margins.
- `cell_eval2/de_compute.py::_to_linear`, `_grouped_means`, and `_lfc_from_means`: the DE path
  normalizes every cell to one million counts and uses the arithmetic mean of those values for
  log2 fold change.
- `cell_eval2/de_compute.py::_apply_cpm_filter`: the >5 CPM tested-gene gate comes from the real
  control alone; emission cannot change the tested universe, only ranks, fold changes, and calls.

`template_emitter.py` tests a different construction:

1. Select 400 real control cells from the target context. Their row sums are the requested row
   margins `r_i`.
2. Convert the model profile to integer column margins `c_g` with Hamilton rounding at
   `S = sum_i r_i`. If bit-identical compatibility with an existing build matters, pass that
   build's exact `c_g` instead and Hamilton-rescale the control row margins to the same `S`.
3. For every gene observed in the selected controls, allocate `c_g` over cells in proportion to
   that gene's control counts. This preserves the observed zero support until integer rounding
   makes a zero unavoidable. A target-positive gene absent from all 400 templates falls back to
   depth-proportional allocation.
4. Repair row margins by moving counts from excess rows to deficient rows **within the same
   gene**. Donor and recipient rows are matched by nearest inverse depth, and transfers first use
   genes already supported in the recipient template.

The result is an integer matrix with both margins exact.

## Why the two margins matter

For emitted counts `X`, the scorer's expression profile is

```text
b_g = log1p(50,000 * sum_i X_ig / sum_ih X_ih).
```

Therefore exact column margins `sum_i X_ig = c_g` make PDS/MSE pseudobulk invariant to every
cell-level rearrangement. The target floating profile can only be represented approximately by
integer counts; Hamilton rounding is optimal in the standard largest-remainder sense and leaves
every column within one count of its quota. Once `c_g` is chosen, the pseudobulk is exact.

The DE path instead ranks

```text
CPM_ig = 1,000,000 * X_ig / r_i
```

and computes log2 fold change from `mean_i CPM_ig`. Exact row margins reproduce the selected
controls' library-size distribution. Allocating each gene according to its control template
approximately preserves the gene-specific correlation between detection/expression and cell
depth. Before row repair, the continuous target for an observed gene is

```text
meanCPM_target,g =
    (c_g / control_total_g) * mean_i(1,000,000 * control_ig / r_i).
```

A one-count repair move from row `i` to row `j` changes that gene's arithmetic mean CPM by

```text
(1,000,000 / n) * (1/r_j - 1/r_i).
```

Matching donor and recipient rows monotonically by `1/r` minimizes the total absolute value of
this disturbance in the unconstrained one-dimensional transport problem. Gene selection then
prioritizes retaining the template support.

## Invariants

For every emitted perturbation block:

- raw integer, finite, non-negative CSR counts;
- `X.sum(axis=1) == row_margins` exactly;
- `X.sum(axis=0) == column_totals` exactly;
- the `cell-eval2` bulk-log-normalized profile is exactly the one implied by
  `column_totals`;
- no cell exceeds the requested row margin or the VCC one-million-count cap;
- support changes are measured explicitly (`introduced_nnz`, `removed_template_nnz`);
- the control profile is a fixed point: if target column totals equal the template totals,
  the output equals the template exactly.

## Metric implications

- **PDS:** unchanged for a fixed `column_totals` vector.
- **MSE plug-in profile:** unchanged for fixed columns. The jackknife correction can change
  because cell layout changes; the corrected metric does not promise exact invariance.
- **NMAE:** can change, because it reads arithmetic mean per-cell CPM. Matching row depths and
  control depth/expression correlation is the intended benefit.
- **Fidelity/reach/Jaccard:** can change through realistic zeros, ranks, and variance. This can
  help or hurt; only leak-free local scoring can settle it.

The emitter is therefore not a free score improvement. It is a cleaner cell-level realization of
the same model pseudobulk and a plausible independent lever for the DE metrics.

## Complexity and memory

Let `n=400`, `G=18,533`, and `E` be nonzeros in the selected control templates.

- Column allocation: `O(E + G*n_fallback)`; fallback is paid only for target-positive genes absent
  from all selected templates.
- Row repair: normally `O(nG)` mask work plus the number of support-priority passes; it moves
  exactly half the pre-repair row-margin L1 error.
- Output conversion: `O(nG)`.
- Peak prototype memory: dense `int32` output `nG*4` bytes (~28.3 MiB), dense Boolean support
  `nG` bytes (~7.1 MiB), plus the selected CSR template and small vectors. A realistic block
  should remain well below 150 MiB. Only one perturbation block is held at a time.

This is deliberately sequential and bounded. It does not load a full submission, donor table, or
reference dataset into memory.

## Verification

The tiny synthetic unit test checks:

- exact identity on an unchanged control profile;
- exact row and column margins after a perturbation;
- exact scorer-style pseudobulk from those column margins;
- integer/non-negative output;
- exact preservation of the template library-size distribution;
- arithmetic-CPM relative L1 error of **1.58%** against the depth-aware template target, versus
  **12.87%** for the old-style even allocation on the same deliberately depth-confounded fixture;
- 6 within-gene count moves repaired the row margins, introducing four nonzeros outside the
  original support and removing none;
- exact margins and non-negativity over 50 seeded random sparse fixtures.

Run:

```bash
uv run python -m unittest discover \
  -s experiments/043-template-emitter -p 'test_*.py' -v
```

## Proposed bounded follow-up

Before production integration, use a small real control slice and compare only cheap diagnostics:
library-size quantiles, per-gene detection rate, arithmetic mean CPM, zero-pattern overlap, and
jackknife correction on a handful of genes/groups. A full `cell-eval2` run should happen only
after that audit, one scorer process at a time under the memory watchdog.

---

## Real-data follow-up: Orion-HCT116 × experiment 040

The bounded diagnostic is implemented in `real_data_diagnostic.py`; complete machine-readable
output is in `real_data_results.json`.

### Artifact audit and target choice

The selected experiment-040 file was:

```text
out/runs/g0.75_q1.0_a1.08_spread0.25/pred.h5ad
```

Its actual shape is 52,400 × 18,080 with 85 perturbations plus 18,400 controls—the H1-test
layout, despite experiment 040's README currently saying those files were built on HCT116.
Orion-HCT116 is 38,200 × 18,533 with 132 perturbations. The files share 18,077 genes and only
three perturbations:

| target | real HCT116 bulk-log effect RMS |
|---|---:|
| MED15 | 0.07350 |
| **MED25** | **0.07609** |
| PLAGL2 | 0.06653 |

MED25 was selected mechanically as the strongest real HCT116 pseudobulk response of the three.
The experiment remains useful as an emitter comparison because both emitted versions use MED25's
identical experiment-040 column totals. It is not a clean model-validation result across cell
contexts.

### Exact pseudobulk

| check | result |
|---|---:|
| Maximum column-total difference | **0 counts** |
| Maximum `bulk_lognorm` difference | **0.0** |
| Full bulk vector array-equal | **yes** |

The template rearrangement therefore cannot move PDS's profile. It also leaves MSE's plug-in
profile unchanged, though MSE's cell-layout jackknife can still move.

### Library sizes

Exact columns fix the grand total at 7,796,544 counts, so absolute template depths must have the
same mean as quasi. The HCT116 control depth distribution was therefore preserved up to one global
scale and integer rounding.

| source | median | CV | normalized q05 | normalized q95 | shape KS vs HCT116 |
|---|---:|---:|---:|---:|---:|
| HCT116 controls | 16,798 | 0.464 | 0.435 | 2.012 | — |
| Existing quasi | 19,497 | **0.0036** | 0.993 | 1.005 | **0.493** |
| Fixed-margin template | 18,047 | **0.461** | 0.438 | 2.004 | **0.0016** |

The quasi block is essentially constant-depth. The fixed-margin block recovers the real control
distribution's shape almost exactly.

### Detection and zeros

| diagnostic | quasi | fixed-margin |
|---|---:|---:|
| Mean detected genes/cell | 7,364 | **4,719** |
| Real MED25 target | \- | 4,427 |
| Real HCT116 control | \- | 5,084 |
| Detection-rate Spearman vs real target | 0.756 | **0.814** |
| Detection-rate MAE vs real target | 0.223 | **0.083** |
| Detection-rate Spearman vs control | 0.756 | **0.806** |
| Detection-rate MAE vs control | 0.205 | **0.072** |
| >5% detected-gene Jaccard vs real target | 0.763 | **0.821** |

The fixed-margin output retained a mean cellwise support Jaccard of 0.781 against its exact source
control cells. Row repair moved 208,628 counts (2.68% of the block), retained 1,723,741 template
nonzeros, removed 310,573, and introduced 164,728.

### Arithmetic mean CPM

| diagnostic | quasi | fixed-margin |
|---|---:|---:|
| Relative L1 vs nominal column CPM | **0.00018** | 0.03933 |
| Relative L1 vs depth-aware template CPM | 0.07306 | **0.06355** |
| Relative L1 vs real MED25 CPM | 0.92938 | **0.92873** |
| Mean-CPM Spearman vs real MED25 | 0.7583 | **0.7638** |
| LFC Spearman vs real | **0.0680** | 0.0647 |
| LFC NMAE proxy | **8.828** | 8.988 |

This exposes the tradeoff cleanly. Quasi nearly equates arithmetic CPM with count-weighted CPM by
flattening depth. The template emitter deliberately gives that up to reproduce HCT116's
depth/expression relationship. It modestly improves raw mean-profile resemblance, but slightly
worsens the perturbation LFC comparison on this one cross-context target.

### Simple DE-rank diagnostic

This is not `cell-eval2`. It uses two-sided SciPy Mann-Whitney U tests against 2,000
library-quantile-stratified real HCT116 controls, BH within the 10,134 control->5 CPM genes, and
excludes MED25 itself.

| diagnostic | quasi | fixed-margin |
|---|---:|---:|
| Predicted BH<0.05 genes | 9,263 | **7,039** |
| Real BH<0.05 genes | 551 | 551 |
| Significant-set Jaccard | 0.0546 | **0.0630** |
| Recall of real significant genes | **0.922** | 0.817 |
| Direction accuracy in significant intersection | 0.557 | **0.624** |
| Top-100 set overlap | 0.000 | **0.100** |
| Top-500 set overlap | 0.028 | **0.174** |
| Top-1,000 set overlap | 0.117 | **0.205** |
| Signed −log10(p) Spearman | **0.0319** | 0.0193 |
| Rank-biserial Spearman | **0.0024** | −0.0153 |

The fixed-margin emitter substantially improves the top-ranked DE set and sign purity while
reducing extreme over-calling. It also loses some real-DE recall and does not improve global rank
correlation. That pattern is directionally favorable for Jaccard/fidelity and ambiguous for reach.

### Resource bound

- Wall time: **5.91 seconds**
- Peak process-tree RSS: **0.98 GiB**
- Swaps during run: **0**
- Outer kill threshold: **4.0 GiB**
- No `cell-eval2`, full scorer, submission build, or external submission was run.

### Decision

The construction is credible enough for a multi-target diagnostic, but this one-target,
cross-context result is not enough to replace quasi. The highest-value next check is the same
bounded analysis over MED15 and PLAGL2, followed by a native HCT116 prediction file if experiment
040 is rebuilt correctly. A scorer run should wait until the emission direction is consistent
across targets.
