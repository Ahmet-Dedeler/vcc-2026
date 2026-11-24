# 285 — PREREGISTRATION: the apples-to-apples reliability benchmark

**Written 2026-09-06, hashed into `work/PREREGISTRATION.sha256` before any expression value was read.**
What had been read at the hash: file *metadata only* — `obs` / `var` / `X`-layout of the local h5ad files
(`work/inv285.json`, `work/inv285.log`), the cached JSON summaries of 267 / 268 / 273 / 280, and the `obs` +
`X.indptr` of the one remote CD4 object probed (`work/probe285-D2_Stim48hr.json`). **No count, pseudobulk,
profile, cosine or `pds` had been computed on any source at the moment of hashing.**

Measurement only. **No model, no build, no submission, no `vcc`, no `cell-eval2`, no price, no model name.**

---

## 0. The question, and the unit

Two measurements of the same perturbation disagree. **Which source of that disagreement dominates at the
challenge's budget, and which part of it can be removed by measuring more?** The deliverable is a decomposition
per source, on MATCHED targets and MATCHED cells per side, in two readouts, with the removable terms
(cell sampling, UMI depth) reported apart from the irreducible ones (guide construct, donor, batch, study, line).

**The unit that is resampled for every interval is the TARGET** (paired target bootstrap). The unit that
defines a *view* is stated per row: cells, UMIs, guide construct, donor, gemgroup/lane/batch, study, line.
Nothing is ever pooled across sources.

## 1. The response object (frozen)

For a view V of target t with summed integer counts `S` over a gene axis `G`:

```
p_V(t)  = log1p( 50,000 · S / ΣS )                     # the scorer's bulk profile, TS = 50,000
δ_V(t)  = p_V(t) − p_V(controls of the SAME view)      # the view's own controls, where controls exist
r_V(t)  = δ_V(t) − mean_{t' ∈ T} δ_V(t')               # 159's panel-centred residual, over the set's targets T
```

then every target-gene column of `T` is zeroed and the cosine is taken on the remaining columns of the shared
gene axis. **Note (declared before measuring): `r_V` is algebraically independent of the control profile**,
because the control term is identical for every target of a view and is removed by the centring — so the
residual-cosine rows do not depend on whether the two views share a control pool. The `pds` rows *do* depend on
it, and use each view's own controls where they exist and a single fixed control profile (named per source)
where they do not.

**Two readouts per cell.**

* **`cos_resid`** — the per-target cosine `cos(r_A(t), r_B(t))`; medians, p10/p90, share ≤ 0.
* **the scorer's own statistic** — `resid218._per_target` imported verbatim (273's KA-2 proved it equals
  cell-eval2's `pds_cosine` to 5.6e-17): view **A as the prediction**, view **B as the truth**, **every target
  of that source** as a distractor, midrank ties, denominator n − 1, on the un-centred `δ`; plus `cos_agg`
  (the energy-weighted aggregate cosine). Reported as `pds(A→B)` averaged with `pds(B→A)`.

## 2. The units (each a separate row; never pooled across sources)

| id | unit | how a view is formed |
|---|---|---|
| **(a)** | same guide(s), **random cell halves** at matched n | a fixed random ranking of a target's cells; side A = ranks [0, n), side B = [n, 2n). **8 independent rankings (draws), averaged per target.** |
| **(a′)** | same cells, **thinned counts** | at fixed n, each side's summed counts thinned independently by the binomial at fraction f ∈ {1, ½, ¼, ⅛} (280's `lib280.thin`, seeded) |
| **(b)** | **different guide constructs**, same target, same line/donor | side A = construct 1's cells, side B = construct 2's, sub-sampled to matched n |
| **(c)** | **different donors**, same guide/target | CD4 only (4 donors) |
| **(d)** | **different batches** — gemgroup / lane / Orion parquet batch | parity of the sorted batch index, sub-sampled to matched n |
| **(e)** | **different study, same line** | K562-GWPS vs K562-essential (Replogle 2022, two independent experiments, same line), on their 2,053 shared perturbations |
| **(f)** | **different lines** | one reference row from 276's cross-line conservation statistic (no new measurement) |
| **(s)** | **different activation state**, same donor+guide | CD4 bank, Rest / Stim8hr / Stim48hr (282's deconfounded contrast) |

## 3. The sources and the rungs (frozen from `work/inv285.json`; a target enters rung n only if it has ≥ 2n cells)

| source | genes | targets | cells/target | rungs n | units available |
|---|---:|---:|---:|---|---|
| **h1-test** (H1 hESC, 10x Flex, `data/harness/h1-test.h5ad`) | 18,080 | 85 | 400 (exact) | **25, 50, 100, 200** | a, a′, b (5 two-construct targets), d (48 lanes) |
| **K562-GWPS** (`ReplogleWeissman2022_K562_gwps.h5ad`, dense f32) | 8,248 | 272 panel (258 / 223 / 95 at n = 25 / 50 / 100) | median 171.5 | **25, 50, 100** | a, a′, b (8 two-construct panel targets), d (267 gemgroups), e |
| **K562-essential** (`…_K562_essential.h5ad`) | 8,563 | 2,053 shared with GWPS, **0 on the panel** | median 121 | 25, 50 | a, d (48 gemgroups), **e** |
| **CD4 GWCD4i** — 280's bank (guide × donor × condition pseudobulks) | 18,129 | 297 panel present | median 54 per row | matched by `n_cells` | b, c, s, a′ |
| **CD4 GWCD4i** — the bounded cell-level extraction (§6), if admissible | 18,130 | 290 panel present in D2_Stim48hr | median 99.5 | 25, (50) | **a**, a′, b, d (24 lanes), and with a second donor **c at cell level** |
| **Broad-16** (268's `work/halves/<line>.npz`, 12 lines) | 34,648 | 98 per line | ~55 | ~25 (fixed by 268's halves) | a, b, f |
| **Nadig 2025** Jurkat / HepG2 | 8,882 / 9,624 | 2,393, **0 on the panel** | 83 / 45 | 25 (50 on Jurkat) | a, b (135 multi-guide), d (55/56 gem groups) |
| **Orion** HCT116 / HEK293T (267's cached half-banks) | 18,106 | 290 / 294 | 70–95 per half | cached only | a, d (reference rows) |

**Cohort control (the correction the third review demands).** Every ladder is reported TWICE: (i) the
all-available cohort at each rung, and (ii) a **fixed cohort** — the targets that reach the largest rung of that
source, measured at every rung. **The fixed-cohort ladder is the primary**; the all-available one is reported to
show the selection.

## 4. The frozen contrasts

1. **The ladder (a)** per source, on the fixed cohort, at every rung, plus its Spearman–Brown extrapolation
   `r(kn) = k·r(n) / (1 + (k−1)·r(n))` computed **from the smallest rung** and **from the largest rung**, and the
   measured value at the intermediate rungs. Discrepancy = measured − SB.
2. **Depth vs cells at a fixed UMI budget.** At equal total UMIs per side, compare `(n = 2m, f = ½)` against
   `(n = m, f = 1)`. If the two agree, only total UMIs matter; if the first is higher, cells carry information
   beyond UMIs. Frozen pairs: (50, ½) vs (25, 1); (100, ½) vs (50, 1); (200, ½) vs (100, 1) where the rung exists.
3. **The unit contrasts at matched n**, on the SAME targets where the design allows (h1-test's 5 and K562's 8
   two-construct targets, CD4's two-guide targets), otherwise on the source's own cohort reweighted to a common
   **response-magnitude stratum histogram** (deciles of the full-depth residual row norm).
4. **The n at which (a) reaches (b) / (c) / (d) / (e)** — read off the measured ladder, and off its SB
   extrapolation where the rung does not exist. Reported per source.

## 5. Known-answer rows — the run aborts and reports FAIL if any of these misses

| row | what | required |
|---|---|---|
| **KA-1** | 267's batch-split residual-cosine medians read from its cached `work/half/truth-<fold>.npz` `rel|*`: +0.0471 / +0.0276 / +0.0282 (K562-GWPS / HCT116 / HEK293T) | **0.0e+00** |
| **KA-2** | 273's h1-test **lane-parity 200/200 residual median +0.2429**, RECOMPUTED here from `data/harness/h1-test.h5ad` with this experiment's own code | \|Δ\| ≤ 1e-4 vs `reliability273.json` |
| **KA-2b** | 159's `resid_disjoint` at 200 cells/side = **0.330039**, read from 159's stored rows | **0.0e+00** |
| **KA-3** | 280's Stim48hr **donor-half median +0.0539** recomputed from `bank280_X.npy` through `lib280` | \|Δ\| ≤ 1e-4 |
| **KA-4** | the union of two halves reproduces the source's full row: `sum(A) + sum(B) == sum(full)` on integer counts, every source I split (267 KA-1's identity) | max \|Δ\| **= 0** |
| **KA-5** | `pds` of a view against ITSELF is 1.0 for every target, and `resid218._per_target` is imported verbatim (no re-implementation) | exactly 1.0 |
| **KA-6** | 268's halves satisfy `A + B == full` and `v1 + v2 ≤ full` per line | max \|Δ\| = 0 / no violation |
| **KA-7** (extraction only) | the extracted cells of one (target, guide) of D2_Stim48hr, summed, equal the corresponding row of 280's `bank280_X.npy` **exactly** | max \|Δ\| = 0 |

## 6. The bounded CD4 cell-level extraction — the decision rule, frozen

Probed object: **`D2_Stim48hr.assigned_guide.h5ad`**, 137,554,355,000 B. Measured layout (metadata only,
`work/probe285-D2_Stim48hr.json`): `X` is **`csr_matrix`, uncompressed** (`data` float32, `indices` int64,
chunked (7245,), compression `None`), 2,863,571 cells × 18,130 genes, 11,417,773,379 nnz; `obs` carries
`perturbed_gene_name` (12,084 labels, **290 of the 300 panel present, 48,055 panel cells, median 99.5/target**),
`guide_id`, `guide_group`, `lane_id` (24), `total_counts`, `low_quality`.

**PROCEED only if** `X` is uncompressed CSR with range-addressable rows **and** the estimated contiguous-span
bytes for the planned rows are **≤ 25 GB**. Rows planned: every cell whose `perturbed_gene_name` is one of the
300 panel targets (capped at 900 cells per target), plus a non-targeting sample of ≤ 2× those cells if
non-targeting cells exist in this object (the probe finds **none** under that label, so the control profile comes
from 280's bank NTC pseudobulks for the same donor × condition, and the residual rows do not need it — §1).
Kept bank: **sparse** (int32 indices + float32 data) ≤ 3 GB, sha256 recorded; parts deleted as read.

**A second donor (`D3_Stim48hr`, 137.8 GB) only if** the first came in under 15 GB transferred **and** ≥ 2 hours
of the 7-hour box remain. It buys the one thing nothing else on disk has: a **cell-level donor contrast at
matched n against a cell-split of the same cells**.

**Otherwise** stop at the plan and use the bank-only CD4 units (b, c, s, a′).

Recorded either way: bytes transferred, number of range requests, sha256 of what is kept, wall clock.

## 7. Analysis rules frozen before measurement

* **Intervals**: paired target bootstrap, **20,000 draws, seed 20260906**, resampling targets within a source.
  Contrasts between two units on the same targets are paired; contrasts on different cohorts are reported as
  stratum-reweighted differences with the same bootstrap over the shared stratum histogram.
* **Sensitivity**: leave-one-donor-out (CD4), leave-one-lane-out (h1-test), leave-one-gemgroup-out (K562),
  leave-one-line-out (Broad-16) on the headline rows.
* **Multiplicity**: the number of (source × unit × rung) cells screened is reported in RESULTS; no cell is
  selected after a look, and no rung or unit is added after the first number is read.
* **Missing/degenerate**: a target with a zero-norm view on either side is dropped from that cell and counted.
* **Reporting**: medians (the distribution is skewed and 267/273/280 all report medians), share ≤ 0, and the
  mean where a mean is what the scorer averages (`pds`).

## 8. What would make this a negative result

If every unit reads the same at matched n — i.e. if the guide construct, the donor, the batch and the study are
indistinguishable from a cell split of the same size — then **there is no irreducible term to find at these
depths and the whole disagreement is sampling**, and that is the result. If instead the removable terms saturate
below the irreducible ones, the benchmark states the n at which more cells of one construct stop helping.
Either outcome is written up. **No arm, no dose, no price, and nothing is proposed for a slot.**
