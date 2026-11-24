# The Orion HCT116 harness is depth-confounded — read this before trusting any DE result measured on it

**Found 2026-08-27.** Verified twice, independently, and reproduced by intervention.

## The defect

In `data/harness/orion-hct116.h5ad` the **control cells are 38% deeper than the
perturbation cells**, and this holds for **every single perturbation**:

| harness | control median UMI | pert median UMI | pert/ctrl | perturbations shallower than control |
| --- | ---: | ---: | ---: | ---: |
| **orion-hct116** | 16,798 | 12,178 | **0.725** | **132 / 132 = 100%** |
| h1-test *(same official geometry)* | 20,190 | 19,878 | 0.985 | 58 / 85 = 68% |

Per-perturbation ratio on orion: median `0.719`, IQR `[0.694, 0.752]`. No biology makes
132 unrelated knockdowns each cut RNA capture by the same 28%.

## Where it comes from

`scripts/build_orion_harness.py`. Its per-batch loop takes control cells greedily,
first-come-first-served, so all 18,400 controls fill from batches 1–20, while
perturbation cells (150 needed per target, only 1–3 available per batch) accumulate
across all 109 batches. Source batch depth varies a lot — median non-targeting UMI is
24,850 in batches 1–20 against 16,589 in 21–60 and 18,007 in 61–109.

**The confound is manufactured by the builder and is not in the source data.** In the raw
Orion parquets, within each batch, panel-target and non-targeting cells have the same
depth: median ratio `0.9896`, IQR `[0.976, 1.012]` across 109 batches.

## What it does to the metrics

The Wilcoxon rank-sum is driven by the per-cell zero/dropout fraction. CPM normalization
fixes the mean but **not** the zero fraction, so uniformly shallower perturbed cells look
uniformly rank-shifted *down*. The reference's DE calls therefore skew massively negative:

| reference | pert/ctrl depth | fraction of reference-significant genes with **negative** lfc |
| --- | ---: | ---: |
| **orion-hct116** | 0.725 | **0.910** [0.897, 0.920] |
| Norman 2019 K562 CRISPRa | 0.885 | 0.797 |
| **h1-test** | 0.985 | **0.538** [0.497, 0.567] |
| **Replogle 2022 K562 essential** | 0.962 | **0.497** [0.382, 0.614] |
| **Replogle 2022 RPE1** | 1.020 | **0.553** [0.494, 0.606] |

Proof by intervention — binomial thinning of raw counts, same DE re-run:

| orion condition | pert/ctrl depth | fraction negative |
| --- | ---: | ---: |
| as shipped | 0.725 | 0.910 |
| controls thinned past the perturbation depth | 1.140 | **0.275** |
| every cell thinned to exactly 7,000 UMI | 1.000 | **0.541** |

**The sign of the skew tracks the sign of the depth mismatch and inverts when you invert
it.** That is not biology.

The effect is graded, not a threshold — Norman at ratio 0.885 already sits at 0.797.

## What it does NOT depend on

- **Not the lfc formula.** Four definitions (cell_eval2's exact arithmetic-CPM ratio,
  difference of mean log1p, library-weighted pseudobulk ratio, median-centred) agree to
  within 0.002 on every dataset.
- **Not the significance threshold.** orion holds 0.910 → 0.920 → 0.933 from `padj < 0.05`
  to `< 1e-4`, and 0.914 in the top 100.
- **Not a compositional drift you can normalize away.** Median-centring changes orion by
  0.000; DESeq2-style median-of-ratios size factors make it *worse* (0.939); removing the
  shared cross-perturbation LFC axis moves it only 0.908 → 0.873. CPM *is* per-cell
  renormalization to a fixed library size and it does not help. Only count-level thinning
  to a common depth removes it.
- **Not the cell line.** HCT116 is not special. A correctly built HCT116 reference would
  land near 0.5, and any reference with unbalanced depth shows the skew whatever the line.

## Cheap diagnostic for any future harness

Compare `mean log1p difference` against the `CPM ratio` direction. On orion, the mean
log1p difference is negative for **84.6%** of *all* genes while the CPM ratio sits at
0.462 — that gap is the depth/dropout signature. On h1-test the same number is 0.542.

Two more tells on orion: pairwise LFC correlation between unrelated perturbations is
`0.145` (h1-test `0.062`), and 8 genes are called significant-down in ≥80% of
perturbations against 0 up (34 vs 2 at ≥50%); h1-test has 0 and 0.

## What this invalidates, and what it does not

**Suspect — any absolute DE reading taken on this harness**, and in particular anything
about direction, `fid`, `reach`, `jac` or `n_real`. That covers experiments 031, 061, 065,
090, 091, 092, 093 and the HCT116 half of 088.

**This is the missing explanation for the v14 transfer failure.** Experiment 093 measured
`--tail-spread 1.0 --keep-top 850` at **+0.0202 avg6** on this harness and the board
returned **−0.00291**, with `fid` arriving at the wrong sign. Its own BOARD-RESULT blamed
cells-per-perturbation. The real cause is that the harness's reference DE set is
depth-skewed, so its `n_real`, its precision curve and its direction distribution are all
wrong at once.

**Still valid — same-harness A/B comparisons where both arms face the identical
confound**, which is the design of exp090 and exp091. What they measure is which operator
predicts *this* reference better; the ranking is meaningful, the magnitudes and the
per-metric attribution are not, and transfer is not guaranteed.

**The official 2026 data looks clean.** `data/challenge/manifest.json` is 3 contexts × 300
perturbations × 400 cells with 18,400 controls and 46 NTC ids; per-NTC-id median UMI in
`context_A/B/C` varies only 1.12–1.19× max/min at exactly 400 cells per id — a
deliberately balanced build. h1-test has the identical official shape and a 0.985 ratio.
So a plan calibrated on the 0.908 number is fitting a harness bug.

The one residual risk that cannot be checked locally: if the held-out 2026 ground-truth
perturbation cells were sequenced shallower than the released controls, a weak version of
this would be real. Size any directional prior to `h1-test`'s 0.538 or RPE1's 0.553 — a
weak negative prior is defensible, `0.9` is not.

## Fix

Rebuild the harness with depth-matched pools: either draw controls from the same batches
in the same proportion as the perturbation cells, or binomially thin both pools to a
common depth before writing. Until then, **use `data/harness/h1-test.h5ad` for anything
touching the DE metrics.** exp088 §5 already established it as the reach-representative
instrument (raw reach 0.23–0.24 against the board's 0.169; orion is 0.015).

---

## Does it reach the MODEL? No — the donor tables are clean

The obvious follow-up: two of our three donors (`orion_hct116`, `orion_hek293t`) are built
from the same Orion source by `scripts/build_orion_donor.py`. Audited 2026-08-27.

**The donor builder has none of the three constructs that confound the harness.**
`build_orion_harness.py` fills its control quota first-come-first-served (`:120`), caps
perturbations at 400 (`:126`), and breaks once both quotas are met (`:137-139`) — so
controls come from the first few deep batches while perturbations accumulate across all
109. `build_orion_donor.py:70-106` is a single unconditional pass over **every** batch with
both pools accumulating from the same cells, and it normalizes per pool
(`b_ctrl = log1p(TS * ctrl_sum / ctrl_sum.sum())`, `:112`), which removes a global
library-size difference exactly.

Measured on 12 batches spanning all 109 (357,300 perturbation + 18,206 control cells):

| | control | perturbation | ratio |
| --- | ---: | ---: | ---: |
| median UMI | 18,144 | 18,088 | **1.0031** |

Control fraction is a flat 4.85% in every batch, so the large batch-to-batch depth swing
(median 13,990 → 25,758 UMI) cancels instead of accumulating.

**Direct test.** Building the pure artifact direction from HCT116 *control cells only* —
split within batch at the median library size, ratio 0.614, harsher than the harness's
0.725 — gives `‖d_depth‖ = 7.06`, as large as a real effect (median `‖delta_p‖ = 7.82`). If
the donors carried it, it would dominate. Projecting the 290/294 panel perturbations onto
it, on gated genes, using the centered `s` the model actually emits:

| donor | cos median | IQR |
| --- | ---: | ---: |
| orion_hct116 | −0.013 | (−0.187, +0.140) |
| orion_hek293t | −0.018 | (−0.229, +0.241) |

No consistent signed alignment; a confound would show one.

`delta.npy` *is* 71%/68% negative overall, but that is entirely one bin: 18,630 near-silent
gene tokens under 1 CPM where a 165,931-cell control pool detects a trace and a ~154-cell
perturbation pool does not, so `delta = 0 − ε`. Median magnitude `0.0018`, 3.7% of the
delta variance, and **every one of them is below the `--min-cpm 10` gate**. Restricted to
genes above 10 CPM the fraction negative is `0.514` (HCT116) and `0.525` (HEK293T) against
`0.510` for `k562_gwps`, which comes from an entirely different pipeline.

**Two caveats to carry forward rather than act on:**

1. The unequal pool sizes do put a real negative bias on sub-1-CPM genes, and the CPM gate
   does **not** fully remove it. The gate reads the **target context's** expression
   (`make_submission_pm_cd4.py:1178`, `b_ctrl = log1p(TS * totals / totals.sum())` from the
   context's own control cells), not the donor's — so a gene that is well expressed in
   A/B/C but near-silent in Orion passes the gate while carrying the donor's `0 − ε`.

   Sized on context A: 9,135 genes pass the 10-CPM gate, 8,962 of them are on the HCT116
   donor axis, and **466 of those (5.1% of the scored axis) sit below 1 CPM in the donor**.
   Their deltas are 76.9% negative — but with median `|δ| = 0.0025` against `0.0435` for
   well-expressed genes, i.e. **17× smaller than a real effect**, and far below the 0.10
   `|lfc|` where most of our call set already lives.

   So: real, correctly identified, and immaterial at the shipped operating point. Not worth
   a fix. It would stop being immaterial if the CPM gate were lowered — exp088b raised it
   5 → 10, moving away from the exposure — or if the Orion donors were ever weighted up
   relative to `k562_gwps`, whose axis was pre-filtered to expressed genes and has no such
   bin.
2. HEK293T's shared component `m` is the most depth-aligned of the three
   (`cos(m, d_depth) = −0.66`). Harmless today because production runs `BETA = 0.0` and
   `MultiDonorModel.effect()` emits only the centered residual — a control-pool artifact is
   identical across all of a donor's perturbations and therefore lives entirely in `m`. **If
   anyone ever turns `beta` up, that vector is partly artifact.**

---

## Repair attempts, measured (2026-08-27)

Two fixes were built and tested. Neither fully restores the file, and the honest
recommendation is unchanged: **use `h1-test`.**

| harness | depth ratio | zero-fraction gap | **frac. of reference-DE genes negative** | `n_real` median |
| --- | ---: | ---: | ---: | ---: |
| greedy (as shipped) | 0.725 | **+0.0340** | **0.897** | 322 |
| `--sampling reservoir` | 0.987 | +0.0012 | 0.873 | **3** |
| **`--sampling greedy --thin-controls`** | **1.000** | **−0.0015** | **0.630** | 74 |
| `h1-test` (clean reference) | 0.985 | +0.0019 | 0.500 | 222 |

The **zero-fraction gap** — the mean per-cell fraction of zero entries, perturbation minus
control — is the quantity the rank test actually reacts to, and it is the cleanest single
diagnostic. Both repairs drive it to h1-test's level from the shipped `+0.0340`.

**`--sampling reservoir` fixes the depth and destroys the harness.** Drawing controls
uniformly over all 109 batches does equalize depth, but Orion's batches differ enough that
a uniformly-drawn control pool is heterogeneous, which widens the reference distribution
and costs nearly all the power: `n_real` falls from 322 to **3**. Its `frac_neg` of 0.873 is
then a small-*n* artefact — with three significant genes per perturbation the set is just
the few strongest effects, and the strongest effects of a knockdown genuinely are
down-regulation. Do not use it.

**`--thin-controls` is the better repair.** It keeps the tight, single-region control pool
the greedy fill produces and binomially thins each control cell's counts so the control
depth distribution matches the perturbation pool's quantile for quantile — which is exactly
how a shallower sequencing run would have produced those cells. Depth ratio lands at
`1.0004`, the zero gap at `−0.0015`, `frac_neg` falls **0.897 → 0.630**, and `n_real` stays
usable at 74.

**A residual skew survives, and it is not depth.** `0.630` against h1-test's `0.500`. The
likely cause is the one thing thinning cannot touch: the control pool still comes from
batches 1–20 while the perturbation cells span all 109, so batch effects remain
unmatched. Fixing that needs batch-matched pools, not a depth correction — and at 150 cells
per perturbation against h1-test's 400, the file would still be the weaker instrument.

`--sampling greedy --thin-controls` is also a **drop-in replacement**: same 133 target
labels, same 18,533-gene axis, same `(38,200, 18,533)` shape, and the same perturbation
cells. The *only* difference from the shipped file is the control counts. So any result
measured on the old harness can be re-run against it as a controlled A/B with exactly one
variable changed.

### Verdict

`data/harness/orion-hct116-thinned.h5ad` is a usable file for *relative* A/B work where
both arms face the same residual, which is what experiments 090 and 091 were. It is still
**not** a place to read an absolute `fid`, `reach`, `jac`, `n_real` or direction number.
For those, use `data/harness/h1-test.h5ad`, and gate any harness with
`scripts/check_harness_depth.py` before trusting it.

---

## The repair is validated against the board

Experiment 091's pruned basis is the one change measured on **all three** surfaces, so it is
the test of whether the repair improved the instrument:

| | predicted Δ avg6 | reach sign | over-statement |
| --- | ---: | :---: | ---: |
| **broken** harness (`orion-hct116.h5ad`) | `+0.000876` | **up (wrong)** | **7.0×** |
| **repaired** harness (`--thin-controls`) | `+0.0002` | **down (right)** | **1.6×** |
| **the board** | `+0.0001257` | down | — |

Per-metric on the repaired harness — `nmae` up, `reach` down, `fid` up, `jac` up, `pds`
flat-to-down — matches the board's signs on **four of five** members. The broken harness had
predicted the gain would land in `nmae`, `reach` **and** `fid` with `pds` flat, and `reach`
went the other way.

So depth-matching the control pool took a 7× over-statement with an inverted `reach` sign
down to a 1.6× over-statement with the signs right. **The repair works**, and
`orion-hct116-thinned.h5ad` is a usable instrument for this class of A/B.

*Caveat on the arithmetic:* the two arms are separate `run.py` invocations rather than a
paired run. The pipeline is deterministic at a fixed seed, so the subtraction is sound here
— but experiment 091's own RESULTS warns that cross-run differencing on this harness once
produced a spurious `+0.046`, so keep arms in one invocation where possible.

**This does not lift the `n_real` limitation.** The repaired harness's `n_real` is 68
against the challenge's thousands, so `fid`'s coverage term, `jac`'s union and `reach`'s
normalisation are still unreadable there — that is a geometry problem, not a depth problem,
and it is what cost a8 its `jac`.
