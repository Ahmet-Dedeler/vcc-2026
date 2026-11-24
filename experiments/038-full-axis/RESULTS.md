# 038 — The operator was reaching 41% of the scored genes. Fixing that is worth +0.022.

**Date:** 2026-08-22 · LODO built on the **official 18,533-gene axis**, as the submission is,
rather than on the donor intersection as every previous fold was. 259 VCC-panel perturbations.

## The observation that had been sitting in the notebook for a day

Experiment 025 audited the gene axis and reported: at the CPM>5 gate the scorer tests 9,923
genes in context A, of which **6,729 are on the whitening operator's axis and 3,131 are not**.
That was filed as a *scaling* problem — global renormalisation was silently down-weighting the
off-axis genes — and fixed with `--whiten-renorm preserve` for +0.006.

The larger question went unasked for a day: **why are a third of the tested genes not being
whitened at all**, under the single largest lever in the project?

The cause is the operator's axis. It is fitted on genes **every** donor measures, and K562-gwps
carries only 7,679 of the official axis where each Orion line carries 18,111. Intersecting three
donors to build a covariance throws away everything two of them can see.

## Result

Folds holding out an Orion line, where the scored space is 17,806 genes and the intersection
operator reaches 7,312 of them — **41%**:

| operator | hct116 | hek293t | macro |
|---|---|---|---|
| none | 0.7030 | 0.7061 | 0.7045 |
| **intersection, gamma 0.85** (what v6/v8 ship) | 0.7728 | 0.7573 | 0.7650 |
| **union, gamma 0.65** | **0.8019** | **0.7716** | **0.7868** |

**+0.0218 macro PDS**, replicated on both folds (+0.029, +0.014). For scale, whitening itself
was +0.070 and the emission lever was +0.0008 on the board.

The gamma optimum moves **0.85 → 0.65**, which is the expected direction: the union covariance is
estimated over 2.4x more genes from the same number of perturbations, so it is less
well-conditioned and tolerates less aggressive inversion.

| gamma | 0.3 | 0.4 | 0.5 | 0.6 | **0.7** | 0.85 |
|---|---|---|---|---|---|---|
| macro | 0.7593 | 0.7727 | 0.7818 | 0.7868 | **0.7875** | 0.7809 |

Broad plateau at 0.6–0.7; **0.65 chosen** as its middle.

## Why the earlier folds could not see this

Every LODO fold before this one defined the feature space as the **donor intersection**, so
"genes the operator cannot reach" was empty by construction. The harness had the blind spot
built into it, and the production pipeline did not. That is the same class of error as the
`gene_norm`-on-the-truth bug caught in experiment 017 — a harness that quietly answers an easier
question than the one being asked.

## Cost, and why it was not tried sooner

An 17,806 x 17,806 eigendecomposition rather than 7,312 x 7,312: ~15x the flops, ~15 minutes and
several GB. Cheap in absolute terms, and the reason it went untested was that no fold had ever
been built on an axis where it would matter.

## Shipped in v9

`--whiten-axis union --whiten-gamma 0.65`, operator fitted on 17,942 genes (union minus panel
targets) from 33,787 donor perturbations. Emitted magnitude preserved: overall effect norm
4.801 → 4.784 against v8's 4.766, so alpha keeps its meaning and the change is isolated to the
operator.

Expected board effect at the LODO→board ratio calibrated on v3→v6 (1.26): raw `pds`
0.8163 → ~0.843, scaled 0.701 → ~0.757, **≈ +0.009 overall**.

---

# 038b — H1 count-level result (quarantined 2026-08-22)

**Correction:** this comparison used a mismatched pre-whitening exclusion mask.
The runner removed the 300 VCC-2026 target genes rather than the 85 targets in
the H1 reference. Because whitening is a dense transform, 80 H1 self-target
coordinates could mix into scored coordinates before `cell-eval2` removed the
target columns. H1 truth did not enter the fit, so this is not truth leakage,
and code inspection alone does not establish the direction of bias. Treat the
numbers below as unfaithful until rerun with the actual H1 targets excluded.

Experiment 038's Orion LODO result above is unaffected and remains the evidence
for v9.

`h1-test` (85 perturbations, 400 cells/pert, 18,400 controls) is the only harness matching the
competition on both cell counts, and it has one calibration point behind it: it predicted v8's
board gain at +0.0008 against +0.00082 delivered.

| | v8 (intersection, gamma 0.85) | **v9 (union, gamma 0.65)** | delta |
|---|---|---|---|
| **avg6** | 0.1791 | **0.1895** | **+0.0104** |
| `pds` raw | 0.8552 | **0.8655** | +0.0103 |
| `nmae` raw | 0.9619 | **0.9561** | −0.0058 (better) |
| `reach` raw | 0.2105 | **0.2373** | +0.0268 (+12.7%) |
| `fid` raw | 0.5057 | 0.5061 | +0.0004 |
| `jac` raw | 0.0644 | 0.0638 | −0.0006 |

**Five of six metrics improve and the sixth moves by 0.0006.** `nmae` improving is the notable
one — whitening trades cosine for discrimination and `nmae` rewards cosine, so every previous
increase in whitening cost it. Extending the *axis* rather than the *strength* apparently does
not: reaching 17,942 genes at gamma 0.65 beats reaching 7,312 at gamma 0.85 on both.

Projected onto the board metric by metric: `pds` +0.020 scaled, `reach` +0.030, `nmae` +0.010,
`jac` −0.002 → **≈ +0.010 overall, to ~0.148**. The board currently reads Ibrahim Mansour
0.1448, VIVAI 0.1412, us 0.1377, over 127 teams.

Two caveats kept explicit. The single calibration point behind `h1-test` came from a change that
moved only `fid`/`reach`/`jac`; v9 moves `pds` and `nmae` as well, and those may not translate
at the same ratio. And LODO put the `pds` gain at +0.022 macro where `h1-test` reads +0.0103 raw
— the two disagree by roughly 2x on the largest component, so the projection is a range, not a
number.
