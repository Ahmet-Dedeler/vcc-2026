# 036 — v8 candidate: v6 + exact-total quasi spread 0.25. Do not submit yet.

**Date:** 2026-08-22 · artefact `submissions/v8.h5ad` (20.35 GB), not packaged, not submitted.

## What it is

v6 config in every modelling knob, one emission change:

| knob | v6 (live, 0.1369) | v8 |
|---|---|---|
| donors | 1:1:1 | same |
| whitening | γ=0.85, preserve | same |
| α | 1.0796 | same |
| CPM gate | 5 | 5 |
| emit | even | **quasi, spread=0.25** |

RNG split: `pick_rng = seed`, `emit_rng = seed+1000003`, so library-size sampling matches v6
exactly. First v8 build shared one stream; quasi consumed extra randoms; picks diverged;
mean profile moved. Rebuilt after the split.

## Verification (this session)

Context A, first 8 perturbations, gene-panel totals v6 vs v8:

**max |Δ| = 0, genes differing = 0, L1 relative = 0.** Bit-identical pseudobulk. nnz dropped
2.675e9 → 2.540e9 because scatter creates more zeros; that is the spread, not a mean shift.

h1-fast (exp 034b), spread 0 vs 0.25:

| | pds | nmae | reach | fid | jac | avg6 |
|---|---|---|---|---|---|---|
| 0 | 0.8506 | 0.9450 | 0.1669 | 0.5135 | 0.0730 | 0.1824 |
| 0.25 | 0.8506 | 0.9449 | 0.1901 | 0.5180 | 0.0694 | 0.1873 |

Claimed board translation ≈ **+0.007 overall**, downside bounded to jac (worth 0.002 to us).

## Why we are waiting

1. Rank uses the **latest** entry. A worse v8 drops us off rank 3. One slot left until 00:00 UTC.
2. Somatic HCT116 confirmation is **in flight** (`out_som_spread`, 4 cell-eval2 workers):
   spread 0 / 0.25 / 0.25+cpm20 / 0.4. Trust this for *which genes carry signal*; trust
   h1-fast for fid yield (400 ref cells).
3. Do **not** bundle spread with CPM gate 20. Gate 20 is a separate +0.003 claim from
   somatic-only data. Ship one variable.
4. Do not ship v7 (sharpening). Somatic nmae 0.96 → 1.26, avg6 nearly halved.

## Submit rule

Submit v8 as `"v8"` only if `out_som_spread` shows spread 0.25 with **pds within seed noise of
spread 0** and **reach up**. If pds moves, do not ship. If only cpm20 wins, that is a
different candidate and still one variable.
