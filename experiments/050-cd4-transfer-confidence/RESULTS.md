# 050 — CD4 effect breadth predicts transfer quality, not condition choice

**Date:** 2026-08-22. Official primary-CD4 supplementary data, restricted to
the 87 VCC-panel targets with direct CD4-versus-K562 comparisons in all three
conditions.

## Result

Every target-condition comparison had positive cross-line effect correlation.
Mean Pearson correlation was 0.339 in Rest, 0.357 in Stim8hr, and 0.355 in
Stim48hr.

The strongest CD4-only confidence feature was response breadth:

| condition | Spearman(`log1p(n DE genes)`, K562 transfer correlation) |
|---|---:|
| Rest | 0.500 |
| Stim8hr | 0.582 |
| Stim48hr | **0.615** |

All three associations are highly significant (`p <= 8.1e-7`). The atlas'
published within-CD4 donor-correlation field was weak on the rows where it was
available (`rho = 0.02–0.21`), and target cell count was uninformative.

## What breadth does not do

It does not select the best condition for an individual target. Choosing the
condition with the most DE genes gave mean transfer correlation 0.3562 versus
0.3572 for fixed Stim8hr:

`delta = -0.00105`, bootstrap 95% interval `[-0.00606, +0.00416]`.

Transferability is already nearly rank-identical across conditions:

| pair | Spearman |
|---|---:|
| Rest vs Stim8hr | 0.9535 |
| Rest vs Stim48hr | 0.9530 |
| Stim8hr vs Stim48hr | 0.9534 |

The oracle condition is worth only +0.0141 mean correlation over fixed
Stim8hr, so condition routing has modest headroom even with truth.

## Decision

After the compact CD4 effects are extracted:

1. Test fixed Stim8hr and the three-condition mean first.
2. If direct CD4 blending helps, test **blend strength as a smooth monotone
   function of CD4 effect breadth**.
3. Do not use effect breadth to choose one condition per target.
4. Do not use the sparse `donor_correlation_mean` field as the primary
   confidence signal.

The K562 correlations are diagnostic truth only. Any production weight must be
computed from CD4 metadata and selected without challenge outcomes.

## Coverage value independent of lineage

The current three donor tables cover all 300 targets, but only 259 are observed
in all three donors. Adding CD4 changes the donor-count distribution from:

```text
1 donor: 6, 2 donors: 35, 3 donors: 259
```

to:

```text
1 donor: 1, 2 donors: 9, 3 donors: 35, 4 donors: 255
```

CD4 directly measures 27 panel targets absent from K562-GWPS. It therefore
adds useful independent coverage even if its lineage-specific advantage for
context A turns out to be weak.
