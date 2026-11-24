# 010 — Leave-one-donor-out validation, scored as PDS

**Date:** 2026-08-21 · 3 donor lines, 259 VCC panel perturbations present in all three.
Each line held out in turn, predicted from the other two, scored with `pds_cosine` implemented
to the vcc2026 spec (panel-scope target-gene exclusion, midrank ties), macro-averaged.

## Why this replaced the previous harness

The old stack tuned on K562 → H1 hESC and scored **mean cosine**. Both were wrong:

- **Wrong biology.** H1 is embryonic; the 2026 contexts are somatic lines from different
  tissues. A +20% transfer cosine on H1 produced a *negative* pds move on A/B/C.
- **Wrong objective.** Mean cosine and 300-way PDS ranking disagree. Averaging donors smooths
  every prediction toward the shared response, raising mean alignment while **reducing
  discrimination** — which is exactly what PDS measures.

Holding out a donor line is somatic-to-somatic transfer across tissues, structurally the same
task as A/B/C.

## The harness predicted the v3 regression

| | LODO | competition |
|---|---|---|
| v1 family (gene_norm, single donor) | 0.7245 | 0.7383 |
| v3 family (EB, 3:1 multi) | 0.7141 | 0.7294 |
| **Δ** | **−0.0104** | **−0.0089** |

Right sign, right magnitude, from data that never saw the answer. Had this existed a day
earlier, v3 would not have been submitted.

## Results (macro PDS over three held-out lines)

| primary share | raw | **gene_norm** | eb |
|---|---|---|---|
| 1:0 (single donor) | 0.7184 | 0.7246 | 0.7183 |
| 20:1 | 0.7262 | 0.7328 | 0.7271 |
| 12:1 | 0.7298 | 0.7364 | 0.7269 |
| 8:1 | 0.7333 | 0.7395 | 0.7259 |
| 6:1 | 0.7349 | 0.7408 | 0.7243 |
| **4:1** | 0.7337 | **0.7412** | 0.7191 |
| 3:1 | 0.7311 | 0.7400 | 0.7141 |
| 2:1 | 0.7243 | 0.7351 | 0.7075 |
| 1:1 | 0.7133 | 0.7264 | 0.6970 |

**1. `gene_norm` > raw > EB at every weight.** v1's choice was right. EB shrinkage is a
regression, not something needing tuning — the +10.6% cosine it bought was measuring the wrong
objective, and it cost pds *and* 87% of v3's total loss via `nmae`.

**2. The weight curve is single-peaked with a broad plateau at 3:1–8:1** (0.7400–0.7412), i.e.
the primary donor should carry **75–86%** of the weight. Both ends lose: 1:1 costs 0.015 (the
smoothing effect), 20:1 costs 0.008 (throws away the second donor).

**3. Multi-donor beats single donor by +0.0166** (0.7412 vs 0.7246) — real, but only with a
dominant primary. Naive averaging (1:1) is *worse than one donor*.

## Negative result: basal-profile similarity does not pick donors

| temperature | raw | gene_norm | eb |
|---|---|---|---|
| basal5 | 0.7114 | 0.7243 | 0.6944 |
| basal20 | 0.7057 | 0.7183 | 0.6878 |
| basal50 | 0.6938 | 0.7068 | 0.6764 |
| basal100 | 0.6822 | 0.6951 | 0.6665 |

**Every basal-weighted configuration loses to every fixed weight, and monotonically worse with
sharper weighting.** The per-fold detail shows why: when K562 is held out, basal weighting is the
single best config for that fold (0.7491) — but for the other two folds it is among the worst.
It concentrates weight confidently on the wrong donor.

**This matters well beyond donor weighting.** The basal control profile is the *only* signal the
Challenge gives us about a context's identity. This says it does not predict which donor's
perturbation responses will transfer there. Any scheme that tries to adapt to A/B/C by matching
their basal profiles to a reference line should expect the same failure.

## Chosen configuration for the next submission

`gene_norm`, no EB, K562-gwps primary at **80% of weight** (8:1:1 across three donors, matching
the 4:1 two-donor optimum in primary share), α re-derived to hold emitted magnitude identical to
v1 so `nmae`/`fid` stay attributable.

## Engineering note

The first implementation used `list(map(str, genes)).index(g)` inside a comprehension — it
rebuilt a 38,606-element list and linear-scanned it per gene, ~270M operations per fold. One run
took over ten minutes and produced no output. With a precomputed dict the same sweep is **11
seconds**, which is what made the exhaustive grid above affordable.
