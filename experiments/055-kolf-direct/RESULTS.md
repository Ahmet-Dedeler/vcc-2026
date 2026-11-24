# 055 — Direct KOLF2.1J donor transfer results

Date: 2026-08-22

## Outcome

The compact KOLF2.1J donor is real and modestly positive, but not yet strong
enough to alter v9 without a whitening-aware confirmation.

The fixed rule selected using only held HCT116 and HEK293T folds was:

```text
KOLF representation: centered raw quantized z-score
KOLF blend fraction: 0.10
```

| evaluation | baseline PDS | + KOLF PDS | delta |
|---|---:|---:|---:|
| held HCT116 | 0.683473 | 0.682322 | -0.001151 |
| held HEK293T | 0.695868 | 0.697995 | +0.002127 |
| **development mean** | — | — | **+0.000488** |
| frozen K562 target-role holdout | 0.736349 | 0.740342 | **+0.003992** |
| **three-fold macro** | — | — | **+0.001656** |

The untouched K562 result is encouraging, but the negative HCT116 fold and
small development mean make this a weak positive rather than a shipping
decision. The next gate is a production-shaped union-whitening test using the
fixed raw/0.10 rule.

## Protocol

- Strict HTTP range requests fetched exactly 273 VCC-panel rows from the
  official KOLF browser matrix: 10,539,438 bytes, not any complete source
  chunk.
- Compact donor: 273 perturbations x 18,106 mapped VCC genes.
- Evaluation intersection: 240 perturbations x 7,312 non-panel genes.
- Every query was ranked against the complete 240-target gallery.
- All 300 VCC target coordinates were removed before centering,
  normalization, and retrieval.
- Candidate selection used only HCT116 and HEK293T. K562 was evaluated once
  after the mode and fraction were frozen.
- Raw and mean-absolute gene-normalized KOLF effects were tested at fractions
  0.10, 0.25, and 0.50.

The source is pinned to repository commit
`8f84531a784f24b14d7e2c9b80c375d7e78fd46c`. Its metadata, matrix builder,
browser loader, and license were SHA-256 checked before decoding. The decoded
effect is exactly `int8_value * 0.04`.

## Memory and validation

- Extraction produced a 30 MB compact donor.
- Five synthetic checks passed.
- Peak evaluator RSS: 2.79 GiB.
- No process remains running.

Because this exceeds the new 1 GiB local-experiment ceiling, do not rerun the
evaluation locally. Any whitening-aware follow-up must run on AWS.

## Decision

**Advance KOLF to experiment 057; do not modify or rebuild v9 from this result
alone.** The direct gain is smaller than the strongest established geometry
gain and has one negative development fold. It must survive union whitening
on the frozen K562 role before it is considered for a later candidate.
