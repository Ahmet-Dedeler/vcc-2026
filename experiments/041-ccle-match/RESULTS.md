# 041 — CCLE matching identifies the three validation contexts

Date: 2026-08-22

## Result

Cross-platform matching against the 1,600-sample DepMap Public 25Q3 CCLE
expression mirror gives the same top identity for every feature/rank variant:

| challenge context | consensus top match | lineage |
|---|---|---|
| **A** | **JURKAT** | lymphoid / T-ALL |
| **B** | **HELA** | cervix |
| **C** | **CAL33** | head and neck squamous |

The runner compared Pearson and Spearman similarity on all common genes and
on 2,000, 5,000, and 10,000 high-variance CCLE genes. The top match was
unchanged across all eight methods.

The nearest alternatives support the same tissue calls:

- A: PF382, KE37, MOLT3, and P12-ICHIKAWA, all T-lymphoblastic lines;
- B: SIHA and SW756, both cervix lines;
- C: CAL27 and several other head-and-neck squamous lines.

## Cross-platform validation

The UMI-control-to-bulk matcher was calibrated on three known controls before
interpreting A/B/C:

| method | K562 rank | HCT116 rank | independent Jurkat rank |
|---|---:|---:|---:|
| HVG-2k Pearson | 1 | 1 | 2 |
| HVG-2k Spearman | 1 | 1 | 3 |
| HVG-5k Pearson | 2 | 1 | 1 |
| HVG-5k Spearman | 2 | 1 | 1 |
| HVG-10k Spearman | 2 | 1 | 1 |

Whole-transcriptome methods were less platform-robust, especially for
Jurkat, and were not used as the primary validation. The consensus identity
nevertheless remained the same under those methods.

## Interpretation and boundary

The result is strong evidence for lineage and likely exact line identity. It
does not by itself authorize a hand-written context-specific prediction.
Public data may condition a learned or predeclared model, but the challenge
rules prohibit revising outputs from literature values.

The immediate scientific use is therefore a predeclared proxy test:

- primary-CD4 effects are a biologically matched auxiliary donor for A/Jurkat;
- the already-negative global HCT116 CD4 test does not answer whether CD4
  helps a hematopoietic target under union whitening;
- experiment 060 tests that exact question on held K562 without using any
  challenge outcome.

No production model or submission changed from the CCLE match alone.
