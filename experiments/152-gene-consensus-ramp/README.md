# Experiment 152 — gene-wise consensus-ramped dual moment

Status: **preregistered before Experiments 150/151 outcome access**

Date: **2026-08-31**

## Question

Experiments 148 and 149 show that changing the arithmetic-CPM moment while
preserving every pooled gene total can improve the DE metrics with PDS exactly
unchanged. Experiment 151 asks whether one lambda per target should be
confidence-gated. This experiment asks the finer question: should the
gamma-0.45 relaxation be applied strongly only to target-gene coordinates that
agree across independent source-only HEK293T and HCT116 predictions?

## Frozen transform

For each target shared by the two source builds, compute the gamma relaxation
in each context and the minimum absolute relaxation magnitude. Within each
target, rank that magnitude by exact midrank percentile.

- disagreeing or lower-half coordinates receive lambda `0.5`;
- upper-half sign-agreeing coordinates receive a smooth lambda ramp from
  `0.5` to `1.0`;
- ineligible coordinates receive lambda `0.5`;
- the 12 HCT-only targets retain uniform lambda `1.0` and are excluded from
  the primary analysis.

The gene-wise desired arithmetic CPM is renormalized to one million and passed
through the existing bounded exact-margin allocator. Every incumbent column
total and scorer pseudobulk must remain exactly unchanged.

The immutable specification is `PREREGISTRATION.v1.json`.

This experiment authorizes neither a challenge build nor a submission.
