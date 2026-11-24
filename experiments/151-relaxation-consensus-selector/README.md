# Experiment 151 — source-only relaxation-consensus selector

Status: **preregistered before Experiment 150 outcome access**

Date: **2026-08-31**

## Question

The fixed dual-moment emitter improved total projected score on both H1 and
held HEK293T, but different DE metrics carried the gain. This experiment asks
whether a tiny source-only confidence rule can make that gain more robust on
held HCT116 without changing any perturbation pseudobulk or PDS value.

## Frozen selector

For each target shared by the held-HEK and held-HCT source-only builds, compare
the relaxation vectors

```text
direction(gamma=0.45) - incumbent(gamma=0.65)
```

gene by gene. Weight each sign-agreement indicator by the smaller absolute
relaxation magnitude in the two contexts. Rank targets by this weighted
agreement, breaking ties alphabetically.

- top 60 of 120 shared targets: dual-moment lambda `1.0`;
- bottom 60: lambda `0.5`;
- 12 HCT-only targets: lambda `1.0`, excluded from the primary comparison.

Everything else is identical to Experiment 150. Both lambda values preserve
the incumbent target-wise column totals and scorer pseudobulk exactly.

The immutable specification is
`PREREGISTRATION.v1.json`; its digest is recorded in
`PREREGISTRATION.v1.sha256`.

## Decision rule

The gated candidate must first pass the same score-level safety gate as
Experiment 150 on the 120 common targets. It may replace the uniform lambda-1
candidate only if its projected average is at least `+0.001` higher and the
paired target-bootstrap lower bound against uniform is at least `-0.001`.

This experiment authorizes neither a challenge build nor a submission.
