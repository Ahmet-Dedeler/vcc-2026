# 062 — CD4 residual addition after source-only whitening

## Question

Experiment 060 showed that a CD4 donor helps before whitening but fails badly
when CD4 itself is passed through a whitener fitted only on somatic donors.
This experiment isolates that mechanism by whitening only the ordinary
two-donor baseline, then adding an unwhitened CD4 mean3 residual in the
original official-gene output axis.

## Frozen protocol

- Main perturbation donors: K562-GWPS, Orion HCT116, Orion HEK293T.
- CD4 donor: equal mean of Rest, Stim8hr, and Stim48hr.
- Target gallery: VCC-panel targets shared by all six banks.
- Development folds: held Orion HCT116 and held Orion HEK293T.
- Candidate grid:
  - modes: `raw`, `gene_norm`;
  - additive fractions: `0.01, 0.02, 0.05, 0.10, 0.20`.
- Selection: maximize the worst development-fold PDS gain, then mean gain,
  then prefer the smaller fraction, then `gene_norm` as a deterministic final
  tie-break.
- Final evaluation: exactly one frozen selected candidate against held K562
  target-role truth, alongside its matching baseline.
- Every fold fits a gamma-0.65 signal-covariance union whitener using only that
  fold's two source donors.
- All 300 VCC target coordinates are excluded from covariance fitting,
  transformation, and retrieval.
- Held truth is centered raw donor effect, never normalized or transformed.
- The baseline uses production median-norm preservation. After residual
  addition, one global scalar restores the baseline's output median norm; this
  cannot alter cosine PDS but keeps output magnitude compatible.

This is an effect-level PDS experiment. It does not emit counts, invoke
`cell-eval2`, create a submission, or call the VCC API.

K562 is a target-role holdout rather than a fully independent donor holdout:
as in the earlier three-line LODO harness, K562 can be one of the two source
donors in the two development folds. Its held-target score is never computed
or used until after the mode and fraction are frozen.

## Safety

- Synthetic local tests run under a 1 GiB process-tree guard.
- Real covariance work runs on one fresh isolated AWS worker.
- The cloud process has both an internal 55 GiB guard and an external
  process-tree 55 GiB guard.
- One experiment process runs at a time.
- Cloud outputs are copied to the experiment's unique S3 prefix and downloaded
  locally before the worker is terminated.

Unique prefixes:

```text
s3://<your-s3-bucket>/input/exp062-cd4-post-whiten/
s3://<your-s3-bucket>/output/exp062-cd4-post-whiten/
```

