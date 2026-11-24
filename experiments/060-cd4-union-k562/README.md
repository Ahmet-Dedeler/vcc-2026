# 060 — fixed CD4 rule under union whitening, held K562

## Frozen question

Does experiment 047's fixed primary-CD4 rule help when composed with the
production geometry in the lineage-relevant held-K562 fold?

- source donors: Orion HCT116 + Orion HEK293T;
- held truth: raw centered K562-GWPS effects;
- candidate: CD4 Rest/Stim8hr/Stim48hr mean, centered mean-absolute
  gene-normalized, frozen fraction `0.50`;
- operator: source-union signal covariance, gamma `0.65`, fit only on HCT116
  and HEK293T perturbation rows;
- target gallery: every VCC target shared by HCT116, HEK293T, K562, and all
  three CD4 conditions;
- retrieval axis: every operator-axis gene observed in K562 truth; and
- target exclusion: all 300 VCC target coordinates are removed before
  covariance fitting, transformation, and retrieval.

K562 target names and gene coverage define the frozen evaluation shape. K562
effect values are not touched until the operator is fitted and both prediction
arms have been transformed. No parameter is selected on K562.

The optional raw-CD4 arm is diagnostic only. The predeclared comparison is the
gene-normalized CD4 candidate versus the HCT116+HEK293T baseline after the same
gamma-0.65 operator, with each arm using production-style on-axis median-norm
preservation.

This experiment emits no counts, invokes no `cell-eval2` scorer, and cannot
submit a challenge artifact.

## Infrastructure

Heavy covariance fitting and eigendecomposition must run on AWS. Local work is
limited to synthetic tests under a 1 GiB process-tree guard.

Unique S3 prefixes:

```text
s3://<your-s3-bucket>/input/exp060-cd4-union-k562/
s3://<your-s3-bucket>/output/exp060-cd4-union-k562/
```

The worker must be terminated after `results.json`, `summary.csv`,
`per_target_pds.csv`, `cloud.log`, and the whitener hashes are verified
locally. The unique S3 output prefix should remain until project-wide cloud
cleanup.
