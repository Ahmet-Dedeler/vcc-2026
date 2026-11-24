# 057 — KOLF under source-only union whitening

## Decision

**Reject KOLF from the union-whitened production model.**

Experiment 055's fixed raw-KOLF blend (`fraction=0.1`) still helps before
whitening, but the gain reverses under every tested nonzero gamma. At the
a-priori production gamma `0.65`, KOLF lowers K562 holdout PDS by `0.006869`:

| gamma | baseline PDS | + KOLF PDS | delta |
| ---: | ---: | ---: | ---: |
| 0.00 | 0.736332 | 0.740202 | **+0.003870** |
| 0.55 | 0.676360 | 0.668305 | **-0.008054** |
| **0.65** | **0.667748** | **0.660879** | **-0.006869** |
| 0.75 | 0.664609 | 0.657758 | **-0.006851** |
| 0.85 | 0.664104 | 0.657148 | **-0.006956** |

The result is not a preprocessing mismatch. On the exact 7,312-gene
experiment-055 shared axis, the unwhitened checkpoint reproduces experiment
055 bit-for-bit:

| | exp055 expected | exp057 observed | absolute error |
| --- | ---: | ---: | ---: |
| baseline | 0.7363493724 | 0.7363493724 | 0 |
| + KOLF | 0.7403417015 | 0.7403417015 | 0 |

At gamma `0.65` on that same shared axis, the KOLF delta is also negative:
`-0.006729`. The rejection therefore survives both the broader
K562-observed axis and the exact exp055 comparison axis.

## Fixed-gamma fold detail

These are five deterministic, disjoint query-row slices. Each query is still
ranked against the complete 240-target gallery.

| query fold | queries | baseline PDS | + KOLF PDS | delta |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 48 | 0.652720 | 0.635722 | -0.016998 |
| 1 | 48 | 0.720014 | 0.722368 | +0.002354 |
| 2 | 48 | 0.644351 | 0.659693 | +0.015342 |
| 3 | 48 | 0.665969 | 0.631276 | -0.034693 |
| 4 | 48 | 0.655683 | 0.655335 | -0.000349 |
| **overall** | **240** | **0.667748** | **0.660879** | **-0.006869** |

Only two of five query folds improve, and the largest loss (`-0.034693`) is
more than twice the largest gain (`+0.015342`). More importantly, the
aggregate is negative at the frozen gamma and remains negative throughout
the entire nonzero-gamma diagnostic sweep. There is no defensible
post-holdout gamma choice that rescues the rule.

## Interpretation

KOLF contributes a modest useful direction in raw space, exactly as
experiment 055 found. That direction does not compose with the covariance
geometry learned from the two predictive sources. Once HCT116/HEK293T union
whitening reweights the feature space, the fixed KOLF addition makes
retrieval worse. This mirrors experiment 054's result for CD4: an auxiliary
donor can help before whitening yet conflict with the dominant production
geometry.

This experiment does **not** show that KOLF is intrinsically useless. It
rejects the specific eligible change tested here: adding raw KOLF at the
already-selected `0.1` fraction before the current union-whitening stage.
Any future KOLF use would need a new source-trained integration rule and a
fresh untouched holdout; tuning the fraction or gamma on this K562 result
would leak the holdout.

## Leakage and protocol audit

- Holdout: K562-GWPS.
- Predictive sources: Orion HCT116 and Orion HEK293T.
- Auxiliary: KOLF2.1J, raw mode, fixed fraction `0.1`.
- Targets: 240 shared VCC-panel perturbations; full 240-way retrieval gallery.
- Operator basis: 12,000 rows, exactly 6,000 sampled from each source donor.
- Operator axis: 17,806 source-union genes on the official VCC axis.
- All 300 VCC panel coordinates were removed before fitting and transforming;
  measured operator/panel overlap is exactly zero.
- K562 rows in the operator: zero.
- KOLF rows in the operator: zero.
- Primary gamma `0.65` was fixed before K562 effects were opened.
- No K562 result selected a mode, fraction, gamma, axis, or candidate.
- The 136 K562-observed non-panel genes outside the source union are absent
  from KOLF as well and remain zero in both predictions.
- No count emission, `cell-eval2`, submission build, or challenge submission.

Hashes:

- target order:
  `60efa65231bcac23892c258e4e757de53ec0b1b3ef6558d7f7f4a678d53da754`
- operator gene axis:
  `280aab34a0e86f27db8ade213e000f9066c8d5f9b89edd3ea966aba061ba9ae9`
- excluded 300-target set:
  `cd3d00f1cd5da660609ad41d85b53e2d5557f7f5edb3e1f64559f78f7fb3cbc4`
- remote whitener:
  `1718c54f230929ac1b8f62f67f721855f4564b35a73f99345730dcd6e537ffa7`

## AWS execution

- Instance: `i-050d8c98b9f09132d` (`r7i.xlarge`, `us-east-1a`).
- Successful SSM command:
  `4c14c879-baea-4d98-a66f-6cb6bbdbe97d`.
- Earlier SSM command that failed closed on S3 authorization:
  `7bb7be4a-7ae8-4ae7-8ba2-9f757cbc2efa`.
- IAM change: the existing inline policy `vcc-scoring-s3` was extended only
  for `exp057-kolf-union-whiten/20260822T210455Z/input/*`,
  `.../output/*`, and the matching list prefix.
- S3 root:
  `s3://<your-s3-bucket>/exp057-kolf-union-whiten/20260822T210455Z/`.
- Runtime: 3 minutes 34 seconds.
- Cloud peak RSS: 3.99 GiB.
- Exit code: 0.

The 1,268,357,498-byte whitener remains only in S3. Compact result artifacts
were downloaded locally:

- `results.json`
- `summary.csv`
- `per_target_pds.csv`
- `whitener_union.json`
- `cloud.log`
- `artifact-manifest.json`
- `exit-code.txt`

Shared AWS instances, volumes, bucket objects, and IAM resources were left in
place as requested.
