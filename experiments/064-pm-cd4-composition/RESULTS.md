# 064 — PerturbMatch × post-whitening CD4 composition results

Date: 2026-08-22

## Verdict

**The two frozen gains compose. Keep the joint PerturbMatch-plus-CD4 direction
as the effect-level candidate for a later count-level confirmation.**

With no retuning, arm D improved both valid held-Orion roles over the incumbent
arm A. The macro gain was `+0.008566434` PDS, with a paired target-cluster
bootstrap 95% interval of `[+0.001532291, +0.015723982]`.

The factorial interaction was small and positive (`+0.000226244`). The result
is therefore mostly additive, with no evidence that one change erases the
other.

## Frozen 2×2 result

| arm | definition | HCT116 PDS | HEK293T PDS | macro PDS | macro delta vs A |
|---|---|---:|---:|---:|---:|
| A | incumbent production union-whitened baseline | 0.814356232 | 0.794261621 | 0.804308926 | — |
| B | A + 50:50 incumbent/PerturbMatch K562 slot | 0.819991773 | 0.797799260 | 0.808895516 | +0.004586590 |
| C | A + post-whitening CD4 mean3 `gene_norm`, fraction 0.05 | 0.819395311 | 0.796729741 | 0.808062526 | +0.003753599 |
| **D** | **both frozen changes** | **0.825092555** | **0.800658165** | **0.812875360** | **+0.008566434** |

### Main effects and interaction

| estimand | HCT116 | HEK293T | macro |
|---|---:|---:|---:|
| PerturbMatch without CD4, B − A | +0.005635541 | +0.003537639 | +0.004586590 |
| CD4 without PerturbMatch, C − A | +0.005039079 | +0.002468120 | +0.003753599 |
| PerturbMatch with CD4, D − C | +0.005697244 | +0.003928425 | +0.004812834 |
| CD4 with PerturbMatch, D − B | +0.005100782 | +0.002858906 | +0.003979844 |
| Combined, D − A | +0.010736323 | +0.006396545 | **+0.008566434** |
| Interaction, D − B − C + A | +0.000061703 | +0.000390786 | **+0.000226244** |

Arm D improved 91/221 targets, worsened 55, and tied 75 on HCT116. It
improved 78, worsened 67, and tied 76 on HEK293T. After averaging each
target's paired delta across both folds, 96 target clusters improved, 67
worsened, and 58 tied.

## Paired bootstrap

- Unit: target gene, paired across A and D and clustered across both held
  folds.
- Replicates: `20,000`.
- Seed: `6402026`.
- Point estimate: `+0.008566434`.
- Percentile 95% interval: `[+0.001532291, +0.015723982]`.
- Bootstrap probability of a positive mean: `0.9924`.

This interval excludes zero in this fixed two-fold effect-level harness. It
does not establish a count-level or challenge-context gain.

The gain is not a one-target artifact. The 5%-trimmed mean target-cluster
delta is `+0.007492005`, and every leave-one-target-out mean remains positive,
from `+0.007469008` to `+0.009442149`. The target-cluster median is exactly
zero because rank PDS has many ties. PerturbMatch-only and CD4-only
target-cluster deltas have low correlation (`0.1085`), consistent with the two
changes contributing mostly complementary corrections.

## Protocol audit

- Valid held roles: Orion HCT116 and Orion HEK293T only.
- K562 was not evaluated as truth because PerturbMatch and the incumbent K562
  bank come from the same Replogle assay.
- Gallery: 221 challenge targets shared by incumbent K562, PerturbMatch, both
  Orion banks, and all three CD4 conditions.
- All 300 panel-gene coordinates were zeroed or excluded from predictions,
  CD4 construction, both frozen operators, and retrieval.
- Operator: the exact experiment-062 source-only union whitener, gamma `0.65`,
  fitted from incumbent K562 plus the non-held Orion donor. PerturbMatch never
  entered covariance fitting.
- HCT116 whitener SHA-256:
  `7cecbf7f36e40d22a6bea1ff3adc86a8788b329def397044fb5b9ba675247cbc`.
- HEK293T whitener SHA-256:
  `97e486e0d1a805e1f24a622f1cefbfb0917e9282c5d0491c4ddb41c8475f2af1`.
- Operator axes: 17,942 genes; retrieval axes: 17,806 genes in both folds.
- Production renormalization was applied independently to A and B after
  whitening. The experiment-062 global preserve scalar was then applied
  independently when adding CD4 to form C and D.
- PerturbMatch was blended 50:50 on 7,075 official non-panel genes shared with
  incumbent K562. Incumbent K562 values were retained outside that shared
  support, which is the union-axis extension of experiment 058's common-axis
  source-slot blend.
- No hyperparameter selection, count emission, `cell-eval2`, challenge
  prediction, or submission occurred.

The absolute PDS values should not be compared directly with experiments 058
or 062 because this joint test necessarily uses their 221-target common
gallery. The within-experiment paired deltas are the estimands.

## Memory, runtime, and infrastructure

- Local direct synthetic harness: 4 passed; peak RSS `0.053 GiB`.
- Final local pytest verification: 4 passed; peak RSS `0.169 GiB`.
- Cloud synthetic tests: 4 passed; peak RSS `0.058 GiB`.
- Worker: Spot `r7i.2xlarge`, instance `i-018a3867545310b8e`.
- Successful heavy wall time: `25.28 s`.
- Internal Python peak RSS: `4.026924 GiB`.
- External process-tree peak RSS: `4.072193 GiB`.
- Guard: `55 GiB`; no kill and no swap.
- Two staging-only attempts stopped before the real runner: first because the
  instance role cannot read the output prefix, then because a provenance path
  named an unstaged launcher. Inputs were restaged under the experiment input
  prefix and the provenance path was corrected before the successful run.
- The instance is verified `terminated`.
- Its delete-on-termination EBS volume `vol-0df50f1185c850671` is verified
  deleted (`InvalidVolume.NotFound`).
- Durable S3 output:
  `s3://<your-s3-bucket>/output/exp064-pm-cd4-composition/`.

## Preservation and hashes

Every downloaded cloud artifact passed
`shasum -a 256 -c cloud-output/manifest.sha256`.

```text
run.py                         39d0f5a29731c5633f3d8fd4035a5a9c3a91cbd435780a201059eaa25a36df81
test_run.py                    bff91bbeb3cb5419e7bdb90ee7e46755038f36c4f48f7dcf5ad7ab094c0ad01b
memory_guard.py                7687db3db528ec7bf724df5645931072f0b355ee7140e75d6a8ffb4241dbb525
aws/run.sh                     21b50545d2727325afb270b6925dfea8dd0739b45315b0f73841506f44b70739
aws/source.tar.gz              5a44d214918d4dfebbace4d9968dbc5181efdee0944d091b94c796408b8335d0
cloud-output/results.json      723c7e94b3aee687b001be435bf7b7a7b55d6a54991404e0ce9dcdd3f59a6113
cloud-output/manifest.sha256   d77b17099a4387818e0764fca267259986b123d92379b1782947a1046562f637
```
