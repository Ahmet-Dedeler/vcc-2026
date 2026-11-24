# 062 — CD4 post-whiten residual results

Date: 2026-08-22

## Verdict

The mechanism is supported. Whitening only the ordinary baseline and adding a
small unwhitened CD4 residual afterward reverses experiment 060's large loss
and gives a consistent, modest PDS gain on all three LODO target roles.

The development-selected rule was:

- CD4 source: equal mean of Rest, Stim8hr, and Stim48hr;
- mode: centered mean-absolute gene normalization;
- additive fraction: `0.05`;
- selection: maximize worst HCT116/HEK293T fold gain, then mean gain, then
  smaller fraction, then `gene_norm`.

| held target role | matching baseline PDS | selected PDS | delta |
|---|---:|---:|---:|
| Orion HCT116 | 0.803744728 | 0.808562806 | **+0.004818079** |
| Orion HEK293T | 0.787747334 | 0.790164946 | **+0.002417613** |
| K562-GWPS, opened after selection | 0.810157402 | 0.814563972 | **+0.004406570** |

Development mean delta was `+0.003617846`; its worst fold was
`+0.002417613`.

On K562's 242 targets, 90 improved, 43 worsened, and 109 tied at rank
resolution. A 20,000-replicate paired target bootstrap with seed `6202026`
gave a mean-delta 95% interval of `[+0.001285964, +0.007612908]`. This is
encouraging effect-level evidence, not yet a count-level or challenge-context
validation.

For comparison, experiment 060's pre-whitening CD4 blend at fraction 0.50
changed held-K562 PDS by `-0.143085`. The post-whiten result supports the
specific interpretation that source-only whitening amplified CD4 directions
that did not transfer, rather than showing that all CD4 information is harmful.

## Full development grid

| mode | fraction | HCT116 delta | HEK293T delta | worst delta | mean delta |
|---|---:|---:|---:|---:|---:|
| raw | 0.01 | +0.000462947 | +0.000291485 | +0.000291485 | +0.000377216 |
| raw | 0.02 | +0.000874456 | +0.000240047 | +0.000240047 | +0.000557251 |
| raw | 0.05 | +0.002383320 | +0.000617263 | +0.000617263 | +0.001500291 |
| raw | 0.10 | +0.002880560 | -0.000325778 | -0.000325778 | +0.001277391 |
| raw | 0.20 | +0.005641096 | -0.003583553 | -0.003583553 | +0.001028771 |
| gene_norm | 0.01 | +0.000943040 | +0.000668701 | +0.000668701 | +0.000805871 |
| gene_norm | 0.02 | +0.002211858 | +0.000960187 | +0.000960187 | +0.001586022 |
| **gene_norm** | **0.05** | **+0.004818079** | **+0.002417613** | **+0.002417613** | **+0.003617846** |
| gene_norm | 0.10 | +0.006944206 | +0.001988958 | +0.001988958 | +0.004466582 |
| gene_norm | 0.20 | +0.008607387 | -0.002417613 | -0.002417613 | +0.003094887 |

The conservative worst-fold criterion selected `0.05`; a mean-only rule would
have selected `0.10`.

## Protocol audit

- Gallery: 242 VCC targets shared by all six banks.
- Gamma: fixed `0.65`.
- Every fold fitted a new union covariance on only its two source donors.
- Operator axes: 17,942 genes in both development folds; 17,806 for held
  K562.
- Retrieval axes: 17,806 genes for HCT116 and HEK293T; 7,312 for K562.
- All 300 VCC panel coordinates were zeroed/excluded before covariance fitting,
  transformation, CD4 construction, and retrieval.
- Every held truth matrix remained raw-centered and untransformed.
- The K562 candidate grid was never evaluated. Exactly the one frozen
  `gene_norm, 0.05` candidate was compared with its matching baseline.
- K562 is a target-role holdout, not a donor-independent outer holdout, because
  it served as a source donor in the two development folds.
- No count emission, `cell-eval2`, VCC preparation, or submission occurred.

## Memory and runtime

- Local synthetic tests: 3 passed; peak process-tree RSS `0.069 GiB`.
- Cloud worker: Spot `r7i.2xlarge`, instance
  `i-023c509eada647e81`.
- Launched `2026-08-22 21:50:23 UTC`; successful service started
  `21:52:58 UTC`; outputs reached S3 at `22:09:14 UTC`.
- Heavy wall time: `15:39.61`.
- Internal Python peak RSS: `19.491940 GiB`.
- External process-tree peak RSS: `19.536476 GiB`.
- Guard limit: `55 GiB`; no kill, no swap, exit code `0`.
- A first cloud-side test invocation used the Amazon Linux Python 3.9 default
  and stopped before heavy work; the wrapper was pinned to Python 3.12 and the
  successful run followed on the same isolated worker.
- Termination was requested at `2026-08-22 22:12:19 UTC`; the instance is
  verified `terminated`, and its 100 GiB delete-on-termination EBS volume
  `vol-0ae17b6c069bede41` is verified deleted.

## Preservation and hashes

Durable S3 prefix:

```text
s3://<your-s3-bucket>/output/exp062-cd4-post-whiten/
```

The complete 3.844 GB output was downloaded to `cloud-output/`. Every local
file passed `shasum -a 256 -c cloud-output/manifest.sha256`.

Key SHA-256 values:

```text
run.py                         16417a9b985d9aad8993bf2f24ce550ea1e6871d9401bb4672d11ec6d0046176
test_run.py                    3acaf3d7c58392f32b27902bf805e13ab3c49cb86fb69eea4a8e0eb4b32f22b2
memory_guard.py                c907bb8ebc8fbd48af4b412bd5d543a04c5b7e9f0975c885826b90b119180776
aws/run.sh                     330156e73244485973c0256dc0118b90525d248d11cc706fba37bef00e73fa0d
aws/source.tar.gz              8c8d01180b3216e09ede83d17b13d0782981ec189015689d799732e175cfdf1d
cloud-output/results.json      3b2195c997d1a14cfd95623d85c83b34e85f29f41f074bb05741df6a7affe5c6
cloud-output/manifest.sha256   f9ff3e15f74808919ecadf8db0d6251966908c60d58f31e8e9a54ece3f5d9e27
held HCT116 whitener           7cecbf7f36e40d22a6bea1ff3adc86a8788b329def397044fb5b9ba675247cbc
held HEK293T whitener          97e486e0d1a805e1f24a622f1cefbfb0917e9282c5d0491c4ddb41c8475f2af1
held K562 whitener             13d988fcc3b129d25d254524db73e8f8253143fa493ee64d9654f0f61bc3ac55
```

All input hashes, including every donor matrix, are embedded in
`cloud-output/results.json`.
