# Experiment 123 results — NO-GO

Date: 2026-08-30

## Conclusion

Exact SLIM contains a real unseen-target directional signal, but its
preregistered production blend is not safe enough to promote.

The selected 5–10% SLIM direction blends improved signed credit and top-20%
sign precision in all three outer donor folds. They also beat both
selection-repeated null families. However:

- macro raw-PDS gain was only `+0.0011904762`, below the frozen `+0.003` gate;
- RPE1 PDS was negative at `-0.0011204482`;
- the cluster-bootstrap PDS lower 95% bound was `-0.0012131020`; and
- the top-20 precision bootstrap lower 95% bound was `-0.0002700849`.

Promotion is therefore false. H1 and Orion were not accessed, and no held
evaluator was created.

## What was reproduced

The standalone model is the exact upstream SLIM family from
`RasmussenLab/SLIM` revision `5a7e9ade5d0a`:

```text
Y_hat = G @ W @ P.T + b
```

- train-mean bias;
- pseudobulk response PCA basis;
- unstandardized pinned 64-D STRING embeddings; and
- the upstream two-sided closed-form ridge solve.

The exact upstream model is reported separately from the production
adaptation. The production candidate removes the fit-target response mean,
equalizes donor scale using fit rows only, and row-norm-preservingly blends
the learned direction into the direct donor-transfer incumbent.

Experiment 087 did not run this model. It regressed STRING features onto an
existing model's residual.

## Frozen design

- Source donors: K562-GWPS, RPE1, Nadig Jurkat.
- Target pool: 600 STRING-covered targets shared by all three donors and not
  in the VCC panel.
- Split: 360 fit / 120 tune / 120 untouched test.
- Response axis: 5,997 exact-symbol genes.
- Removed coordinates: all 300 VCC targets and all 600 development targets.
- Outer validation: leave one whole donor out.
- Complete gallery: all 120 test targets.
- Grid: `K={5,10,20}`, ridge `{0.01,0.1,1,10}`, blend
  `{0.05,0.10,0.20,0.35,0.50,1.0}`, plus exact incumbent identity.
- Nulls: 32 target-pairing and 32 Gaussian-embedding draws, each repeating
  complete selection.
- Uncertainty: 20,000 fixed reciprocal-STRING-cluster bootstrap draws.

## Frozen selections

| outer held donor | selected source-only candidate |
| --- | --- |
| K562-GWPS | `K=10`, ridge `10`, blend `0.10` |
| RPE1 | `K=5`, ridge `10`, blend `0.05` |
| Nadig Jurkat | `K=20`, ridge `10`, blend `0.10` |
| deployment selection | `K=20`, ridge `10`, blend `0.05` |

Every selection was sealed before any source test row was read.

## Measured results

| outer held donor | incumbent PDS | candidate PDS | PDS delta | signed-credit delta | top-20 precision delta | exact upstream PDS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| K562-GWPS | 0.86029412 | 0.86211485 | +0.00182073 | +0.00264560 | +0.00207970 | 0.56015406 |
| RPE1 | 0.87857143 | 0.87745098 | **-0.00112045** | +0.00185758 | +0.00052674 | 0.53697479 |
| Nadig Jurkat | 0.86519608 | 0.86806723 | +0.00287115 | +0.00427962 | +0.00235117 | 0.55658263 |
| **macro** | — | — | **+0.00119048** | **+0.00292760** | **+0.00165254** | — |

The exact standalone model is above chance in every donor, establishing that
the upstream family learned some unseen-target biology. It remains far below
the direct donor-transfer incumbent, so only a small blend was competitive.

## Uncertainty and nulls

Cluster-bootstrap 95% intervals:

| metric delta | lower | median | upper |
| --- | ---: | ---: | ---: |
| PDS | -0.00121310 | +0.00118509 | +0.00364374 |
| signed credit | +0.00070649 | +0.00285619 | +0.00507194 |
| top-20 precision | -0.00027008 | +0.00161103 | +0.00342728 |

Null p95:

| null family | PDS | signed credit | top-20 precision |
| --- | ---: | ---: | ---: |
| target pairing | +0.00085668 | +0.00031005 | +0.00027469 |
| Gaussian embedding | +0.00075163 | +0.00029582 | +0.00015543 |

The observed macro values exceed every null p95. This supports a real STRING
signal, but it does not rescue the failed cross-donor and uncertainty gates.

## Failed frozen gates

1. Macro PDS delta was below `+0.003`.
2. RPE1 PDS delta was negative.
3. PDS cluster-bootstrap lower 95% was not positive.
4. Top-20 precision cluster-bootstrap lower 95% was not positive.

No post-outcome hyperparameter, fraction, axis, target subset, null, or gate
may be tuned against this result.

## Verification

- Ruff: passed.
- Unit tests: 13 passed.
- Independent arithmetic reconstruction from `result-arrays.npz`: passed.
- `run.py verify`: passed.
- `cell-eval2`: not run.
- Counts: not generated.
- Challenge build/submission: not run.
- Commit/push: not run.

Key SHA-256:

```text
preregistered_protocol.json
c976bb0eebe37a6c89ab677a432f37ffd7758590a9d8a64627690eeea8210047

SOURCE_ONLY_AUTHORIZATION.json
912b0ac8e4b4a397776f721ff327cc467308d67893e29b8ce7038c9b93ff9918

run.py
30f7d9e48b8906c794dd3ceb00b29622e2cea85778d9365d3e5bfc5357bd7aca

test_run.py
fee41f0f99fb4338d2e5ff3a86221c600df8aa7dcb6b2f8cfbf912cf97fe38a1

plan-manifest.json
4cefff0751184d23f661f7c759755462b5bc0b95dd5b03fceb52dfe38138b91e

selection.json
a22ecbc70646c545863c43efab382e9d6e71ad3ae4364ad8cc8a3790ac760e50

source_factors.npz
64471bbf6b42b33e7b633bbefd65b6a7869481cfb00fb59bf13cfff06d840944

selection freeze-manifest.json
105ffc3c4efaac4687f44c65f45b63013bd3385c58bda8c0542c8e7ff490d21f

source-test-access.lock
b564f5507fa158dc722cd8eab469f058e4bdc041d451283cd8d700c2d1f716e4

results.json
9f60d2c882668a33254e68c6d5609672cf10cc78f43f73ddfdce29ac4363a3d0

result-arrays.npz
d631f17ed1249a7f7677c08ebf741ebc2ab54bf81e51d3af9e6e674fea114bfd

result-manifest.json
12d0c3ab53671f72d3a49ba8a5b2fd94ff19f35fc48d7d056513065b433947b8
```
