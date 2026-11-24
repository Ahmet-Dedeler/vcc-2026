# Experiment 122 Phase 1 Results

## Outcome

`equal_mean_agreement` was selected. It passes every preregistered
leave-one-family-out gate and both complete-selection null gates, so it
**earns a future independent DLD1 confirmation**. DLD1 values were not opened
in this phase.

| Held biological family | Targets | Signed credit | Top-20% sign precision |
| --- | ---: | ---: | ---: |
| Replogle K562 | 267 | 0.022005007 | 0.535544645 |
| Orion | 280 | 0.052960833 | 0.542422605 |
| CD4 | 209 | 0.014930095 | 0.534259555 |
| Macro | — | 0.029965312 | 0.537408935 |
| Worst family | — | 0.014930095 | 0.534259555 |

The selected worst-family signed credit exceeded both complete-selection null
thresholds:

| Null | Draws | p95 worst-family credit | Monte Carlo p |
| --- | ---: | ---: | ---: |
| Target-label derangement | 200 | 0.003142138 | 0.004975 |
| Response-label derangement | 200 | 0.001373382 | 0.004975 |

## Candidate comparison

| Candidate | Eligible | Worst credit | Worst precision | Macro credit | Macro precision |
| --- | --- | ---: | ---: | ---: | ---: |
| equal_mean_magnitude | yes | 0.014930095 | 0.533925567 | 0.029965312 | 0.536831452 |
| **equal_mean_agreement** | **yes** | **0.014930095** | **0.534259555** | **0.029965312** | **0.537408935** |
| equal_median_magnitude | yes | 0.014930095 | 0.533925567 | 0.029965312 | 0.536831452 |
| equal_median_agreement | yes | 0.014930095 | 0.534259555 | 0.029965312 | 0.537408935 |
| sign_consensus_median_min2 | no | 0.009571097 | 0.363796939 | 0.016116238 | 0.409501168 |

Mean and median candidates are identical in LOFO because each prediction has
two training families. The frozen candidate order resolves their exact tie;
agreement-weighted confidence wins over magnitude-only confidence on the
secondary selection statistic.

## Coverage and invariants

- K562 family: 269 valid target rows, 61,913 valid off-diagonal panel pairs.
- Orion family: 295 valid target rows, 88,205 valid panel pairs.
- CD4 family: 209 QC-valid target rows, 59,356 valid panel pairs.
- Final atlas: all 300 target rows supported; 89,541 nonzero panel pairs.
- All target-self residual and confidence entries are exact zero.
- Every non-panel entry in both 300 × 18,533 expanded arrays is exact zero.
- All arrays are finite and the panel expansion is bit-exact.

The frozen source-effect comparison excludes all 300 panel genes from PDS, so
the source-only validation itself is PDS-null.

**Production correction, 2026-08-30:** direct addition to the production
log-effect matrix is not exactly PDS-null. Production exponentiates and
normalizes the complete profile, so changing panel genes changes the
normalization denominator and therefore every non-panel abundance. Experiment
124 measured this on all A/B/C controls. Its first profile-level alternative
was later rejected because it used ordinary log-fold units and failed to
preserve the **summed emitted integer** panel budget after per-gene flooring.
Experiment 128 is the clean replacement: correct `log1p(50000*p)` units plus
exact integer-budget conservation. Its mechanical audit passed, but no
production integration is authorized unless an independent biological
confirmation also passes.

## Source recoverability

No compact source had its panel coordinates zeroed:

- K562 compact: 269 panel target rows and 231 panel response columns.
- Orion HCT116: 290 target rows and all 300 response columns.
- Orion HEK293T: 294 target rows and all 300 response columns.
- CD4 conditions: 283–287 target rows and 285 response columns.

Raw reconstruction is feasible but unnecessary for Phase 1. The local K562
H5AD is present, and the local Orion metadata plus all 109 HCT116 and 223
HEK293T parquet batches are present. The existing builders retain panel
coordinates.

## Verification and hashes

Tests: **11/11 passed** with the standard-library `unittest` runner.

Key SHA-256 hashes:

| Artifact | SHA-256 |
| --- | --- |
| `preregistered_protocol.json` | `446950dbf9ffe5f139df60a353849f993ef66a35a352feb232dab7c576756236` |
| `SOURCE_INVENTORY.json` | `cb407c3ced26c4055f076147826216bb46b80ac8111704913de93d23b614c82b` |
| `SOURCE_ONLY_AUTHORIZATION.json` | `25bc4e404280ecf47dc64ec69e85c8770991e28f953c2e93613b61c4b61e060d` |
| `PREREGISTRATION.sha256` | `e8122e4386e93f60f935df7e04157177366b1443d73a39a343e713b37516fa50` |
| `IMPLEMENTATION_FREEZE.json` | `b3d6046277ef6641f30b18fe39ae866d92088a8c26160f26fe31dd4afaed899f` |
| result `artifact-manifest.json` | `f4b1a42bfffde3c2499130ac378bf2c2cfbf5b07e58dca11a1046d4cce661472` |
| `result-receipt.json` | `bfd137cba001ff43b617ab245b4cd08676652a829e1c44194139ee2ba75d6d24` |
| `selected-panel-residual.npy` | `6ad50b0dbd134b6074da69bf42e94f75aa4cd999abce0d614e7b6cc73f70ec63` |
| `selected-panel-confidence.npy` | `01f2884408d4fde02549abd9b40fe4204961cf711567ab4b32d08c2fc56d6b65` |
| `selected-full-residual.npy` | `518a6e28d0acb6ab2bde9e1f2558e8e6f00bd12f08ac671aefaf110ae1d22972` |
| `selected-full-confidence.npy` | `3e6085a2df8a659a935786f6e284cb85869cd06839924170a10818e27e56f99e` |

The result manifest verified all 18 files with zero hash or size errors. The
result directory is read-only. No challenge truth, held harness/outcome,
A549, Pan H1, DLD1 matrix, `cell-eval2`, count generation, challenge build,
submission, commit, or push was used.
