# 080 — Leakage-safe listwise metric adapter

Status: **terminal NO-GO; cloud result recovered and verified**.

## Measured result

The frozen non-panel pre-gate selected
`listwise:r=8:beta=0.02:temp=0.03:reg=0.0001` and passed every predeclared
source-only condition. Its untouched non-panel macro PDS gain was `+0.023376`
against a complete-selection permutation-null p95 of `+0.021217`.

That gain did not transfer to the held panel:

| fold | raw-PDS delta |
| --- | ---: |
| Orion HCT116 | **-0.0006884** |
| Orion HEK293T | `+0.0006435` |
| KOLF2.1J confirmation | `+0.0006868` |
| primary macro, HCT116 + HEK293T | **-0.00002245** |

The primary cluster-bootstrap 95% interval was
`[-0.001863, +0.001803]`, with probability positive `0.4975`. The observed
macro delta was below the 100-replicate permutation-null p95
`+0.0001197`. HCT116 also failed the nonnegative-fold and
leave-one-cluster-out gates. The angular and target-coverage safeguards
passed.

`promotion_gate.passed=false`. Do not integrate, retune, or run a production
or count-level confirmation for this adapter.

The complete AWS recovery contains 331 files. All 328 files governed by
`artifact-sha256.txt` verified, as did the pinned source-code hashes, all 303
serialized model files, all nine query files, all three cluster files, the
selection freeze, and the frozen permutation plan. The EC2 worker
`i-08300b395f88f6f55` was stopped and verified in state code 80.

- result SHA-256:
  `9b12cdcb440081bd3910c72d00ba3e41f1a137b2024119c1c1a258c6ddb8ccb2`
- selection-freeze SHA-256:
  `8a7335cf21ff6c576c5ce47eaf3a85d14d6e212a62e3680c8cddcf57204f8524`
- artifact-list SHA-256:
  `971602b00c44dc2bf7747ff34c7ff94131f50c3a947417aa09fc0bc40d228312`

The implementation now satisfies the independent protocol review:

- exact experiment-064 arm-D panel-baseline reconstruction from raw compact
  banks using the original all-mapped-coordinate normalization followed by
  VCC-coordinate blocking, with an independent saved-query equality check at
  maximum absolute error `<=1e-6`;
- explicit separation between that fixed inference baseline and the
  leakage-safe train path: centering, norm calibration, restricted whitening,
  adapter bases/weights, selection, nulls, and clusters never use either
  target panel's coordinates;
- source clustering now consumes a separate panel-row `cluster_query.npy`
  built with the retained-coordinate-safe preprocessing and excludes all VCC
  and historical target coordinates; the restored-coordinate
  `panel_query.npy` is used only for frozen-model outer evaluation;
- the declared objective is now exactly
  `listwise_loss + 0.5 * lambda * ||UV^T||_F^2`; the undeclared factor penalty
  was removed from both serial and batched optimizers;
- candidate selection now exactly maximizes eligible macro tune PDS, with a
  deterministic complexity ordering only for exact floating-point ties; the
  prior near-tie window was removed;
- permanent serial/batch equivalence coverage now spans ranks 4/8/16, both
  teachers, frozen label permutations, and forced gradient clipping;
- separate non-panel and final-panel feature axes, so final PDS excludes the
  300 VCC coordinates but retains historical target coordinates;
- principal-restricted dense whitening for non-panel rows, with an explicit
  algebraic proof and synthetic dense-multiplication regression test;
- complete 600-target galleries for tune and pretest queries;
- complete candidate-grid reselection in every permutation-null replicate;
- exact final real/null weight serialization before panel-truth
  authentication, followed by reload-only evaluation;
- experiment-064 arm-D equality audit on its validated 221-target intersection
  plus explicit fallback-extension audits;
- fold-specific source-query clustering and independent within-fold cluster
  resampling;
- a separate `<10°` median angular guard in HCT116, HEK293T, and KOLF;
- required untouched KOLF confirmation using the pinned three-donor operator;
- source-archive/runner hash pinning, available-memory checks, and broad
  heavy-process collision detection in both launcher and worker.
- the exact 17,806-gene final score axis and the real 49/600 coordinate-less
  historical-target case;
- coherent frozen null plans across all fitting stages;
- nine reusable fit caches plus serial-equivalent same-rank batching;
- exact Python 3.12.12 and NumPy 2.5.2 pins and S3 RUN_ID collision rejection.
- blocked-coordinate-invariant centering and magnitude calibration;
- HCT116/HEK293T-only primary inference, with KOLF used solely as a required
  nonnegative confirmation.

Synthetic and real-shape validation currently contains 24 tests. It covers the
analytic gradient, exact product-only regularization behavior through
serial/batch equivalence, full-gallery PDS, exact macro-PDS selection,
complete-null-selection audit, frozen-model roundtrip, pre-authentication fit
boundary, dense-operator restriction proof, independent exp064
reproduction/fallbacks, fold-specific bootstrap, per-fold angle guard, KOLF
fail-closed behavior, and AWS safety contracts.

The frozen source archive was launched once on AWS as run
`20260830T1228Z-codex-lead`, service
`vcc-exp080-20260830T1228Z-codex-lead.service`, instance
`i-08300b395f88f6f55`. It completed at 2026-08-30 19:03 UTC.

No count generation, `cell-eval2` run, challenge artifact, or submission has
occurred.

## Local readiness evidence

- Guarded test suite: **24/24 passed** under Python `3.12.12` and NumPy
  `2.5.2`.
- Peak process-tree RSS for the full suite: **0.136 GiB** against a hard
  **1 GiB** guard; return code `0`.
- Bounded 36-candidate rank-16 batch benchmark: **0.157 s** for 5 epochs,
  **0.413 GiB** peak RSS. Linear extrapolation remains approximately
  **1.27 hours** for the preregistered 1,212 rank batches and is not an AWS
  runtime guarantee.
- Dependency hash verification, Python compilation, Ruff lint/format checks,
  Bash syntax checks, and JSON parsing all pass. `shellcheck` is unavailable
  on this machine, so Bash received syntax validation only.

Pinned implementation hashes:

```text
run.py                  f472fde750c18f414a2dce9cfb5c8aa1af4d7b732e19130e73eea601cabc304f
stage_cache.py          30dc5778d5bf93b77bdfeb17eaae094192007e25097006caad211353871396a4
test_run.py             f0d85aaec94d37be71075fdf53a2016fc0b266b783fb070be485332329944d6a
artifact_manifest.json  a811e325ce80bda36e3225a454379f1ec3b0872b7037c454265f36d4119380cc
aws/launch.sh           004ca5ed46fa3de986e20425bcacb4c07beb6621fed2acb7549639adb134b344
aws/run.sh              40b287343065554eabe0f6fa5c73652a3df0721b78f823708f67e8e10e6f58de
```
