# 103 — residual-only a15-direction confirmation on KOLF2.1J

**Verdict: NO-GO.**

The candidate improved complete-gallery raw PDS on KOLF2.1J by `+0.011716`,
but the target-bootstrap interval crossed zero and one of the 20 frozen
source-label nulls was slightly larger. It therefore fails the preregistered
promotion gate and does not justify a count build or challenge submission.

## Primary result

| quantity | value |
| --- | ---: |
| targets | 273 |
| score genes | 17,806 |
| exact staged-a15 baseline PDS | 0.620408 |
| residual-only candidate PDS | 0.632124 |
| candidate minus baseline | **+0.011716** |
| 95% target-bootstrap interval | **[-0.003138, +0.027365]** |
| conservative 20-null p95 | **+0.011757** |
| add-one one-sided null p | **0.0952** |

The observed gain was only `0.0000404` below the largest frozen null draw, but
the decision rule was fixed before KOLF access: exceeding a linearly
interpolated percentile would not be enough. With 20 nulls, a 5% one-sided
test requires beating all 20.

## Promotion checks

| condition | result |
| --- | --- |
| mean PDS delta at least `+0.003` | PASS |
| target-bootstrap lower 95% above zero | **FAIL** |
| observed delta above conservative null p95 | **FAIL** |
| add-one one-sided null p at most `0.05` | **FAIL** |
| source coverage unchanged | PASS |
| all provenance and reconstruction checks | PASS |
| overall | **NO-GO** |

## What was frozen

- Baseline: exact staged a15 direction (`alpha=0.982`, pruned union whitener,
  PerturbMatch 50:50 K562-slot replacement, CD4 target ratio `f=0.10`).
- The exact Linux a15 effect matrix was re-exported from its original pinned
  build tree. Its array hash is
  `2980d70655f2a9d4a81aa0d8cdb06f2389fca20dddf586d64c132cb05885e5f2`.
- The local Apple-silicon reconstruction matched that matrix elementwise at
  `rtol=1e-12, atol=1e-12`; maximum absolute difference was
  `1.3323e-14`.
- The frozen scored baseline is the preserved Linux matrix, not the local
  reconstruction.
- Candidate: subtract each real donor's projection onto its own non-panel
  donor-mean axis before the unchanged production whitening and CD4 stages.
- PerturbMatch retention is gene-specific: only PM-observed genes halve the
  surviving real-K562 projection; 363 scored K562 coordinates absent from PM
  retain the full real-K562 slot.
- Twenty source-wise target-label permutation nulls, all candidate arrays,
  source axes/loadings, protocol, implementation, experiment-102 evidence,
  and non-KOLF inputs were sealed before KOLF access.

Peak freeze RSS was `20.34 GiB`; scoring peak RSS was `1.28 GiB`.

## Interpretation

The sign and size of the KOLF effect are compatible with experiment 102's
post-hoc three-surface signal (`+0.02050` macro), so donor-axis residualization
remains scientifically interesting. The preregistered evidence is still
insufficient: uncertainty over targets includes zero and the observed gain
does not beat the source-label null.

Do not tune projection strength, axes, target subsets, or null seeds on this
KOLF outcome. A future test needs a different held biological context and a
new frozen protocol. This result concerns a source-donor-axis surrogate, not
official COMPASS. It is effect-level evidence only; KOLF has no biological
control profile, so it cannot authorize count generation or challenge
submission.

Machine-readable results:
`work/results.json`, with per-target values in
`work/per-target-results.npz` and the complete sealed pre-truth bundle under
`work/frozen/`.
