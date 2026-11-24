# 326 — Cap dose: local improvement unresolved; interaction distinguishable

Completed September 9, 2026 by Codex. **COMPLETE_DESCRIPTIVE_ONLY.**
No production edits, packaging, submission, reserve access, or external model
call. This is the frozen diagnostic of the already-built 319d3, not a new model.

## Decision

Cap30 minus cap15 gives **+0.0001924 H1-local overall**, paired 95% interval
**[-0.0046193, +0.0031757]**. It does not establish an advantage.

The cap dose response differs from the plain dose response:
`(cap30-cap15) - (plain30-plain15)` is **+0.0061907**
**[+0.0017846, +0.0107100]**. On this reused H1 panel, cap avoids much of the
loss from increasing the plain exponent. This is a descriptive interaction
between complete constructions, not proof of a donor-noise mechanism or an
A/B/C gain.

**Important clamp:** cap15 and cap30 both score zero on H1-scaled MSE at the
point estimate, since both raw values exceed the local baseline .9706906.
Cap30 nevertheless worsens raw MSE by **+.0151952**, with interval entirely
above zero. The near-neutral local overall conceals this cost at its point
estimate; it must not be exported as an A/B/C score forecast.

## All six actual-count metrics

| metric | cap15 | cap30 | cap30 minus cap15 | paired 95% interval |
|---|---:|---:|---:|---|
| PDS (higher) | .8896358543 | .8897759104 | +.0001400560 | [-.0071428571, +.0082633053] |
| MSE (lower) | .9749462419 | .9901414572 | +.0151952153 | [+.0029295985, +.0290382227] |
| NMAE (lower) | .9508364630 | .9515994309 | +.0007629678 | [-.0000029298, +.0015815815] |
| fidelity (higher) | .4975545182 | .4969975816 | -.0005569366 | [-.0014259500, +.0002891363] |
| reach (higher) | .2543472636 | .2566977860 | +.0023505224 | [-.0002869711, +.0070692513] |
| Jaccard (higher) | .0747598330 | .0747968027 | +.0000369697 | [-.0000917164, +.0001621249] |

Local overall: cap15 **.2809531467**, cap30 **.2811455658**.
These use H1's own anchors, not the challenge server's A/B/C anchors.

## Gates and provenance

- Cap15 pre/PM/final effect hashes reproduce 323 exactly.
- Cap15 emitted-count and panel-total digests reproduce 323 exactly.
  **All six raw metric differences are exactly 0.0.**
- A/B/C cap30 source-side hashes match the actual 319d3 build metadata:
  pre `5271106106afe873560e79cc7ecc010ce6771414f6c4d186f559e3b14d86d04d`,
  PM-corrected `b8464100b1aa339b1bf970406e215d55f12df3c0cac05d916d8a0c836c0882f6`.
  This is construction binding, not A/B/C validation.
- Correct H1-target-excluding operator
  `a46fd3239a74752d5c72f5c6ae7169d354793f04f7c22b9b3fd39fe6f74b9741`.
  H1 union OFF; mass fix unset. No held H1 donor loaded.
- Cap30 minus its target derangement: **PDS +.3875350**
  **[+.3032213, +.4712885]**. Negative-control gate passes.
- Five unit tests pass. Aggregate/per-target identities pass; maximum MSE
  identity error 3.78e-15 and variance-budget identity error 2.63e-14.
- 10,000 paired target draws, seed 20260909; fixed models, gallery, and anchors.
  The panel MSE ratio and global claim cap are recomputed on every draw.
  No nonfinite resamples.

The code/protocol/manifest were committed before measurement (`58b566cf`).
Input manifest SHA-256:
`732815b6b0c79f1620b5d0706047987629c31806fadd84cf0257fef243de87e7`.
Machine-readable evidence: `work/results.json`, `work/known-answer.json`,
`work/compose.json`, per-arm metrics/statistics/digests and `work/scores.csv`.

## Resources and limits

Composition took 13.66 s (process peak 10.72 GiB; sampled guard peak 10.64).
Each emitted arm took about 143 s; sampled process-tree peak 4.69 GiB.
The report took 5.87 s. All processes exited zero, serialized under the guard.
The guard's normal command mode does **not** write `--report-file`; resource
evidence is the captured command output/logs, not nonexistent guard JSON files.

This is one repeatedly used 85-target H1 context. The intervals condition on
this context and are not new-context uncertainty or independent confirmation.
323 remains unqualified as a standalone small-arm selector. No dose tuning
followed these results. 319 remains the best recorded A/B/C artifact;
319d3 remains unsubmitted and unarmed.
