# 105 — TeloHAEC source-transfer gate

**Verdict: NO-GO before residual-model fitting.**

The repaired TeloHAEC extraction was tested as an independent source on a
frozen 300-target non-panel gallery shared by TeloHAEC, K562-GWPS, HCT116, and
HEK293T. The gallery and the 6,937-gene response axis excluded both the complete
VCC-300 panel and all 300 gallery target coordinates before source centering,
normalization, or scoring.

Predictions and 20 target-label null permutations were frozen before either
held HCT116 or HEK293T effect matrix was opened.

| held context | K562 PDS | TeloHAEC PDS | TeloHAEC - K562 | equal blend PDS | equal - K562 |
| --- | ---: | ---: | ---: | ---: | ---: |
| HCT116 | 0.69904 | 0.55799 | **-0.14105** | 0.68079 | **-0.01825** |
| HEK293T | 0.66237 | 0.53583 | **-0.12654** | 0.64506 | **-0.01731** |

The equal-blend 95% target-bootstrap intervals were
`[-0.04155,+0.00407]` on HCT116 and `[-0.04191,+0.00575]` on HEK293T.
Macro equal-blend delta was **-0.01778**. It did not beat the frozen
target-permutation null maximum on either context.

All four promotion gates failed:

- equal blend nonnegative on every surface;
- macro gain at least `+0.003`;
- positive bootstrap lower bound on every surface;
- gain above the null maximum on every surface.

This closes direct TeloHAEC transfer and the preregistered follow-on residual
adapter for now. The source may be reconsidered only with a new mechanism that
explains why a model can recover useful information despite the large,
consistent direct-transfer deficit.

Artifacts:

- `work-source-gate/manifest.json`
- `work-source-gate/frozen.npz`
- `work-source-gate/results.json`

The first freeze attempt stopped before held truth because the extracted source
uses `extraction_report.json` rather than `meta.json`. Its frozen predictions
were preserved under `work-source-gate-failed-20260830T034213Z/`; no held
effect value had been read.
