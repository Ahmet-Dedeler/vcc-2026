# Experiment 140 — batch-matched controls

## Decision

**Terminal NO-GO.**

Target-total-UMI-weighted batch-matched controls did not improve the frozen
production-donor geometry on the 613-target development split. The selected
final eta was exact identity (`eta = 0`), so the 410-target confirmation set
remains sealed and no production integration, count build, scorer run, or
submission is authorized.

## Primary result

Five-fold nested selection chose:

```text
fold eta: 0, 0, 0.25, 0, 0
final eta: 0
```

Nested out-of-fold macro deltas:

| endpoint | delta |
| --- | ---: |
| raw PDS | `-0.0000879634` |
| signed credit | `-0.0000017071` |

The fixed final candidate is identity, so both fixed deltas are exactly zero.
Only 5 of 30 frozen checks passed. The nested bootstrap intervals were
`[-0.000281220, +0.000073152]` for PDS and
`[-0.000016813, +0.000013334]` for signed credit. Both frozen null p95 values
were zero because null selection also overwhelmingly chose identity.

Biological-family nested PDS deltas:

| family | delta |
| --- | ---: |
| K562 | `+0.0000159933` |
| RPE1 | `-0.0001999168` |
| Jurkat | `-0.0000799667` |

## Frozen component diagnostics

The diagnostic arms reused the primary fold-specific and final eta values
without retuning:

| arm | nested PDS | nested signed credit |
| --- | ---: | ---: |
| all donors | `-0.0000879634` | `-0.0000017071` |
| Orion only | `-0.0000852978` | `-0.0000046390` |
| K562 only | `-0.0000168819` | `+0.0000029312` |

The genuinely new Orion mechanism is negative. The K562 processing-replicate
arm is essentially null and overlaps a15's existing PerturbMatch correction.

## Integrity

- All 613 development targets were supported in all three production donors.
- The 4,992 response coordinates exclude all 300 challenge targets and all
  1,023 gallery targets.
- K562's matched and incumbent control pools reproduced bit-exactly.
- HCT116 and HEK293T legacy controls reproduced bit-exactly while excluding
  154 and 233 NEGR1 cells, respectively, from candidate controls.
- Eta-zero identity passed at raw, centered, pre-whitening, and
  post-whitening stages.
- The `union_pruned` bases retained exactly 21,610 rows for K562-held and
  29,394 rows for non-K562 evaluation.
- The 410-target test was not accessed.

Key SHA-256 values:

- runner:
  `3bea7cb91e7b5d86cb41d2921b156de164e79d9c3dc0837c20f416c4d3fc1a93`
- protocol:
  `2368f5652893d47501bee04d2f81330775b1d5ccfacf44b6e06082c6d1bb6fce`
- source manifest:
  `70c860e308bdf433fbd8e6812e95c5e24b090173f2cc3ac4e6d689ab9195a149`
- preparation manifest:
  `90439ff15cbb18d2e283fdee247ba5cb75d471dd07930a8ce49d9448367df4b2`
- identity manifest:
  `1967a66f463ba78b02c674c09274c9c399fb68b405ff9b0deb045ec4772a5649`
- development result:
  `4160c8d319a0894b1784d0c98436c0d05e8c48ce67044a9810048de1f6f7a3c1`
- development arrays:
  `955036d541d9b8a01120c9859c1d61a8631e7c5be6ec3516b0932599f05472d6`
- development artifact manifest:
  `e40654cfe609e3b24aeb495f8efbd1ddcc0f2ea18e41b3a1fad33b932427f6f6`

Do not open the 410-target confirmation split, retune eta, integrate this
correction into a15, build counts, or submit it.
