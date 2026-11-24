# 339 — small retrieval signal, insufficient against controls

September 9, 2026. **DOES_NOT_EARN_FURTHER_DEVELOPMENT.**
Freeze **64b20a9d**. Twelve outer-line/target-block fits completed.

The zero-anchored empirical marginal map improves retrieval over the
source mean, clipping-only and gene-shuffled curves, but does not beat
the learned diagonal or Gaussian marginal control. Its gain misses the
frozen .005 floor and its cosine loss versus the diagonal exceeds .005.
No retuning, blend, count-stage test or production change is earned.

## Measured results

Six reused CRISPRi lines, 1,024 targets, 512 response coordinates.
Unit-RMS raw source directions, not VCC MSE/NMAE or an A/B/C score.

| arm | absolute error | squared error | own cosine | retrieval | decoy cosine |
|---|---:|---:|---:|---:|---:|
| source mean | .807217 | 1.063200 | .195532 | .790095 | .046059 |
| diagonal | **.770353** | **.971669** | **.198735** | **.794614** | .046296 |
| marginal | .806778 | 1.062885 | .191103 | .793028 | .043090 |
| Gaussian | .806502 | 1.061925 | .194656 | .793186 | .044452 |
| clipping | .806760 | 1.062242 | .193289 | .786969 | .046335 |
| shuffled curves | .808474 | 1.066907 | .189926 | .785352 | .046084 |

Marginal minus comparator retrieval; 10,000 paired resamples:

| comparator | difference | target 95% interval | crossed line/target 95% interval |
|---|---:|---|---|
| source mean | +.00293255 | [+.00193574,+.00390466] | [+.00013215,+.00580927] |
| diagonal | -.00158624 | [-.00282536,-.00035957] | [-.00431455,+.00098835] |
| Gaussian | -.00015799 | [-.00105836,+.00070754] | [-.00170962,+.00138946] |
| clipping | +.00605905 | [+.00531509,+.00682958] | [+.00318240,+.00922550] |
| shuffled curves | +.00767551 | [+.00663291,+.00874183] | [+.00473383,+.01120800] |

Retrieval improves over the source mean on five of six lines. Cosine
falls on all six: macro -.00442932 versus source mean and -.00763216
versus diagonal, the latter crossed interval [-.01008199,-.00574308].
Both lower/higher cell-number halves have positive retrieval contrast;
this is not an adjustment for depth or a study-independent result.

This is evidence of a small correspondence-dependent marginal effect,
not evidence that nonlinear calibration improves on simpler useful
estimators. The primary construction stays closed under its frozen gates.
All six competition metrics are unpriced and potentially movable.

## Verification and provenance

- Ten pure tests pass.
- Five historical raw-loader arrays reproduce exactly.
- All twelve source-mean and diagonal prediction replays are exact.
- Held-line poisoning leaves fit curves and predictions unchanged.
- Independent audit reconstructs every fit curve, diagonal and prediction
  from raw sources using alternate interpolation and accumulation code.
  Maximum prediction discrepancy: **1.89e-14**.
- **36,864** endpoint rows reproduce; retrieval/absolute/squared error
  are exact, cosine differs by at most 2.22e-15.
- Fifteen endpoint/comparator bootstrap pairs reproduce using an
  independent counts-matrix reduction, maximum difference **6.94e-17**.
- All registered inputs and prediction/curve receipts verify.

Run **17:07:26.566–17:07:43.508 UTC**, 16.94 seconds before reporting;
driver max RSS 1,195,851,776 bytes, sampled tree peak **1.14 GiB**.
Audit peak sample .62 GiB. Both processes exited.

The initial freeze's wrong historical output path was caught before
fitting; the unchanged selection and correction to `recovery-work/`
are disclosed in the registration. No outcome-based recovery occurred.

- Input manifest SHA:
  `0f62a420a8cb2c6e7f2b2d31113c72b318a1cb87733bd66c7e78fbc7cd922ccf`
- Result SHA:
  `500d69c8319ae8372f39a930d744961ab726d6ccd821f6224474dbf03461886f`
- Audit SHA:
  `581e54ff1dc5665dc7e8eccb8c61778c5c0a678058af85ec8d6b3a44db4c9c4c`

No H1 response, new reserve, model download, paid service, package or
submission. Same-source reconstruction is not independent replication.
