# 332 — absolute error improves; retrieval and recipient gate fail

September 9, 2026. **DOES_NOT_EARN_FURTHER_DEVELOPMENT.**
Main freeze **87adc5db**; depth addendum **f08983a9**, disclosed before any
new accuracy result. All 24 outer-line/target-block splits completed.
No retuning, count realization, package, production change or submission.

## Measured result

Twelve working Broad-16 lines, 82 targets, 512 source-control-selected genes
per fold. Each fit uses 2,048 sampled pairs from eleven training recipients
and the opposite target block. Both held-line responses and the held
target-block labels are absent from fitting; query-target source effects
remain valid inputs. All historical baseline checks and actual-data
held-response poison tests pass.

Losses below are raw response errors divided by source-only RMS (squared
for squared error). **They are not VCC MSE or DE NMAE.** Retrieval is an
82-target, 512-coordinate source-stage statistic, not an A/B/C score.

| arm | absolute error ↓ | squared error ↓ | cosine ↑ | retrieval ↑ |
|---|---:|---:|---:|---:|
| source mean | .758537 | 1.054425 | .108710 | .732673 |
| zero | .729747 | 1.013777 | .000000 | .500000 |
| **TabICL median** | **.727251** | **1.006566** | **.083522** | **.675662** |
| same-fit mean | .729136 | 1.006155 | .086631 | .670117 |
| no recipient | .727255 | 1.005007 | .084944 | .686302 |
| wrong recipient | .728961 | 1.008038 | .076755 | .679088 |
| shuffled training labels | .731361 | 1.018832 | .010321 | .501644 |
| ridge | .729361 | 1.003267 | .102972 | .692738 |
| absolute-error boosting | .732469 | 1.016036 | .062654 | .656818 |

### Paired intervals (10,000 draws)

Positive loss improvement means comparator minus primary. Both intervals
condition on the fitted models; the crossed bootstrap resamples lines and
targets, not studies, cells or a complete refitting procedure.

| primary contrast | point | target 95% interval | crossed line/target 95% interval |
|---|---:|---|---|
| absolute improvement vs source mean | +.03128591 | [.02795385,.03419924] | [.02541088,.03669618] |
| absolute improvement vs zero | +.00249576 | [.00053713,.00474071] | [.00038781,.00491866] |
| absolute improvement vs no recipient | **+.00000416** | **[-.00026110,.00028280]** | **[-.00064268,.00062722]** |
| absolute improvement vs wrong recipient | +.00171010 | [.00144832,.00196004] | [.00112961,.00230346] |
| absolute improvement vs shuffled labels | +.00410996 | [.00224083,.00620061] | [.00188905,.00658120] |
| absolute improvement vs absolute boosting | +.00521750 | [.00480159,.00564985] | [.00439983,.00606367] |
| absolute improvement vs same-fit mean | +.00188458 | [.00147794,.00225819] | [.00106698,.00279288] |
| squared-error improvement vs source mean | +.04785920 | [.03796473,.05687415] | [.03425897,.05920820] |
| cosine change vs source mean | **-.02518719** | **[-.03288359,-.01776639]** | **[-.03735541,-.01356621]** |
| retrieval change vs source mean | **-.05701094** | **[-.07620697,-.03768851]** | **[-.08854073,-.02834143]** |

The median does what its estimand suggests on absolute error, including a
small resolved advantage over zero and the same-fit mean. That does **not**
establish useful target discrimination. Both cosine and retrieval worsen,
and the full model does not resolve an advantage over the recipient-free
version. Beating a wrong-recipient control alone would misread this result.

## Gates and depth

Original failures: absolute-error advantage over no-recipient, retrieval,
cosine non-inferiority, and line robustness. Only **7/12** lines improve
absolute error over no-recipient (required 8). Integrity gates pass.

The pre-outcome depth addendum also fails. Primary absolute-error
improvement over no-recipient is **+.00021738** on the lower-depth six and
**-.00020906** on the higher-depth six. Descriptive correlation with log
total UMIs/KO is Spearman **-.4825**, p=.1121, n=12. This is not a
depth-adjusted biological conclusion; assay batch and line remain confounded.

This closes the registered small-data, fixed-model construction, not all
conditional-median or recipient-response learning. It is reused Cas9
development data, with 82 targets and limited coordinates, not new-context
confirmation for the competition. All six competition members remain
unpriced and could change.

## Verification and execution

- Ten pure unit tests pass.
- 12/12 historical baseline cosine/retrieval arrays reproduce exactly.
- Held-response poison leaves every training/query feature and label
  unchanged; training target blocks and recipient/source indices audited.
- Independent `audit332.py` reconstructs **8,856** per-target arm rows.
  Absolute/squared-error and retrieval differences are exactly zero;
  maximum cosine difference **1.22e-15**.
- Independent paired/crossed bootstrap contractions reproduce the four
  main endpoint intervals within **4.17e-17**. The orthogonal-cosine
  descriptive diagnostic is not independently reproduced.
- Every frozen input and all 24 prediction receipts/hash bindings pass.

Start 13:01:44 UTC, finish **14:01:05 UTC**; runner **3,561.40 seconds**.
Guard sampled peak tree RSS 2.17 GiB. Driver maximum RSS was
2,577,432,576 bytes (~2.40 GiB), and the largest fit process maximum was
1,227,603,968 bytes (~1.14 GiB). The guard sample missed the brief driver
peak and is not a complete maximum-memory measurement.

Hashes:

- Input manifest:
  `67e8ddf18455284ed7b4e0ad06f9d4a72b0f26afaca0e70ae6a3a87953aaa621`
- Main result:
  `598fd615fbcd1d864045d341e395ad16f9dbd3e2794f60ccf5f24057136ac2a0`
- Depth result:
  `d54f1b0da4a5a93a9014bbb1d3b23d465ae20583c7bd5b326c28869fe3f90689`

No spent Broad reserve was read. No Claude worker or built-in subagent,
paid model API, cloud compute, public disclosure or submission.
