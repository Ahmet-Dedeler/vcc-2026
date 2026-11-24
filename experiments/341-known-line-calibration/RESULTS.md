# 341 — own-line training labels improve signs but fail retrieval safeguards

September 9, 2026. **DOES_NOT_EARN_FURTHER_DEVELOPMENT.**
Frozen at **feba70a5**. Every sign gate passes; both retrieval
noninferiority gates fail. No retuning, cross-study/count-stage advancement,
production change, package or submission.

This is **known-line, held-target development**, not unseen-line
generalization. A recipient's training-target responses are available;
its test-target responses are excluded throughout fitting. The full
cross-line and same-line models share the same architecture and features.
Repeated training targets across source views are not independent replicates.

## Measured result

Six reused CRISPRi source banks, 1,024 targets, 512 downstream coordinates,
complementary target halves. Responses are unit-RMS raw directions. These
metrics are NOT the official DE scores, MSE/NMAE or A/B/C overall.

| arm | own cosine | retrieval | weighted sign | top-64 weighted sign |
|---|---:|---:|---:|---:|
| source mean | .195532 | .790095 | .589398 | .642517 |
| cross-line diagonal | .198735 | .794614 | .589398 | .642517 |
| cross-line full ridge | .198113 | .759749 | .590195 | .643030 |
| same-line diagonal | .213020 | .789080 | .589398 | .642517 |
| **same-line full ridge** | **.268228** | **.746316** | **.621947** | **.682761** |
| shuffled-label full ridge | .151432 | .599445 | .567841 | .612525 |

Primary minus comparator, 10,000 paired draws:

| comparator / endpoint | difference | target 95% interval | crossed line/target 95% interval |
|---|---:|---|---|
| cross-line diagonal / weighted sign | +.03254973 | [+.02995022,+.03519118] | [+.01553181,+.05421251] |
| cross-line full / weighted sign | +.03175260 | [+.02949554,+.03407009] | [+.01768180,+.05115912] |
| cross-line diagonal / top-64 sign | +.04024402 | [+.03700934,+.04349701] | [+.01811010,+.06847171] |
| cross-line full / top-64 sign | +.03973051 | [+.03672841,+.04281615] | [+.02170877,+.06353616] |
| cross-line diagonal / retrieval | **-.04829784** | [-.05445425,-.04226641] | [-.07344306,-.02750191] |
| cross-line full / retrieval | **-.01343338** | [-.01962721,-.00712844] | [-.03801936,+.00948270] |

The retrieval guard required both lower bounds above -.005. Both
comparators fail; the unresolved crossed contrast against full ridge
does not establish noninferiority.

Weighted sign improves on all six lines and in both fixed cell-number
halves. On Jurkat alone, weighted sign improves +.03380218 versus the
diagonal, target interval [+.03013512,+.03762417], and +.02800431 versus
cross-line full ridge, [+.02478061,+.03144280]. No challenge cell identity
or cross-study transfer is established by that recipient check.

The information regime carries a sign signal in this fixed estimator.
It also increases decoy similarity and reduces target discrimination.
That is a scoped negative for this construction, not evidence that all
use of own-line labels fails. All six competition members remain
unpriced and potentially movable. No candidate is cleared.

## Verification

- Nine pure tests passed before launch.
- All twelve source-mean and cross-line-diagonal prediction replays are
  bit-exact against 334.
- The executed held-target poison checks pass; a fixed linear map of
  the five-source mean equals the mean of the four-source-subset maps
  within 8.88e-15.
- The independent audit reconstructs raw source rows using a different
  RMS reduction, all split/permutation identities, 36 full ridge models
  and 18,432 output optima.
- Maximum normal-equation residual per training row: **3.75e-15**.
  Dense-operator prediction replay difference: **1.60e-14**.
- **36,864** endpoint rows reproduce. Retrieval, both sign proxies and
  absolute/squared errors are exact; cosine error <=2.34e-15.
  Top-64 ordering is independently reconstructed with explicit index
  tie breaking.
- **37** endpoint/comparator interval pairs, including both Jurkat
  contrasts, reproduce via counts-matrix reductions within **3.34e-16**.
- Frozen inputs and model/prediction/result hash receipts verify.

The audit does not repeat upstream donor extraction or re-execute the
original poison test. It verifies fitting and reporting on the same
data, not independent biological replication. Reused targets, correlated
genes, study/line confounding and the selected coordinate axis remain
limitations.

Fit/execution interval **18:35:04.055–18:35:14.994 UTC**, 10.94 seconds
before reporting. Driver peak RSS 1,458,159,616 bytes; run sampled tree
peak 1.36 GiB. Independent audit finished **18:56:22.741 UTC**, sampled
peak 1.05 GiB. Both processes exited successfully.

- Inputs SHA:
  `6a8683cafa510533c28f9ac7b10495ff726aebd384a6740a8d4b54c61db7f990`
- Result SHA:
  `3f4dbdd9eb943d0d200195401a8e19c09d1cce7f31a5d5cf68d2a7cad772531e`
- Independent audit SHA:
  `75288e3da57415bb7ada0d8e44554ca6b0fbb0f44a720fb6b30d4bf267d83bd9`
