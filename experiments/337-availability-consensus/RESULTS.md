# 337 — MSE improves, but DE harm and unresolved overall: no advancement

September 9, 2026. **DOES_NOT_EARN_FURTHER_DEVELOPMENT.**
Freeze **0bd137ef**. Four actual-count arms completed. No retuning,
production edit, package or submission.

The test normalizes a zero-filled donor consensus by available cap-weight
mass, with two-source minimum support and pre-PM operator-subvector norm
restoration. It adds no new biological information and does not establish
that missing low-expression genes should be treated as random missingness.

## All-six-metric result

85-target reused H1 development panel, current-family cap15 baseline,
H1 union OFF. NMAE has 58 eligible targets under the scorer's unchanged
real-side gate (three empty and 24 below ten reference-DE genes).
These are actual-count local results, not A/B/C forecasts.

| arm | raw PDS | raw MSE | raw NMAE | raw fidelity | raw reach | raw Jaccard | H1-local overall |
|---|---:|---:|---:|---:|---:|---:|---:|
| cap15 | .889636 | .974946 | .950836 | .497555 | .254347 | .074760 | .280953 |
| available | .886415 | .925904 | .953945 | .496500 | .242938 | .074768 | .284175 |
| count-only | .885294 | .922214 | .954349 | .496480 | .246325 | .074873 | .284974 |
| shuffled | .895378 | .977377 | .952592 | .496292 | .249170 | .074685 | .281311 |

Available minus cap15, paired 10,000-target-resample intervals:

| endpoint | difference | 95% interval |
|---|---:|---|
| H1-local overall | **+.00322188** | **[-.00555539,+.00643962]** |
| raw PDS | -.00322129 | [-.01050420,+.00350140] |
| raw MSE (lower better) | **-.04904229** | **[-.07363237,-.01321133]** |
| raw NMAE (lower better) | **+.00310870** | **[+.00134941,+.00495995]** |
| raw fidelity | **-.00105456** | **[-.00199135,-.00013596]** |
| raw reach | **-.01140948** | **[-.02343516,-.00222149]** |
| raw Jaccard | +.00000832 | [-.00030702,+.00032509] |

Overall versus shuffled is +.00286416 [-.00654159,+.00788962].
Versus count-only it is -.00079943 [-.00202448,+.00071504].
Neither required positive overall interval resolves; NMAE worsens and the
count-only point is better. Only the permissive +.005 NMAE upper-bound
gate passes. The fixed construction is rejected.

The positive overall point comes with a substantial MSE improvement:
its local scaled MSE rises from the zero clamp to .0465544. That does
not erase resolved NMAE, fidelity and reach harm.

## Geometry and interpretation limits

Availability affects 844,400 partial-support coordinates across 80
targets. Raw multipliers span 1.1672–2.2896, median 1.5; the factor-4
guard is never active.

Pre-PM operator norms are matched, but delivered pseudobulk norms are
not: available/baseline norm ratio is **.9034 / .9404 / .9670**
(minimum/median/maximum). The emitted non-target Frobenius change is
.2310 relative to baseline; 32,539 nonzero-baseline signs change.
The shuffled median norm ratio is .9963. Thus this is not an
amplitude-isolated missingness effect, nor proof that availability
normalization improves biological direction.

No universal closure of missingness-aware estimators follows. This one
policy does not earn advancement, and this reused H1 panel does not
supply independent-context or A/B/C clearance.

## Verification and provenance

- Eight pure mask/zero/cap/norm/permutation fixtures pass.
- Cap15 pre/PM/final hashes exactly match 323.
- Baseline emitted-count and panel-total digests exactly match.
- All six baseline metric differences are exactly zero.
- Independent sufficient-statistic scoring reproduces every per-target
  row exactly and aggregate metrics within 1.12e-16.
- An independent **dense gene-wise weighted-variance** MSE bootstrap,
  rather than the report's Gram-matrix shortcut, reproduces all three
  overall/raw contrast intervals within 5.33e-13.
- All registered inputs and saved effects/statistics hashes verify.

The audit does not independently regenerate raw count matrices or
reconstruct the source factors. It validates scorer readback and
uncertainty arithmetic, not a second biological replication.

Run **15:46:05.753–15:56:13.856 UTC**, 608.10 seconds.
Driver peak RSS 13,190,365,184 bytes; sampled tree peak 12.27 GiB.
Independent audit sampled peak 3.85 GiB. Both processes exited.

- Input manifest SHA:
  `aeb1ee4625176e2b2b18b17501c565bf15be2977def0b51e0359d1316bf6b37c`
- Result SHA:
  `1aa63c14ffcd68e819b17bbd813da4ce64199d81afead06a5f269361091f0251`
- Independent audit SHA:
  `2f2ee3b3b227be19e9e50bac6d77d269a004c52b886472cfa50a43992657a063`
