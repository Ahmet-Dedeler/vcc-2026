# 327 — Coordinate median: no useful advantage over controls

Completed September 9, 2026. **DOES_NOT_EARN_FURTHER_DEVELOPMENT** under the
frozen gates. Codex only; no external model, reserve, production edit, build,
submission, or paid resource.

## Decision

The 25% coordinate-median correction gives **+.0008450 H1-local overall**
with paired 95% interval **[-.0020704, +.0032102]**. Its raw MSE cost is
distinguishable; its intended NMAE improvement is absent. Neither magnitude-only
nor shuffled corrections are beaten. Do not retune this construction to pass.

This is a scoped negative for this estimator on the reused H1 panel, not a
claim that every robust estimator or cross-context response model is impossible.

## Actual-count results

| arm | PDS ↑ | MSE ↓ | NMAE ↓ | fidelity ↑ | reach ↑ | Jaccard ↑ | local overall |
|---|---:|---:|---:|---:|---:|---:|---:|
| cap15 | .889636 | .974946 | .950836 | .497555 | .254347 | .074760 | .280953 |
| median25 | .892297 | .978914 | .951291 | .497699 | .254218 | .074687 | .281798 |
| magnitude25 | .892577 | .978730 | .951120 | .497718 | .254222 | .074712 | .281967 |
| power25 | .898739 | .964103 | .967507 | .497180 | .249772 | .074792 | .279775 |
| shuffled25 | .893277 | .976238 | .950920 | .497599 | .261480 | .074708 | .283676 |

Primary median25 minus cap15:

| metric | raw delta | paired 95% interval |
|---|---:|---|
| PDS | +.0026611 | [-.0021008, +.0072829] |
| MSE | **+.0039680** | **[+.0024183, +.0064442]** |
| NMAE | +.0004544 | [-.0004297, +.0014779] |
| fidelity | +.0001446 | [-.0007923, +.0010685] |
| reach | -.0001298 | [-.0110346, +.0122769] |
| Jaccard | -.0000727 | [-.0001972, +.0000438] |

The local MSE score is clamped to zero for cap15, median25, magnitude25 and
shuffled25, but not power25. Thus a positive local overall point estimate
does not erase the measured raw-MSE cost.

## Controls and mechanism

| median25 minus... | local overall delta [95% interval] | raw NMAE delta [95% interval] |
|---|---|---|
| magnitude25 | -.0001687 [-.0007021, +.0002989] | +.0001705 [-.0001217, +.0005793] |
| power25 | +.0020235 [-.0058203, +.0101644] | -.0162163 [-.0402228, +.0013737] |
| shuffled25 | -.0018776 [-.0057339, +.0010096] | +.0003713 [-.0003558, +.0011606] |

The magnitude-only control has identical pre-PM magnitudes to the candidate
and retains the baseline's pre-PM signs. It is slightly better on every raw
metric at the point estimate, although their differences are unresolved.
The candidate therefore supplies no positive evidence that its altered
pre-PM signs help. The power control buys larger PDS/MSE gains at a substantial
NMAE cost; it is not a candidate selected from this experiment.

The shuffled correction also has a larger overall point estimate. It is one
frozen negative-control realization, **not an estimated permutation-null
distribution**; no permutation p-value is claimed.

## Construction and delivered geometry

597,280 target–gene coordinates across **80/85 targets** have all three
measured source values and belong to the H1 target-excluding operator axis.
All other pre-PM coordinates remain exactly baseline. Missing source genes
are not treated as zero observations.

Each arm preserves the eligible pre-PM subvector norm to 1.78e-15 absolute.
Median25 changes **24,124** nonzero-baseline signs pre-PM. After the full
composition it changes **51,313**, and in the emitted non-target log-pseudobulk
it changes **30,451**. These are counts of changes, not counts of corrections.

The candidate's emitted non-target effect row-norm ratio has min / median /
max **.985860 / 1.004567 / 1.013715**. The ostensibly sign-preserving magnitude
control changes **28,752 emitted signs** and has ratio
**.990566 / 1.004182 / 1.013583**. Dense whitening and nonlinear rendering
explain why the pre-PM identities must not be called delivered identities.

## Verification and limits

- Protocol, code, tests and inputs frozen at **caf418d7** before construction.
- Eight synthetic tests pass.
- Baseline pre/PM/final hashes reproduce 323 exactly.
- Baseline emitted-count and panel-total digests reproduce 323 exactly;
  **all six raw differences are 0.0**.
- H1-excluding operator unchanged; H1 union OFF; mass fix unset.
- All aggregate and MSE-budget identities pass. 10,000 paired target
  bootstrap draws, seed 20260909; no nonfinite draws.
- Composition: 14.63 s, sampled process-tree peak 10.71 GiB. Five emitted
  arms: approximately 143 s each, sampled peak ≤4.69 GiB. Report: 6.13 s,
  sampled peak 3.47 GiB. All serialized and exited zero.
- `work/execution-*.json` records per-stage clocks and process RSS;
  `work/compose.json`, `work/known-answer.json`, `work/results.json`,
  `work/scores.csv` and per-target tables contain the primary evidence.

This is one reused H1 context, fixed gallery/anchors/models, not independent
confirmation or A/B/C score prediction. All four advancement conditions tied
to overall/NMAE and its matched controls fail; only the requirement that
power25 not beat the candidate's overall point estimate passes.

Additional history check during scoring: 263–265's “majority” drops a whole
donor per target according to pairwise row similarity. It is not this
gene-coordinate median, although its frequency-control lesson is relevant.
