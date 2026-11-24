# 323 — Current-family H1 validation: engineering PASS, small-effect selection unresolved

Completed September 9, 2026 by Codex alone. No subagents, external model
calls, cloud rental, new primary-data acquisition, reserve access, production
change or submission. Source code and registrations were committed before
measurement; three implementation repairs are recorded separately.

## Decision

**NOT QUALIFIED AS A STANDALONE SELECTOR for this arm family.**
All four known board comparison signs are correct on the repaired H1 panel,
but **0/4 paired 95% intervals exclude zero** in the expected direction.
The large target-shuffle negative control passes decisively.

This is useful regression machinery, not a new winning model or an independent
confirmation. One reused H1 hESC context cannot establish cancer-line transfer.
The H1 union is necessarily OFF throughout.

Primary machine-readable result: `work/results.json`. Raw/scaled arm table:
`work/scores.csv`. Every metric is computed from actual emitted count blocks,
not from design profiles alone.

## 1. Two validation problems repaired

### Historical H1 whitening used the wrong exclusion set

170/273 used the correctly labelled **A/B/C** whitening coordinates on H1.
That is not an arbitrary gene-axis permutation, but it leaves **80 of the
85 H1 target genes inside the dense transform**. The proper H1 axis has
18,157 columns, versus 17,942 for A/B/C.

With historical unset H1 inputs, the additive contribution of these excluded
coordinates to other, scored coordinates has Frobenius norm **0.077267 times**
the total transformed output norm; median per-row ratio **0.062875**, maximum
**0.214592**. These are norm ratios **before preserve normalization and
emission**, not variance explained, score gains, or proof of held-H1-response
leakage. No held response was used in that measurement.

The replacement operator was fitted from the same three non-H1 donors,
essential-row pruning, 12,000-row donor cap, seed 0, and official gene axis
minus the 85 H1 targets:

- Path: `work/whitener_h1_pruned.npz`
- SHA-256: `a46fd3239a74752d5c72f5c6ae7169d354793f04f7c22b9b3fd39fe6f74b9741`
- Axis SHA-256: `2a9c18730b9f77dc47990e77223d6d41de2da616d92fbdaa7cefdc9ba14d88b2`
- Basis: 29,394 rows × 18,157 columns.

Production's loader rejects the old operator on H1 and accepts the replacement.
A target-coordinate-only input produces **exactly zero** on scored coordinates
through the replacement operator.

Evidence: `inputs.json`, `work/compose.json`,
`work/target-coordinate-mixing.json`, `work/whitener_h1_pruned.json`.

### A floating-point normalization order changed a DE call

The first fast evaluator reproduced emitted totals, PDS, MSE and fidelity,
but failed the other three metric tolerances. For MLLT1, FGFR1's adjusted
p-value was **0.04999425835** through Scanpy versus **0.05000068924**
through the initial helper. Multiplication by a reciprocal versus division
by the scaled library total changed rank ties.

The evaluator now follows Scanpy's arithmetic order. The production emitter
still uses its original control block: the emitted count digest is unchanged
between the failed and repaired known-answer passes.

On the full MLLT1 reference comparison, corrected p-values agree with actual
Scanpy to 3.33e-16 and adjusted p-values to 8.88e-16, with zero changed calls.
No tolerance was relaxed.

Evidence: `EVALUATION-AMENDMENT-3.md`, `work/de-parity-diagnosis.json`,
`work/metrics-known-r2-failed.json`, `work/metrics-known.json`.

## 2. Construction and scorer gates

| gate | observed |
|---|---|
| Production donor preparation → historical H1 pre-whitening matrix | max absolute difference **0.0** |
| Cap15 A/B/C pre-whitening and PM-corrected hashes → actual 319 build metadata | **both exact** |
| Re-extracted CD4 complete-mean3 coverage | **77/85**, matching 273 |
| H1 donor union loaded | **none** |
| Historical rc_callset emitted panel-total digest | **exact** `453cdc8d…` |
| Historical emitted count digest after DE evaluator repair | **unchanged** `16b0f69d…` |
| Full-scored historical PDS/MSE/fidelity/reach/Jaccard aggregates | max absolute difference **6.67e-16** |
| Historical NMAE aggregate | absolute difference **3.86e-8**, below frozen 1e-6 |
| Local-anchor constructor vs official full pipelines, synthetic fixture | baseline six metrics exact; replicate max difference **1.11e-16** |
| MSE bootstrap vs each actual full-panel score | max absolute difference **2.67e-15** |
| Bootstrap variance budget vs official `_across_pert_budget` | max absolute difference **5.27e-14** |

Ten unit/regression tests pass, including actual Scanpy DE on a synthetic
fixture and the bootstrap budget on explicitly duplicated target rows.
The synthetic anchor check is additional to those ten tests.

Current-family settings: alpha .982; CD4 f .34 and k=5; re-extracted H1-target
CD4; gamma .65/lambda 0; rc_callset; **mass fix unset**; H1 union OFF.
The historical known-answer adapter alone uses its historical minnorm setting.
No H1 truth enters donor fitting, operator fitting, residual sizing or emission.

## 3. Actual six-member results

Raw metrics on **H1**, not A/B/C:

| arm | PDS ↑ | MSE ↓ | NMAE ↓ | fidelity ↑ | reach ↑ | Jaccard ↑ |
|---|---:|---:|---:|---:|---:|---:|
| uniform | .889356 | .984517 | .950676 | .496694 | .247506 | .074627 |
| plain15 | .885574 | .960565 | .951287 | .496952 | .251786 | .074860 |
| plain30 | .880952 | .989701 | .952546 | .496002 | .241778 | .074854 |
| cap15 | .889636 | .974946 | .950836 | .497555 | .254347 | .074760 |
| parity15 | .886134 | .971604 | .951255 | .496656 | .250302 | .074880 |
| shuffled cap15 | .504062 | 1.111045 | 1.002851 | .483573 | .055946 | .072639 |

The reference has 85 × 400 perturbed cells and 18,400 controls. Its own-target-
excluded DE-count distribution has median **89**, quartiles **7 / 1,019**,
range **0–5,484**. NMAE has **58** reference-eligible targets; reach has **82**.
These are measured cohorts, not a uniform DE-count replacement.

H1-local anchors were rebuilt with the official dispersed baseline constructor,
five official seeded disjoint half splits, separate controls per half, and
NMAE's full-reference gate. They are not server-attested A/B/C anchors.
See `work/h1-anchors.json` and `ANCHOR-NOTE.md`.

## 4. Retrodiction and uncertainty

Local overall means the six **H1-anchor-scaled** members, with the specified
MSE clamp and NMAE floor. These deltas are not predictions of board deltas.

| contrast | known board sign | local overall delta | paired 95% interval |
|---|---|---:|---:|
| plain15 − uniform | positive | +.001136 | [−.003737, +.006631] |
| cap15 − plain15 | positive | +.000483 | [−.001928, +.003895] |
| plain30 − plain15 | negative | −.005998 | [−.013787, +.000234] |
| parity15 − plain15 | negative | −.001876 | [−.004622, +.001477] |

**4/4 point signs agree, 0/4 resolved.** PDS by itself would even give the
wrong sign for plain15 versus uniform; all six members matter.

Cap15 minus its frozen target derangement has raw PDS **+.385574
[+.303081, +.468211]**. The negative control passes.

The 10,000 paired resamples use targets, seed 20260909, and condition on the
fixed reference gallery, H1 anchors and constructed models. MSE remains a
panel ratio; its global noise-credit cap is recomputed, not held at one.
Indeed the cap coefficient falls as low as .898 in some resamples despite
being one in every full-panel arm. No nonfinite resamples occurred.

## 5. Eligibility and what comes next

`work/eligibility.json` and `work/coverage.csv` consolidate 267/273/278/285's
already-used metadata. Coverage rows count targets present in raw metadata,
not just targets retained by a donor bank's minimum-cell filter. For example,
K562 has 272 raw panel targets versus 269 retained donor rows.

H1 is useful for engineering regression, large-effect rejection and descriptive
six-member comparisons. It is **not qualified to pick small norm-weighting
changes alone**, and cannot test the H1 residual on its own recipient.
The Orion/K562 references have different depth/DE problems; Nadig and Broad-16
have zero overlap with the current 300-target panel. No perfect replacement
reference has been established here.

Next bounded research: a source-only feasibility design for independently
repeatable downstream-response residual labels. Inventory actual
study/line/construct units first, compare the proposed label and estimator with
132/161/224/268/269/283/285, and reject a renamed scalar-reliability or basal-ridge
rerun. Preserve every sealed/reserved target and file. A new model is **not**
claimed launched by this next-step note.

## 6. Operations and reproducibility

- Operator fit: **352 s**, process peak **12.75 GiB**.
- Six actual-count arm evaluations: **936 s**, process peak **5.04 GiB**.
- H1 anchor construction: **122 s** computation / **131 s** guarded command,
  process peak **17.71 GiB**.
- Paired report: **4.6 s**, peak **0.43 GiB**.
- Retained experiment scratch: **4.6 GiB**, within the 25-GiB bound.
- Every heavy stage ran serially under the memory gate. No research worker
  remains running; the pre-existing board snapshot process was preserved.
- Authenticated status at **08:09:50 UTC September 9**: 319 remains published,
  **0.21865973944346428**, rank **5**. Rank was 6 at the earlier takeover check;
  no submission was made to change it.
- The old paused `VCC research watch` was updated, not duplicated: hourly
  continuation in this task, no subagents, no automatic rollback/reset upload,
  no paid resources/reserve access without separate authority.

Reproduction commands are the staged modes of `run323.py`, `evaluate323.py`,
`anchors323.py`, `diagnostics323.py`, and `report323.py`. Frozen input manifests
and named repair commits remain part of the record. Binary statistics and the
operator stay local; compact results, provenance, tests and logs are retained
in Git. Production code and the best package are unchanged.
