# Experiment 129 results

Status: **source-only NO-GO at selection; untouched source test not opened**

Date: **2026-08-30**

## Scientific boundary

GSE295214 is Pan2025 H1, already used by Experiments 117 and 118. Those
experiments close direct Pan/H1 residual integration, fractions, target masks,
and Pan-derived confidence rankings.

Experiment 129 therefore evaluated a materially distinct family: Pan effects
were used only as a cross-target neighborhood geometry. The query target's Pan
effect was never added to its prediction. Pan-nearest training targets
retrieved cross-context tangent residuals learned from RPE1 and Nadig Jurkat,
with strict opposite-context teacher isolation.

No held HCT116, HEK293T, KOLF, Pan panel, DLD1, HCASMC, H1-2025, or challenge
outcome was opened. Orion HCT116 and HEK293T were used only as source
predictors and whitening inputs.

## Source reproducibility

The original preregistration and source inventory remained unchanged:

```text
protocol SHA-256:
05af0bb9da6321cdff81818209e0befa1f74552e96964562d37caef378bac998

source inventory SHA-256:
23cb28ba364a94ca4b73b728508d44e825454e43bcfe3ddb2fbc91b834f3f6a4

plan manifest SHA-256:
113ced458a702a8825f659505c212e548db77b2a91707f26d46df3efc4bea570

original authorization SHA-256:
193a613903c754710e5b612e3cf007b64ae99f21c113e2700f908df72013fed4
```

The first guarded selection attempt failed before metrics or selection because
the loader incorrectly required every Pan row to measure every response
coordinate. Recovery receipt 002 was written before changing the
implementation. On the already-authorized 100 selection rows:

- 200 of 573,000 row-coordinate entries were unmeasured;
- every row had exactly two unmeasured coordinates;
- the same two response genes were unmeasured in all 100 rows;
- every unsupported compact-donor delta was bit-exact zero.

The recovery kept the frozen 5,730-gene response axis and explicitly applied
the row-specific measured mask, requiring unsupported deltas to be exact zero.
No targets, splits, candidates, angle, folds, nulls, bootstrap draws, metrics,
or gates changed.

```text
recovery 002 SHA-256:
25d2f4458496f5c19554fb14162dc09336ac644c2916496677821a32da177308

recovery authorization v2 SHA-256:
e9380f14f186e1221443923a126ef2701bf8eb26f03503ff1d3632a513556f6b

repaired implementation SHA-256:
2018f7db99117c7020ca2ea4b979d10daa8b4c19c773ddc8acae628df62da19e

synthetic test SHA-256:
0a4923e61d3386d7ebc09f0ee1027869aaa4cbda3ea255cc509cf06cd6be0d2e
```

Ruff passed and all 13 synthetic tests passed. The source-only selection ran
under the repository memory guard. The external guard sampled a 2.96 GiB peak;
the process's internal peak-RSS record was 3.616 GiB.

## Source-only selection

Five whole-target folds produced out-of-fold predictions for all 100 frozen
selection targets. Identity was compared with the three fixed five-degree Pan
neighborhood rotations.

| Candidate | RPE1 PDS delta | Jurkat PDS delta | Macro PDS delta | RPE1 signed-credit delta | Jurkat signed-credit delta | Eligible |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| identity | 0 | 0 | 0 | 0 | 0 | yes |
| Pan k=5 | -0.00848485 | -0.01353535 | **-0.01101010** | +0.00956154 | +0.00557905 | no |
| Pan k=15 | -0.01656566 | -0.01686869 | **-0.01671717** | +0.01229346 | +0.00668915 | no |
| Pan k=30 | -0.02353535 | -0.01818182 | **-0.02085859** | +0.01360005 | +0.00748545 | no |

All nonidentity candidates improved directional signed credit in both source
contexts, but all reduced complete-gallery raw uncentered PDS in both
contexts. The PDS loss increased with neighborhood size. Because eligibility
required both endpoints to be nonnegative in each context, the frozen selector
chose exact identity.

The immutable selection artifacts are:

```text
selection freeze manifest SHA-256:
7e011515ee51e16873dedc419b3fa9119a5371fc9f0a418fa3d8d0cc466a9d40

selection JSON SHA-256:
88360910c9ad62b237aa539eeb3a2383ac5e3965fe7525fdf75750d217e796e8

geometry SHA-256:
f1f40ac49cc37c7308e986e40217f50997b4262262901794e5082d6e0e2da2d6

selection cache SHA-256:
39068daa6dcb6e9359b548f54a0c17981b7ea1d6fd82e2602a6213c7068dd269
```

## Untouched source test

The 55-target source test was **not opened**. No source-test access lock or
source-test result directory exists.

This is required by the frozen protocol: identity at source-only selection is
a terminal NO-GO. Consequently, the frozen test bootstrap and target-label
permutation arrays remain unused rather than being consumed to characterize a
candidate that already failed selection.

## Decision

**NO-GO.** Broad Pan response neighborhoods contain a reproducible
cross-context directional signal, but the fixed low-capacity learner converts
that signal into a substantial PDS loss. It is not a deployable improvement
over the frozen a14/arm-D-style source baseline.

Do not retune neighborhood size, angle, target masks, fractions, or confidence
rankings against these now-consumed 100-target outcomes. A future attempt
would require a materially different mechanism and newly preregistered,
unconsumed source-only selection and test splits.

No held evaluation, `cell-eval2`, count generation, challenge build,
submission, commit, or push was performed. Even a source-only pass would have
authorized only a separately reviewed independent confirmation, not
production; this experiment did not reach that stage.
