# Experiment 132 — source-only result

Status: **terminal NO-GO at the frozen selection gate**

Date: **2026-08-30**

## Bottom line

The exact low-capacity DepMap-feature model family is scientifically distinct
from Experiments 097 and 114, but it is decisively worse than direct
same-target source averaging.

The frozen selector therefore chose `identity_source_mean`. The untouched
410-target source-test split was not opened, no target-block bootstrap was
run, and no held outcome or production action is authorized.

## Distinctness

- Experiment 097 embedded targets with PCA of a donor perturbation-response
  matrix and learned a donor-response-to-response map.
- Experiment 114 used DepMap cosine neighbors to replace 41 weak rows and fit
  no response map.
- Experiment 132 uses authenticated DepMap target features as input to the
  upstream PCA-50/Ridge-10 multi-output regression over pooled source
  contexts.

The input representation, estimator, target scope, and hypothesis are
different from both prior experiments.

## Frozen design

- Upstream commit:
  `a15214780619736d393f40240e56ba992fd416a3`
- Upstream model:
  `PCA(50)` followed by `Ridge(alpha=10, fit_intercept=True)`
- DepMap feature artifact:
  `7f42b585f788db95cdc345d44b44926bd76ce0fc5ae94c544fd2f6929b9ddc51`
- Common non-panel target universe: 1,023
- Selection targets: 613
- Untouched test targets: 410
- Response genes: 4,992
- Removed before fitting and scoring:
  all 300 VCC coordinates and all 1,023 gallery-target coordinates
- Effects: raw and uncentered
- Scored outcomes:
  K562-GWPS, K562-essential, RPE1, and Nadig Jurkat
- Training-only contexts:
  Pan H1, Orion HCT116, and Orion HEK293T
- K562 exclusion:
  both K562 screens were omitted whenever either K562 screen was scored
- Candidates:
  exact source-mean identity and one PCA-50/Ridge-10 model
- Nulls:
  49 frozen target-label derangements

## Selection result

| outcome | identity PDS | Ridge PDS | PDS delta | identity signed credit | Ridge signed credit | credit delta |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| K562-GWPS | 0.851155 | 0.640592 | **-0.210563** | 0.220175 | 0.202860 | **-0.017315** |
| K562-essential | 0.851963 | 0.650388 | **-0.201575** | 0.226771 | 0.190195 | **-0.036575** |
| RPE1 | 0.859765 | 0.644844 | **-0.214921** | 0.232726 | 0.272392 | **+0.039666** |
| Nadig Jurkat | 0.836921 | 0.643783 | **-0.193138** | 0.217759 | 0.197930 | **-0.019829** |
| biological-family macro | — | — | **-0.204710** | — | — | **-0.002369** |

Every required PDS context gate failed. Three of four signed-credit context
gates failed, and macro signed credit was negative.

The real target features did outperform randomized target-feature labels:

| null endpoint | null median | null p95 | observed |
| --- | ---: | ---: | ---: |
| macro PDS delta | -0.307668 | -0.304054 | -0.204710 |
| macro signed-credit delta | -0.026919 | -0.025649 | -0.002369 |

Both add-one null p-values are `0.02`. This establishes that DepMap features
carry response information in this model. It does **not** establish usefulness:
the learned regression discards much of the much stronger same-target
cross-context response signal already available to production.

## Gate result

| gate | pass |
| --- | :---: |
| PDS delta nonnegative in every scored dataset | no |
| Signed-credit delta nonnegative in every scored dataset | no |
| Macro PDS delta at least +0.003 | no |
| Macro signed-credit delta positive | no |
| Macro PDS above target-label-null p95 | yes |
| Macro signed credit above target-label-null p95 | yes |

Selection verdict: `identity_source_mean`.

## Mechanical recovery

The first frozen selection attempt failed before fitting or scoring because
the v1 response axis used Pan gene-name presence instead of Pan's measured
coordinate mask. It produced only an effect-access receipt and did not access
the untouched test split.

`IMPLEMENTATION_RECOVERY_001.json` froze the exact error and authorized only
the mechanical correction: intersect the response axis with Pan coordinates
measured in every row. This changed the response axis from 5,165 to 4,992
genes without changing targets, candidates, hyperparameters, nulls, or gates.
The v1 artifacts were preserved.

## Interpretation

The Ridge model has substantial absolute retrieval signal—raw PDS is
approximately 0.64–0.65, well above chance—but the direct source mean is
approximately 0.84–0.86. DepMap geometry is informative as a generic
target-response prior, yet it is far weaker than observing that target's
response in other contexts.

This closes the exact upstream PCA-50/Ridge-10 replacement family for the
current production setting. Reopening it would require a materially different
mechanism, such as a preregistered residual model that preserves the direct
same-target source response rather than replacing it. No such follow-up is
authorized by this result.

## Verification

- Synthetic tests: 7/7 passed under the memory guard.
- Selection peak RSS: 1.054 GiB internally; 1.08 GiB process-tree peak.
- Source authorization v2 SHA-256:
  `9bca6e620359e44b5440d0b8588412f40f274304ae67d9f55fd7e7d03341e197`
- Selection result SHA-256:
  `db703007b366267949194ca2272ec43c295328aebd5ba2b07b6c74e8413aa1d5`
- Selection artifact manifest SHA-256:
  `822f624f4de6c6e5b14a0fefe16fbe865d343c17fcf14f52fe854acf19e0b7ed`
- Candidate freeze SHA-256:
  `d1b14368f0aab9ee76798b3f298441af8c50cb6d0b5a8f5e8865e70953acc373`
- Untouched source-test access lock: absent
- Untouched source-test result directory: absent
- No `cell-eval2`, counts, challenge build, submission, commit, or push.
