# Experiment 142 — batch-disjoint donor reliability weighting

Status: **corrected donor-wide plan passed; Polars streaming freeze pending**

Date: **2026-08-30**
Runtime repair: **2026-08-31**

## Question

Does within-donor repeatability across disjoint technical/GEM batches identify
which production donor should receive more weight for a particular target?

This experiment does not alter any donor response coordinate. It estimates one
target-level reliability value per production donor, converts that value to a
within-donor midrank, and uses it only to reweight the existing a15 donor slots.

## Frozen source-only construction

The production donors are:

- `k562_gwps`;
- `orion_hct116`;
- `orion_hek293t`.

The target and response axes are reused byte-for-byte from Experiment 132:

- 613 development targets from `selection_targets.txt`;
- 4,992 safe response genes from `response_genes.txt`.

No Experiment 132 410-target effect value is read by any pre-outcome command.

The first metadata plan used a different batch partition for each target. A
pre-freeze audit rejected it because complete-gallery PDS requires every row
in one split half to share the same donor-wide batch partition. That plan is
preserved as superseded and was never used to read source expression.

The corrected plan uses one partition per donor and split:

1. aggregate total UMI over all 613 targets for every batch;
2. construct 512 deterministic donor-wide candidate partitions for each of
   eight frozen split seeds;
3. choose the candidate that first maximizes supported targets, then the
   worst, 1st-percentile, 5th-percentile, and median target-level UMI balance;
4. apply that one batch partition to every target in the gallery;
5. require at least five target cells and four occupied batches in each half;
6. require at least six valid split pairs per target and eight distinct
   donor-wide partitions.

The corrected v2 plan passed on 2026-08-30:

| donor | all 613 targets valid in all 8 | minimum cells/half | minimum batches/half | minimum UMI balance | median UMI balance |
|---|---:|---:|---:|---:|---:|
| K562-GWPS | yes | 6 | 6 | 0.5278 | 0.8968 |
| Orion HCT116 | yes | 5 | 4 | 0.4245 | 0.8690 |
| Orion HEK293T | yes | 6 | 5 | 0.3983 | 0.8587 |

The first K562 source build then completed as a reconstruction pilot:
all 613 targets had eight valid pairs and pooled raw log-pseudobulk reproduced
the compact donor within `1.20e-7` maximum absolute error. Before any source
evaluation or Orion expression scan, an independent numerical audit found
three downstream issues: bootstrap resamples were not reranked, target nulls
were downcast before a strict comparison, and later-stage directories were
not atomic. `PREFREEZE_AMENDMENT_002.json` supersedes that execution chain.
The v1 K562 artifact is preserved as provenance but is not accepted into the
repaired v2 chain.

A second independent control-plane audit then found that the corrected clock
timestamp was bound in the authorization but not enforced by `run.py`, and
that only source builds—not source evaluations—held an exclusivity lock.
`PREFREEZE_AMENDMENT_003.json` repairs both boundaries. The next accepted
chain is authorized only by `SOURCE_ONLY_FREEZE_AUTHORIZATION.v3.json`;
v2 remains immutable superseded evidence.

The v3-authorized K562 build and evaluation then passed, but two HCT116 build
attempts failed closed before publishing an HCT artifact. Whole-file Parquet
materialization exceeded the frozen memory envelope: first the 17.5 GiB
internal limit, then the 18 GiB process-tree guard under an alternate Arrow
allocator. `PREFREEZE_AMENDMENT_004.json` replaces that execution chain with
256-row filtered Parquet scanning. The arithmetic and filters are unchanged;
only peak materialization is reduced. The accepted replacement chain uses
`SOURCE_ONLY_FREEZE_AUTHORIZATION.v4.json` and `work/*-v3`.

A final read-only audit caught that later artifact stages still permitted a
fallback to the unreliable host clock. `PREFREEZE_AMENDMENT_005.json` removes
that fallback: every plan/build/evaluation/gate artifact now requires an
explicit timestamp resolving to August 31 UTC. The accepted authorization is
`SOURCE_ONLY_FREEZE_AUTHORIZATION.v5.json`; v4 was superseded before freeze.

The v5 K562 chain then reproduced exactly, but the first real HCT116 run showed
that Arrow still decoded the one-row-group nested Parquet columns too eagerly:
process-tree RSS reached 22.18 GiB and the 18 GiB guard killed it before any
HCT artifact was published. A Polars streaming probe on the largest HCT116
file decoded the same nested columns at 3.83 GiB peak. The final bounded-memory
reader therefore uses Polars lazy streaming with 256-row output batches.
On a real HCT116 file, its response sums, cell counts, and total UMI were
bit-identical to the former whole-file reader.
`PREFREEZE_AMENDMENT_006.json` records the replacement; the accepted chain is
`SOURCE_ONLY_FREEZE_AUTHORIZATION.v6.json` with fresh `work/*-v4` directories.

The source build streams one donor at a time. Orion files are additionally
streamed through Polars lazy Parquet scans in filtered 256-row record batches.
For each half it forms

```text
b[s,h,t] = log1p(50000 * sum_counts[s,h,t] / total_umi[s,h,t])
c[s,h,t] = b[s,h,t] - mean_t(b[s,h,t])
```

where the mean is over the frozen source targets with nonzero support in that
half. This centering makes the repeatability statistic invariant to a shared
control-profile translation. Split-derived rows are never used to fit or
refit a whitener. Orion's 22 duplicated gene symbols follow Experiment 140's
last-occurrence coordinate rule. Before a donor build is accepted, its pooled
raw log-pseudobulk must reproduce the frozen compact donor's `delta + b_ctrl`
on the same coordinates within `2e-5` maximum absolute error.

## Fixed geometries and reliability

Experiment 140's frozen `gamma=0.65` operators are reused:

- `k562_gwps` is evaluated under `whitener-non_k562`;
- both Orion donors are evaluated under `whitener-non_k562` and
  `whitener-k562`.

For each split, both A-to-B and B-to-A complete-gallery retrieval PDS values
are computed with exact midranks. Donor-target reliability is the lower
quartile of the 16 directional values. Orion must pass under both applicable
operators.

Each target must also exceed its own 95th percentile under the frozen
depth/batch-count/norm-matched target-pairing null. Null-level targets receive
`z=0`; survivors are converted to an exact midrank `z in [-1, 1]`. The only
candidate weights are

```text
w[d,t](eta) = 1 + eta * z[d,t]
eta in {0, 0.25, 0.50}
```

and are normalized across the donors available in the relevant outer
context. There is no coordinate gate, sign flip, target removal, or hard donor
selection. The `eta=0` arm is exact equal donor weighting, and any incomplete
donor row falls back to exact equal weighting.

Four non-promotable ancestry controls use the identical weighting rule with:

- total target-cell count;
- total target UMI;
- occupied-batch count;
- pooled donor-effect norm under the same fixed operator.

For future a15 integration, the K562 reliability weight applies to the
complete existing K562 slot after its exact 50:50 PerturbMatch construction.
PerturbMatch is not treated as a second biological replicate.

## Source-only pre-gate

Before any 613-target held-context truth can be opened, every applicable
donor/operator role must satisfy:

1. reliability rankings from split pairs 1–4 versus 5–8 have a strictly
   positive target-bootstrap lower 95% bound;
2. observed same-target split retrieval is strictly above the 95th percentile
   of 49 depth/batch-count/norm-matched target-pairing nulls;
3. at least 95% of targets have all eight valid pairs;
4. every target has at least six valid pairs.

The matching null forms deterministic nearest-feature blocks from rank-scaled
log target-cell count, occupied-batch count, and fixed-operator pooled-effect
norm, then independently deranges labels within those blocks for every draw.
The same 49 draws also produce a per-target reliability null distribution;
only targets strictly above their own null p95 receive a nonzero reliability
rank.

The corrected plan also freezes 49 per-target donor-assignment nulls for each
outer donor set. Those nulls preserve each target's donor reliability values
while breaking their donor identities, and every later draw must repeat
nested eta selection.

## Fail-closed stages

```bash
# Lightweight metadata-only plan.
uv run python experiments/142-batch-reliability-weighting/run.py plan

# Pre-outcome source-only freeze. It requires the additive source-only
# authorization and atomically binds implementation, parent manifests,
# fixed operators, corrected donor-wide metadata arrays, and raw-source hashes.
uv run python experiments/142-batch-reliability-weighting/run.py freeze

# Future heavy scan, explicitly one donor at a time and under an 18 GiB guard.
uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib 18 --min-free-percent 15 -- \
  uv run python experiments/142-batch-reliability-weighting/run.py \
  source-build --donor k562_gwps --max-rss-gib 17.5

uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib 18 --min-free-percent 15 -- \
  uv run python experiments/142-batch-reliability-weighting/run.py \
  source-evaluate --donor k562_gwps --max-rss-gib 17.5

# Repeat source-build/source-evaluate for the two Orion donors, then:
uv run python experiments/142-batch-reliability-weighting/run.py source-gate
```

Every source build, source evaluation, and source gate also acquires one
experiment-wide exclusive lock. A concurrent heavy stage fails before opening
its source or stage arrays.

`develop` is an interface boundary, not an authorized stage. It refuses
before opening any held-context outcome unless a separately added
`DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json` exists, declares the 613-target
stage authorized, binds the exact passing source-gate manifest, and keeps the
410-target split, counts, challenge build, and submission disabled. No command
in this directory can create that authorization file.

## Reserved independent confirmation

Two untouched sources are reserved only for an optional confirmation after a
613-target candidate has been frozen:

- `data/external/scperturb/XuCao2023.h5ad` — HEK293 CRISPRi, approximately
  185 shared non-panel targets;
- `data/external/scperturb/TianKampmann2021_CRISPRi.h5ad` — neuronal CRISPRi,
  approximately 141–166 shared targets depending on the final eligibility
  contract.

Their expression matrices are not read by the plan, source build, source
evaluation, or source gate. Confirmation requires a separate additive
authorization and candidate-freeze artifact. Xu/Cao is a HEK293-family
surface, so any Xu/Cao confirmation must exclude Orion HEK293T from donor
construction, reliability estimation, scaling, and operator fitting. These
reserved surfaces do not block or alter the current pre-outcome experiment.

## Prohibited

- reading Experiment 132's sealed 410-target effect values;
- reading Xu/Cao or Tian/Kampmann expression values before a separately
  authorized frozen-candidate confirmation;
- running the 613-target development outcome stage without a later additive
  authorization;
- fitting a whitener on split-derived rows;
- `cell-eval2`;
- count generation;
- challenge build;
- submission;
- writes outside this experiment directory.
