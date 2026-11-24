# 112 — Orion control-label correction

Status: **NO-GO by exact algebra; expression remains unopened**

Date: **2026-08-30**

## Question

Does correcting one deterministic Orion control-label bug improve the frozen
a14 family?

Both current Orion builders classify controls with a pattern ending in
`^neg`. That pattern incorrectly treats the real perturbation target `NEGR1`
as a control. A metadata-only audit found:

| line | exact eligible `Non-Targeting` | eligible `NEGR1` falsely included | legacy control total |
| --- | ---: | ---: | ---: |
| HCT116 | 165,777 | 154 | 165,931 |
| HEK293T | 218,838 | 233 | 219,071 |

`NEGR1` is not in the VCC 2026 300-target panel and is absent from both current
Orion donor perturbation axes. The target sums used by those donor artifacts
are therefore unchanged. Only the control baseline, every donor delta's
additive origin, and whitening operators fitted from those deltas are affected.

This experiment is deliberately binary:

- **incumbent:** the exact frozen a14-shaped pipeline with current Orion donors
  and current operators;
- **corrected:** the same pipeline with only the exact control-label correction
  and mandatory operator refits.

There is no candidate grid, tuning parameter, blending weight, or post-held
choice.

## Exact correction

For a donor target `p`, the current artifact stores

```text
old_delta[p] = target_log_profile[p] - old_b_ctrl
```

The corrected baseline is derived from raw rows satisfying both exact
predicates:

```text
pass_guide_filter == 1
gene_target == "Non-Targeting"
```

No regex or prefix match is allowed. The corrected donor delta is then:

```text
corrected_delta[p] =
    old_delta[p] + old_b_ctrl - corrected_b_ctrl
```

Therefore:

```text
corrected_delta[p] + corrected_b_ctrl
    = old_delta[p] + old_b_ctrl
    = target_log_profile[p]
```

No perturbation-expression row is needed to perform the correction. The
existing target sums, target order, gene order, target support, normalization
constant (`50,000`), and minimum-cell threshold remain unchanged.

## Source-only audit and extraction boundary

The first real-data stage is metadata-only. It may read only:

```text
gene_target
pass_guide_filter
```

It must bind every raw Parquet file by path, byte size, row-group count, and
SHA-256, and it must reproduce the four frozen counts in the table above. It
must also prove:

- the only accepted corrected label is exact `Non-Targeting`;
- the legacy predicate's excess is exactly `NEGR1`;
- `legacy_control_count = exact_control_count + NEGR1_count`;
- `NEGR1` is absent from `data/challenge/pert_counts.csv`;
- `NEGR1` is absent from both current Orion `perts.txt` files;
- the current `delta.npy`, `b_ctrl.npy`, `genes.txt`, `perts.txt`, and
  `meta.json` artifacts are all hash-bound.

Only after that report is frozen and reviewed may a later stage project
`gene_token_id` and `gene_expression`, and then only through two exact
predicate filters:

1. eligible `Non-Targeting` rows, to derive `corrected_b_ctrl`;
2. eligible `NEGR1` rows, solely to reconstruct the legacy mixed baseline.

Before any donor correction is accepted, the baseline reconstructed from
`Non-Targeting + NEGR1` raw sums must be bit-identical to the existing
float32 `b_ctrl.npy`. This proves that the source files, eligibility filter,
gene-token mapping, normalization, and legacy artifact all agree. The new
baseline is computed from exact `Non-Targeting` rows alone.

This scaffold does not implement or run either real-data stage.

## Corrected donor contract

For each Orion line, a future frozen correction artifact must satisfy all of
the following:

1. `genes.txt` is byte-identical to the incumbent.
2. `perts.txt` is byte-identical to the incumbent.
3. `NEGR1` is absent from both perturbation axes.
4. `delta.shape`, dtype, and row order are unchanged.
5. `b_ctrl.shape` and dtype are unchanged.
6. Every corrected row is produced only by
   `old_delta + old_b_ctrl - corrected_b_ctrl`.
7. Reconstructed target log profiles agree with the incumbent to the frozen
   float32 round-trip tolerance (`max_abs <= 2e-6`).
8. No target sum, target count, support mask, or non-Orion donor artifact is
   regenerated or modified.

## Whitening refits

Reusing an incumbent operator is forbidden. The control correction translates
every Orion donor row, so it changes the uncentered signal covariance used by
whitening.

All operator settings remain frozen at the incumbent values:

- union response axis;
- all 300 VCC panel coordinates excluded before covariance fitting;
- signal covariance, not donor-difference covariance;
- `gamma = 0.65`;
- `lambda = 0`;
- incumbent row selection, pruning, donor weighting, and deterministic seed.

The following operators must be refitted from scratch:

| use | donor rows allowed in operator |
| --- | --- |
| held HCT116 | K562-GWPS + corrected Orion HEK293T |
| held HEK293T | K562-GWPS + corrected Orion HCT116 |
| held H1 | K562-GWPS + corrected HCT116 + corrected HEK293T |
| production candidate | K562-GWPS + corrected HCT116 + corrected HEK293T |

The held line's own Orion donor is excluded from both its query and its
operator. PerturbMatch, CD4, alpha, donor weights, emission-independent
post-processing, and every other a14 component remain bit-identical.

## Frozen evaluation

The corrected donors, operators, query effects, exact axes, and complete
candidate galleries must be sealed and hash-bound before any held
perturbation truth is opened.

Compare incumbent versus corrected on:

- HCT116, strict context-held-out;
- HEK293T, strict context-held-out;
- H1 test, target-held confirmation.

For every surface:

- use the exact incumbent target axis and order;
- require the complete gallery, with no target dropping or intersection chosen
  after seeing results;
- use raw, uncentered reference effects;
- remove all 300 VCC panel-gene coordinates from the query, reference, operator
  fit, and PDS distance;
- use paired per-target PDS deltas and target-block bootstrap intervals.

HCT116 direction/DE values are **diagnostic-only** because the Orion HCT116
harness remains depth-confounded. They cannot pass, fail, rescue, tune, or
rank the candidate. Directional authorization uses only H1 and HEK293T.

## Conservative GO / NO-GO gate

Every source invariant above must pass, then every effect-level condition
below must pass:

1. complete-gallery raw PDS point delta is nonnegative on HCT116, HEK293T,
   and H1;
2. equal-surface macro raw PDS delta is at least `+0.001`;
3. the macro target-block-bootstrap 95% lower bound is above zero;
4. each surface's target-block-bootstrap lower bound is at least `-0.001`;
5. H1 and HEK293T linear directional-credit deltas are individually
   nonnegative and their equal-surface macro is at least `+0.001`;
6. H1 and HEK293T top-20% weighted sign-precision deltas are individually
   nonnegative;
7. every panel coordinate remains exactly zero and every unsupported
   coordinate remains bit-identical to the incumbent;
8. target and gene axes, target availability, row norms before whitening, and
   complete-gallery cardinalities match the frozen incumbent contracts;
9. no result depends on HCT116 direction metrics, post-held tuning, partial
   galleries, or a changed non-Orion component.

Any failure is **NO-GO**. A GO only authorizes a separate, reviewed production
patch/build experiment. It does not authorize modifying shared builders,
building a challenge submission, committing, pushing, or submitting.

## Current authorization

Authorized now:

- this preregistration;
- pure correction/algebra contracts;
- synthetic tests;
- read-only inspection of existing text metadata and code.

Not authorized now:

- reading any raw `gene_expression`;
- running the real metadata scan or source extraction;
- opening held perturbation truth;
- refitting a whitener;
- running heavy work or `cell-eval2`;
- editing shared scripts or production artifacts;
- building or submitting a challenge file;
- committing or pushing.

## Result

No real-data execution is warranted. Correcting the control baseline adds the
same vector `c = old_b_ctrl - corrected_b_ctrl` to every stored donor row:

```text
delta_corrected[p] = delta_old[p] + c
```

The production donor path immediately centers rows within each donor:

```text
delta_corrected - mean(delta_corrected)
  = delta_old + c - mean(delta_old + c)
  = delta_old - mean(delta_old)
```

The union-whitener builder applies the same centering before covariance
fitting. Production uses the centered donor effects with `beta=0`, so both
the donor query and refitted operator are mathematically unchanged. Any
measured difference would be floating-point roundoff rather than recoverable
biology.

The metadata bug is real, but it is irrelevant to the current production
predictor. No expression column, held truth, or heavy process was opened.

## Lightweight tests

```bash
uv run python experiments/112-orion-control-label-correction/test_correction.py

uv run ruff check \
  experiments/112-orion-control-label-correction/correction.py \
  experiments/112-orion-control-label-correction/test_correction.py

uv run ruff format --check \
  experiments/112-orion-control-label-correction/correction.py \
  experiments/112-orion-control-label-correction/test_correction.py

python3 -m json.tool \
  experiments/112-orion-control-label-correction/preregistered_protocol.json \
  >/dev/null
```

The tests are synthetic and perform no filesystem data discovery, Parquet
reads, held-truth access, operator fitting, or model execution.
