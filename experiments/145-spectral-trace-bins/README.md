# Experiment 145 — spectral trace-bin attenuation

Status: **terminal NO-GO**

Date: **2026-08-31**

## Question

The incumbent applies one power law across the entire source-covariance
spectrum:

```text
w_i = lambda_i^(-0.65 / 2)
```

Whitening is the largest validated retrieval gain in the project, but that
same power law amplifies the lowest-variance directions most strongly. This
experiment asks whether the tail is over-amplified.

It keeps the exact frozen Experiment-142 prediction and Experiment-140
eigenvectors. It performs no new covariance fit. The ordered source spectrum
is divided into four bins at cumulative trace fractions 25%, 50%, and 75%.
The upper half of source variance is unchanged; the lower half is smoothly
attenuated.

## Frozen transform

For incumbent post-whitening effect `B`, frozen eigenvectors `V`, and beta:

```text
Z = B @ V
multiplier by descending-variance trace bin:
    [1, 1, exp(-beta/2), exp(-beta)]
C = (Z * multiplier) @ V.T
candidate = ||B|| * normalize(C)
```

`beta` is shared across the K562 and non-K562 roles.

Frozen grid:

```text
beta = 0, 0.25, 0.50, 0.75, 1.00
```

At `beta=0`, the implementation returns a direct byte-exact copy of the
incumbent. Positive-beta rows preserve the incumbent L2 norm.

The trace boundaries are frozen from the already fitted source-only
eigensystems:

- K562 role, descending one-based ranks:
  `[1,93] [94,806] [807,2006] [2007,4992]`
- non-K562 role:
  `[1,198] [199,994] [995,2252] [2253,4992]`

## Leakage and access boundary

- `prepare` reads only frozen pre-outcome geometry, eigenpairs, folds, and
  target labels. It cannot load development truth.
- The untouched 410-target values and all reserved confirmation expression
  remain sealed.
- `prepare` freezes 99 spectral-assignment nulls and 99 within-fold
  target-alignment nulls before any development outcome access.
- Development truth requires a separate hash-bound authorization and an
  irreversible access record.
- There is no count, scorer, challenge-build, or submission command.

### Audit/refreeze history

The first prepared package, `work/preoutcome-v1`, remains preserved as
failed-audit evidence. Its candidate arrays and null plans were scientifically
valid, but its development validator would have recomputed SHA-256 over whole
donor `delta.npy` containers that also hold the sealed 410-target rows.
Development access was therefore never authorized.

The corrected `work/preoutcome-v2`:

- validates frozen truth-ledger digests and file sizes without hashing or
  sequentially reading whole truth containers;
- writes the irreversible outcome-access lock before importing truth loaders;
- loads only the authorized 613 target rows and 4,992 response columns through
  Experiment 132's memory-mapped `load_effect_rows`;
- enforces an 8 GiB process-tree RSS ceiling and 25% system-free-memory floor;
- uses one exact comparator for observed and bootstrap selection;
- validates fold/block containment and bootstrap ranges before preparation;
- reproduces v1's candidate arrays and all 198 frozen null plans exactly.

Its source-only preparation peaked at `0.99 GiB`. The corrected package
subsequently passed independent audit and one authorized 613-target
development evaluation. The final selector chose identity, with nested PDS
`-0.00006842` and signed credit `+0.00203827`; the exact family is closed.
See `RESULTS.md`.

## Nested selection

For each of five outer target folds:

1. score every beta on the other four folds against the inner-only gallery;
2. require nonnegative PDS and signed-credit gain in every biological family;
3. select by macro PDS, then macro signed credit, then smaller beta;
4. evaluate the selected beta on the held fold against the complete
   613-target gallery.

The full-development selector uses the same rule. Identity is always eligible.

## Uncertainty and nulls

The 20,000 frozen target-block bootstrap draws repeat all five fold selectors:
each draw reweights the precomputed inner-fold per-target scores, reselects
beta, and then aggregates the corresponding held-fold complete-gallery
scores. It does not merely resample one already-selected OOF vector.

If and only if the primary development gate passes, a separate null stage may
run from the already-frozen plans:

- **spectral assignment:** randomly assign the four multiplier bins to
  eigenvectors while preserving exact bin sizes;
- **target alignment:** derange candidate-minus-incumbent increments among
  targets inside each outer fold, using the same permutation across roles and
  beta arms.

Both nulls must repeat nested and final candidate selection. A primary pass
therefore means only `PASS_FOR_FROZEN_NULL_STAGE`.

## Primary promotion gate

- final beta is nonzero;
- nested and fixed macro PDS gains are each at least `+0.003`;
- nested and fixed macro signed-credit gains are each at least `+0.003`;
- reselected-bootstrap lower 95% bound is positive for both endpoints;
- every family and every family-by-fold cell is nonnegative for both;
- no family, target block, or outer fold contributes more than 50% of
  positive gain;
- every candidate value is finite and row norms are preserved;
- beta-zero is byte-identical to the incumbent.

Any failure is terminal for this exact family. No beta refinement or
family-specific routing is allowed after the result.
