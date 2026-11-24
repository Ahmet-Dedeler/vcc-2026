# Experiment 124 results

Status: **NO-GO after independent audit**

## Finding

Experiment 122's frozen residual is not exactly PDS-null when naively added to
the production log-effect matrix. All 300 residual rows have nonzero sums; the
median absolute row sum is `1.45441235` and the maximum is `9.22657609`.
After production exponentiation and whole-profile normalization, every one of
the 900 context/target profiles changed at least one non-panel abundance.

At dose `1.0`, the worst non-panel abundance changes were:

| Context | Maximum absolute change | Minimum cosine to incumbent non-panel effect |
| --- | ---: | ---: |
| A | `8.32748e-6` | `0.9997959600` |
| B | `2.40553e-5` | `0.9987323310` |
| C | `7.16860e-6` | `0.9998635803` |

The numerical effect is usually small, but it is not the exact invariant
claimed by direct effect-space integration.

## Profile-level construction

The panel-simplex arm starts from the normalized a15 profile, applies the
Experiment 122 residual as a multiplicative tilt only over expressed,
off-diagonal panel genes, and renormalizes only that mutable panel slice back
to its original compensated sum.

Across contexts A/B/C, all 300 targets, and fixed doses
`0.10/0.25/0.50/1.00`, every required invariant passed:

- non-panel relative abundance is bit-identical;
- target-self abundance is bit-identical;
- panel and total mass are exactly conserved;
- all profiles are finite and nonnegative;
- non-panel log effects are bit-identical;
- each non-panel gene's pre-emission integer total is bit-identical at median
  library sizes 19,946, 20,034, and 20,109.

At dose `1.0`, median residual-sign retention after the required panel
renormalization was `0.96564` in A, `0.93284` in B, and `0.95413` in C. The
worst individual rows retained `0.67698`, `0.54478`, and `0.65074`,
respectively. This confirms why the transformed candidate itself—not the raw
Experiment 122 effect matrix—must be independently validated.

## Independent-audit correction

Two defects invalidate the claimed exact production mechanism:

1. Production quasi emission floors every gene total independently. Although
   every non-panel gene retained its own integer total, redistributing panel
   mass changed the **sum** of panel floors and therefore the emitted all-gene
   total in 282/300 A rows, 289/300 B rows, and 278/300 C rows at dose 1.0.
   The scorer normalizes by that all-gene total before excluding panel genes,
   producing maximum non-panel log-profile shifts of roughly `2.1e-6` to
   `3.0e-6`. Exact PDS invisibility is false.
2. Experiment 122's residual is in
   `log1p(50000 * relative_abundance)` coordinates. This experiment used
   `q * exp(residual)`, which interprets it as an ordinary log fold-change.
   The unit-preserving transform is
   `expm1(log1p(50000*q) + residual) / 50000` before a panel-only repair.

The preserved v1/v2 recovery files also do not hash to the historical hashes
written in their recovery records, so the failure chain is not independently
authenticatable.

The only valid conclusion is **NO-GO**. Experiment 128 started from a clean
freeze and subsequently passed exact emitted-integer and scorer-space
invariants while using the correct source units. Biological confirmation
remains mandatory.

## Verification

- Synthetic tests: 6/6 passed.
- Ruff and Python compilation passed.
- Peak RSS: 0.88 GiB.
- Immutable (but scientifically superseded) result SHA-256:
  `ea728c151f90788bc6d8bb5975344d0a7c78a5f822fabffe1a17f8f3a61bb881`.
- Artifact-manifest SHA-256:
  `ad911c5dbfe9a46e385f8cee638865ce0fe812507a64f95d409b02adc9a73f8f`.
- Final implementation-freeze SHA-256:
  `5a3fe0711c543fe2d5dc3677392658168dfbf383a0ca1ea2179e29f9002ffbb2`.

No challenge truth, held perturbation outcome, `cell-eval2`, count build,
submission, commit, or push was used. The independent audit was read-only.
