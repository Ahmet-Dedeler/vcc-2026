# Experiment 128 results

Status: **mechanical PASS; production promotion NOT authorized**

## Verdict

The corrected integer-budget transport closes both defects identified in
Experiment 124:

1. Experiment 122 values are applied in their frozen
   `log1p(50000 * relative_abundance)` units:

   ```text
   a_j = clip(expm1(log1p(50000*q0_j) + eta*r_j), 0, inf) / 50000
   ```

2. Transport occurs after the production per-target median library `L` is
   known. Baseline quasi-emission totals are computed as
   `floor(400 * L * q0_g)`. Every non-panel, target-self, and ineligible panel
   total is copied exactly. Only the existing mutable-panel integer budget is
   redistributed, by deterministic largest remainder with canonical panel
   order breaking ties.

All required integer and scorer-space invariants passed for:

- contexts A, B, and C;
- all 300 targets;
- diagnostic doses `0.10`, `0.25`, `0.50`, and `1.00`;
- fixed `L = 19,946`, `20,034`, and `20,109`;
- the exact replayed production library picks for every context/target.

That is 14,400 complete context/target/dose/library-scheme evaluations with no
failure.

## Exact PDS-visible invariants

For every evaluated row:

- every non-panel emitted integer pseudobulk total was bit-identical;
- target-self and every ineligible panel total were bit-identical;
- the mutable panel integer budget was exactly conserved;
- full panel and grand emitted totals were exactly conserved;
- scorer-normalized non-panel linear profiles were bit-identical;
- scorer-normalized non-panel `log1p` profiles were bit-identical;
- scorer-normalized non-panel effects against the same control were
  bit-identical;
- scorer-normalized target-self linear and `log1p` values were bit-identical;
- all totals were nonnegative `int64`.

The production library-pick RNG was replayed without emitting cells. The
per-target median-library ranges were:

| Context | Control median | Replayed target-L min | median | max | unique |
| --- | ---: | ---: | ---: | ---: | ---: |
| A | 20,109 | 17,713 | 20,117 | 22,340 | 276 |
| B | 19,946 | 18,779 | 19,949.5 | 21,213 | 273 |
| C | 20,034 | 18,485 | 20,003 | 21,545 | 275 |

An independent implementation reread the controls, replayed all 900 production
library selections, reconstructed every one of the 14,400 baseline/candidate
integer-total rows, and matched every frozen row-stream SHA-256.

## Clipping and sign retention

The continuous source-unit transform retained the residual sign on every
nonzero eligible coordinate before integer-budget redistribution.

Preclip negatives were rare but real:

| Context | eta 0.10 | eta 0.25 | eta 0.50 | eta 1.00 |
| --- | ---: | ---: | ---: | ---: |
| A | 0 | 0 | 1 | 5 |
| B | 1 | 1 | 1 | 9 |
| C | 0 | 0 | 0 | 9 |

There were 27 negative coordinate occurrences across 24
context/dose/target rows. They were clipped to zero as preregistered. No
degenerate all-zero mutable-weight fallback occurred.

At `eta=1.0` under the replayed production `L`, post-integerization residual
sign retention was:

| Context | minimum | median | maximum |
| --- | ---: | ---: | ---: |
| A | 0.76632 | 0.95876 | 0.99313 |
| B | 0.63433 | 0.93843 | 0.99628 |
| C | 0.80515 | 0.95221 | 0.99634 |

Lower doses retain fewer signs after integerization because small intended
changes quantize to zero. The doses remain diagnostics only and were not
selected or ranked.

## Production boundary

This experiment validates integer pseudobulk totals, not a count artifact. A
future production implementation must pass the allocated totals directly into
the exact-total quasi scatter. Converting them back to a floating profile and
calling the existing `floor(400 * L * profile)` path again is not equivalent
and is not authorized by this result.

The mechanism is mechanically valid, but production promotion remains
**NO-GO** until a separately frozen biological confirmation evaluates this
exact source-unit transform plus integer allocation. No dose is promoted.

## Verification

- Preregistered before Experiment 128 accessed released control count values.
- Freeze SHA-256:
  `67df1b5c5e4f1960fc4ec7a1cc06d7f0754b485443b4e1b7478daefee3ac1426`.
- Synthetic tests: 8/8 passed.
- Ruff: passed.
- Python compilation: passed.
- Runtime: 12.32 seconds.
- Peak RSS: approximately 0.95 GiB.
- Results SHA-256:
  `660967723f49fde57815294b60f7f3e9820932d5ac269b78ab5888dc4b750019`.
- Artifact-manifest SHA-256:
  `70c35f37fe768aa99b58b4d7e88525e7d7979bac559840c9886e2331ecea360e`.

All declared inputs, Experiment 122 provenance links, staged-a15 metadata, and
production implementation hashes matched. The immutable artifact manifest and
all six immutable output files were independently hash-verified.

No held outcome, challenge perturbation truth, `cell-eval2`, count H5AD/VCC
build, submission, commit, or push was used.
