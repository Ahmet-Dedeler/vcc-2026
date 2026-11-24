# Experiment 137 — high-MOI Jurkat tangent cone

Status: **terminal NO-GO on development; untouched test remains sealed**

Date: **2026-08-30**

## Question

Can Experiment 107's fixed high-MOI Jurkat source effect improve the exact a14
incumbent more reliably when it supplies only a small spherical direction
change, rather than the original Euclidean chord blend?

Experiment 107 established a real but fragile signal:

- six-gallery nested PDS delta `+0.00424563`;
- bootstrap lower 95% `+0.00127271`;
- above all 99 target-label permutations;
- one gallery slightly negative;
- directional credit negative on H1 and HEK293T.

This experiment does not refit the source model, choose a target subset from
outcomes, or revisit Experiment 107's consumed 1,800 targets. It tests one
new geometric mechanism around the exact incumbent.

## Fresh-target boundary

The inherited Experiment 107 target order is:

```text
sort(source targets intersect Nadig Jurkat targets minus 300 VCC targets,
     key = SHA256("exp107|" + target_symbol))
```

The first 1,800 entries were consumed by Experiment 107. Exactly 394 remain.
To avoid introducing a new post-hoc sort:

- untouched test = positions 1,801–1,879 in the inherited one-based order
  (`79` targets);
- development = positions 1,880–2,194 (`315` targets);
- development folds = five contiguous blocks of `63`, preserving that order.

All 79 test targets pass Experiment 107's frozen source support and
transferability gates. Three development targets (`ZNF658`, `LPA`, `NOMO3`)
are unsupported and must remain bit-exact incumbent fallbacks.

The Nadig `delta.npy` file was loaded by Experiment 107, so "fresh" here does
not mean that the file's bytes have never entered a process. It means these
394 target rows were excluded from every Experiment 107 gallery and were not
used in its selection, statistics, or result report. Experiment 137 may not
parse any of those rows until its implementation and independent review are
frozen.

## Current guardrails

- fixed source:
  `experiments/107-jurkat-high-moi/work/frozen-source-model/source-model.npz`;
- fixed source SHA-256:
  `03759a7b2cafa86d2643a58669235aca7fc38a33599aa403c67439f8979a7e26`;
- one operator/incumbent geometry excludes all 394 fresh targets and all 300
  VCC targets from both coordinates and donor-basis rows;
- unsupported targets and unsupported coordinates are exact incumbent
  fallback;
- no HCASMC, TeloHAEC, Pan, Experiment 132 test, challenge outcome, count
  generation, `cell-eval2`, submission, commit, or push.

`freeze_metadata.py` is the only authorized first stage. It reads axes,
support masks, and hashes, but no perturbation-effect matrix.

## Result

The frozen development procedure selected a final cap of `6.0°`, but failed
the promotion gates:

- nested macro PDS delta: `-0.000921659`;
- fold PDS deltas:
  `+0.000256, -0.000256, -0.005632, +0.001024, 0.000000`;
- macro signed-credit delta: `+0.000429369`, with one negative fold;
- cluster-bootstrap 95%:
  `[-0.00430108, +0.00593958]`;
- observed PDS gain did not exceed the 99-permutation maximum
  (`+0.000256016`);
- exact identity, exact fallback, and the six-degree geometry cap all passed.

The 79-target untouched test was not opened. This family is terminal and must
not be retuned against the development outcomes. See `RESULTS.md`.
