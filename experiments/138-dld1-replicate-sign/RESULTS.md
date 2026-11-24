# Experiment 138 result — terminal NO-GO

Date: **2026-08-30**

## Decision

The DLD1 replicate-consistent sparse-sign gate found a real but negligible
signal and failed the frozen development promotion gates. The untouched
410-target test remains sealed and must not be opened.

## Frozen development result

- selected final arm: `k=200`, `beta=0.05`;
- nested macro raw-PDS delta: `+0.0000995142`;
- nested macro signed-credit delta: `+0.0000602495`;
- 20,000-draw target-block bootstrap 95%:
  - PDS: `[+0.0000521183, +0.0001540774]`;
  - signed credit: `[+0.0000460219, +0.0000761047]`;
- both endpoints exceeded both frozen permutation-null p95 values;
- RPE1 PDS was negative in held folds 1 and 4;
- PDS family positive-gain concentration was `0.508929`, above the frozen
  maximum of `0.5`;
- promotion: `false`.

The effect is statistically distinguishable from the frozen nulls, but the
PDS gain is about 30 times smaller than the `+0.003` promotion threshold and
the signed-credit gain is about 50 times smaller. This route is not large
enough to matter for production.

## Integrity

- all 14 source-safe unit tests passed;
- Ruff and diff checks passed;
- peak development RSS was `1.0924 GiB`;
- unsupported targets remained bit-exact incumbent fallbacks;
- the response axis excluded all 300 VCC targets and all 1,023 gallery
  targets;
- inner arm selection used inner-only PDS galleries, while held-fold
  evaluation used the complete 613-target gallery;
- no test authorization was created and no 410-target value was read.

## Artifact bindings

- prevalue authorization SHA-256:
  `0ab0fed573c36c02343f913a6063f03f277e9cef401d5cba2d88d5e2b026f6c0`;
- development result SHA-256:
  `9a3cb07d8b6ca0ac5f009dc050a449a19d23237049ad1ce7e97d444fef89c77f`;
- development statistics SHA-256:
  `933d1aedd000a3a9216ca6d2ddd51df9b103783ea31654695e75fc0b7c896579`;
- development artifact manifest SHA-256:
  `549ad4bd29aa68a999a593945c4505bda40f5eb467615089cbed6805f91012e6`;
- terminal candidate-freeze SHA-256:
  `66bf3a3cef3c2d1e055ca0a9e2be8ce3da0a1c830fa63020bcbb0530d2e20178`.

No count build, `cell-eval2` run, challenge submission, commit, or push was
created by Experiment 138.
