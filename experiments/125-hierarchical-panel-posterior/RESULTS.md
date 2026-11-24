# Experiment 125 results

Decision: **NO-GO**

The strictly nested hierarchical random-effects estimator improved the
top-20% truth-magnitude-weighted sign precision but materially reduced signed
credit versus Experiment 122 `equal_mean_agreement`. It therefore fails the
frozen promotion gates and must not be promoted.

## Primary deltas versus Experiment 122

| Endpoint | Macro delta | Target-cluster bootstrap 95% interval |
|---|---:|---:|
| Signed credit | **-0.015049899010796045** | [-0.02007071489989843, -0.01014371039201259] |
| Top-20 weighted sign precision | **+0.012174572820953597** | [+0.006403098611907555, +0.01791314457069202] |

Family deltas:

| Held family | Signed credit | Top-20 precision |
|---|---:|---:|
| Replogle K562 | -0.01214610823191394 | +0.008779778202899283 |
| Orion | -0.026082415627716488 | +0.0206978642163938 |
| CD4 | -0.006921173172757709 | +0.007046076043567707 |

All three families gained precision, but all three lost signed credit.

## Complete-procedure nulls

| Null | Credit p95 | Credit p | Precision p95 | Precision p |
|---|---:|---:|---:|---:|
| Target-label, 200 draws | +0.0026947370161614047 | 1.0 | +0.015527271023240793 | 0.3333333333333333 |
| Response-label within expression bin, 200 draws | +0.003153097401725953 | 1.0 | +0.013516444470243838 | 0.12935323383084577 |

The precision gain has a positive bootstrap interval but does not exceed
either complete-procedure null p95. The credit result is decisively negative.

## Selection and invariants

- `eb_broad_prior__abundance_exponential`: selected in 14/15 outer
  target-fold evaluations.
- `eb_conservative__abundance_exponential`: selected in 1/15.
- Production-log1p transport arms: selected in 0/15.
- Supported-pair retention: **1.0** (`231,888 / 231,888`).
- Transport invariants: **1,800 / 1,800 rows passed**.
- Expanded non-panel effects: exact zero.
- Target-self effects: exact zero.
- Preclip-negative coordinates in selected candidate plus comparator:
  **0**.

Precision gain concentration also failed the frozen family-share gate:
Orion supplied `0.5666965218078893` of positive family gain, above the
maximum `0.5`.

## Interpretation

The posterior confidence is useful for ranking signs, but the posterior mean
shrinks effect magnitude enough to lose signed credit. This repeats across all
three held families, so it is not a one-family accident. The frozen inner
selection also preferred abundance-exponential transport in every outer fold;
the exact production-log1p arm was not competitive under this estimator and
comparator pairing.

Do not integrate this estimator into production. A future experiment could
test a direction-only confidence/ranking layer while retaining Experiment
122's effect magnitudes, but that is a new hypothesis and is not authorized by
this result.

## Verification

- Ruff: passed.
- Synthetic tests: 11/11 passed.
- Immutable files verified: 16/16.
- Independent recomputation from saved per-target arrays exactly matched both
  reported macro deltas.
- No held outcomes, `cell-eval2`, counts, challenge build, submission, commit,
  or push were used.

Key SHA-256:

```text
protocol             383e4f857746ee4957b35550b03be99ef0378bc606e07f8ff9df545103147509
preregistration      ac26a83f1a8a3340e17c9a0ae77f7261e402656ebf6b09d123199fd40cd5b873
plan freeze          51d06dd0d6275dbb6ae1471a9a5f9319cde16a19039dbe01c08999e3c4d4b5e4
implementation       4e33671056c2ec6641e49b21c8acc4c4ad730734ad96e9badcb74e2994feef25
artifact manifest    b3efbd1db8a39ab933c4e4fe0bd998917216866df6157441a33a0e1560cdc63a
results              635701eefa74271e9209489b7593c045a22b573576c5a6e4319c7e7f99621010
result receipt       c7659a615b9ac4a7fb88ba53fce7002ac9e894e34f27de6738d2e6d7acc4aa86
```

Two additive implementation-recovery records preserve the original frozen
implementations. Both repairs were limited to exact compensated-sum arithmetic
and occurred before any accepted observed metric or result artifact existed.
