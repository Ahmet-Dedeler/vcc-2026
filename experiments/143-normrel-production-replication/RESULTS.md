# Experiment 143 — normrel production replication

Status: **terminal NO-GO**

Date: **2026-08-31**

## Decision

The fixed leakage-safe `normrel613^-1.5` candidate produced a small,
statistically positive retrieval gain, but it missed the frozen practical
threshold, harmed K562, concentrated most of its positive PDS gain in Jurkat,
and reduced signed directional credit in every biological family.

Do not retune, subset, blend, integrate, confirm, build, score, or submit this
candidate. The untouched 410-target split remains sealed.

## Primary result

| endpoint | macro delta | paired bootstrap 95% |
| --- | ---: | ---: |
| raw PDS | `+0.00279128` | `[+0.00112214, +0.00458825]` |
| signed credit | `-0.00093903` | `[-0.00165080, -0.00021027]` |

The PDS interval excluded zero, but the effect missed the preregistered
`+0.003` macro gate. Signed credit was significantly negative. Only 3 of 12
frozen promotion checks passed.

## Family result

| family | PDS delta | signed-credit delta |
| --- | ---: | ---: |
| K562 | `-0.00052911` | `-0.00061226` |
| RPE1 | `+0.00202582` | `-0.00163074` |
| Jurkat | `+0.00687714` | `-0.00057410` |

PDS positive-gain family concentration was `0.77246`, above the allowed
`0.5`; most of the gain came from Jurkat. Four of 15 PDS family-by-fold cells
were negative. Thirteen of 15 signed-credit family-by-fold cells were
negative.

## Interpretation

The old Experiment 097 normrel observation contains a narrow retrieval
signal, but it is not a robust production-weighting improvement. In the
corrected same-row, same-axis, PerturbMatch-aware geometry, inverse relative
norm shifts weight toward a Jurkat-specific PDS gain while weakening broad
directional accuracy. That is the wrong trade for a15, whose remaining
leaderboard deficit is dominated by downstream gene-direction metrics.

The candidate also cannot be rescued through its unmeasured CD4 interaction:
the development panel has zero targets with complete three-condition CD4
support, and the frozen protocol allowed a production-frame CD4 audit only
after the development gate passed. It did not.

## Integrity and access state

- The one authorized 613-target development evaluation ran once under the
  memory guard; peak observed process-tree RSS was `0.17 GiB`.
- The irreversible development-outcome access lock was written before truth
  file verification or loading.
- The 410-target confirmation values were not loaded or interpreted. The
  containing donor files were streamed bytewise after the access lock solely
  to verify their preregistered SHA-256 hashes.
- Reserved confirmation expression was not opened.
- No count generation, `cell-eval2`, challenge build, submission, commit, or
  push occurred.
- The post-run audit validates the prepared artifact, development
  authorization, result manifest, and outcome-access lock without reopening
  truth arrays.
- All 14 Experiment 143 tests pass.

Key SHA-256 values:

- evaluator:
  `ac41752ba81b545c4b61c6d4438cce0b0c9efefc64b736bc4e22f4d249df77c2`;
- protocol:
  `032ce2d73cac34d4c9032705b8d4e9987538f6493a548b844e045db156c5a68b`;
- development authorization:
  `94fe544a3cb5b00b532a16c5c46a7aabf91a2bfa8467e0979fced2970597d3c3`;
- preoutcome manifest:
  `546def6f75bb92b5da59bafbb1ed2e6e8425b77c0d646a7cb203cb997630dda5`;
- development result:
  `7f1bc6fbce304cd7fbd556e385478918d39c3e4e3f95ecbd9e288381ae012e89`;
- development result manifest:
  `93f4b09333a84b771b5c8c8ecd0981f292b18f5b40f21b2cc7a08431311f2576`;
- outcome-access lock:
  `a5ef3e87891607ced8bae138586fcd6957bb36b9d8d90a67c6603f64da3b901e`.
