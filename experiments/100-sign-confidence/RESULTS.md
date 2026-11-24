# 100 — cross-fitted downstream sign confidence

**Verdict: NO-GO.**

The rank-8 model contains real cross-context signal, but it does not beat the
donor-agreement baseline on the preregistered weighted-precision endpoint and
it loses PDS to identity on HCT116. Rank 16 does not repair the failure.

No count-level scorer, AWS worker, challenge artifact, submission, or commit
was used.

## Protocol actually run

- Four deterministic whole-target folds; no random `(target, response-gene)`
  split.
- Strict outer donor holdout:
  - H1 excluded from every H1 fit;
  - HCT116 excluded from every HCT116 fit;
  - HEK293T excluded from every HEK293T fit.
- Training labels came only from the other two reference-DE caches.
- Each training context was predicted without both its own truth donor and the
  outer held donor.
- All 300 VCC panel genes were zeroed/excluded as response coordinates.
- Models: agreement-only ridge logistic, rank-8 STRING ridge logistic, and
  preregistered rank-16 sensitivity model.
- Fixed intervention: positive confidence reweighting with `beta=0.15`, no
  sign flips, and exact per-row L2 norm restoration.
- Twenty complete target-label permutation nulls.
- Two thousand target-block bootstrap draws.

The three real candidates and every null model were frozen before the held
surface's perturbation values were read. The source sets and target-fold
hashes are recorded in each `work/<surface>/freeze-audit.json`.

## Weighted downstream-sign precision

Precision is measured on the top 20% of each target's eligible
reference-significant genes, using experiment 098's linear-credit weights.
Identity is shown over all eligible pairs because it has no confidence
ranking.

| surface | DE-scored targets | identity, all | agreement, top 20% | rank 8, top 20% | rank 8 − agreement, 95% target bootstrap |
| --- | ---: | ---: | ---: | ---: | ---: |
| H1 | 56 / 85 | 0.60355 | **0.69552** | 0.68053 | **−0.01500** `[-0.03821, +0.00544]` |
| HCT116 | 127 / 132 | 0.67506 | **0.74765** | 0.73646 | **−0.01119** `[-0.03531, +0.01588]` |
| HEK293T | 125 / 132 | 0.48876 | 0.57824 | **0.59528** | **+0.01704** `[-0.01018, +0.04440]` |
| macro | — | — | — | — | **−0.00305** `[-0.01752, +0.01090]` |

The rank-16 sensitivity arm was also negative against agreement on H1
(`−0.02126`) and HCT116 (`−0.01336`), while positive on HEK293T
(`+0.03261`). Its macro delta was `−0.00067`
`[-0.01534, +0.01315]`.

This fails the primary endpoint on two of three surfaces and its macro lower
bound is not above zero.

## Experiment-098 linear NMAE proxy

Lower is better:

| surface | identity | agreement | rank 8 | rank 8 − identity | rank 8 − agreement |
| --- | ---: | ---: | ---: | ---: | ---: |
| H1 | 0.93062 | 0.91656 | **0.91362** | −0.01700 | −0.00294 |
| HCT116 | 0.87011 | 0.85163 | **0.84863** | −0.02147 | −0.00300 |
| HEK293T | 0.99512 | 0.98663 | **0.98644** | −0.00868 | −0.00020 |

The rank-8 macro linear-credit gain over identity was **+0.01572**
`[+0.00862, +0.02372]`. Its incremental gain over agreement was only
**+0.00205** `[-0.00022, +0.00483]`, below the project's `0.003`
resolution rule and with a confidence interval crossing zero.

Rank 16 was slightly stronger: macro gain over identity `+0.01638`
`[+0.00927, +0.02404]`, and over agreement `+0.00271`
`[+0.00015, +0.00596]`. It still fails the precision and HCT116 PDS gates.

## Complete-gallery PDS

| surface | identity | agreement | rank 8 | rank 8 − identity | rank 8 − agreement |
| --- | ---: | ---: | ---: | ---: | ---: |
| H1 | 0.86443 | 0.87227 | **0.87465** | +0.01022 | +0.00238 |
| HCT116 | **0.79436** | 0.78932 | 0.79268 | **−0.00168** | +0.00335 |
| HEK293T | 0.72282 | 0.73387 | **0.73583** | +0.01301 | +0.00197 |
| macro | — | — | — | **+0.00719** `[+0.00187, +0.01290]` | +0.00257 `[-0.00061, +0.00570]` |

The macro is positive, but the preregistered rule required a nonnegative
identity-relative delta on every surface. HCT116 is negative and the absolute
change is below the `0.003` resolution threshold.

These PDS values use every target in each gallery and exclude the union of the
300 VCC panel coordinates and that harness's target coordinates. They are
within-experiment comparisons under this stricter axis, not replacements for
the repository's official production PDS values.

## Permutation null

The observed rank-8 macro precision delta over agreement was `−0.00305`.
The 20 target-label permutation nulls had:

- median: `−0.08416`;
- 95th percentile: `−0.05736`;
- one-sided empirical `p = 0.0476`.

So the added features are not behaving like arbitrary target labels. They
carry real target/response structure. The problem is practical: that
structure is weaker than the existing donor-agreement ranking on H1 and
HCT116 and does not satisfy the all-surface safety gate.

## Decision table

| preregistered condition | result |
| --- | --- |
| rank-8 precision no worse than agreement on all surfaces | **FAIL** |
| rank-8 linear proxy no worse than identity on all surfaces | PASS |
| rank-8 PDS no worse than identity on all surfaces | **FAIL** |
| macro precision bootstrap lower bound above zero | **FAIL** |
| observed precision gain above null p95 | PASS |
| rank 16 does not reverse on any surface | **FAIL** |
| overall | **NO-GO** |

## Interpretation

The experiment found a weak but coherent **continuous reweighting** signal:
both STRING models improved linear sign credit over identity on all three
surfaces, and improved PDS over the agreement-only arm on all three. That is
different from demonstrating a better sign-confidence ranking. The stronger
agreement baseline still selected more accurate top-quintile genes in H1 and
HCT116.

The useful residue is narrow:

- STRING/context features may help a future objective trained directly on
  signed linear credit rather than binary correctness.
- Any follow-up must retain identity as a hard per-surface PDS constraint;
  averaging the positive H1/HEK293T gains would hide the HCT116 loss.
- This arm should not be built or submitted as-is.

## Numerical correction during scoring

The first H1 scoring pass produced `NaN` only in the linear proxy because
float32 reload roundoff made a few reconstructed zero-expression probabilities
slightly negative before `log2`. The frozen predictions were not changed.
The scorer was corrected to clip reconstructed probabilities at zero, exactly
matching the profile construction. The invalid first report is retained as
`work/h1/metrics.invalid-preclip.json`; all tables above use the corrected
`work/h1/metrics.json`.

## Reproduction

```bash
uv run python experiments/100-sign-confidence/test_run.py
uv run python experiments/100-sign-confidence/run.py prepare

uv run python experiments/100-sign-confidence/run.py freeze --surface h1
uv run python experiments/100-sign-confidence/run.py freeze --surface hct116
uv run python experiments/100-sign-confidence/run.py freeze --surface hek293t

uv run python experiments/100-sign-confidence/run.py score --surface h1
uv run python experiments/100-sign-confidence/run.py score --surface hct116
uv run python experiments/100-sign-confidence/run.py score --surface hek293t

uv run python experiments/100-sign-confidence/summarize.py
```

Primary machine-readable result:
`work/summary.json`.
