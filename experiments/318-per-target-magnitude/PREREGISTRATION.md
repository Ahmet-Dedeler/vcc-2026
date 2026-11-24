# 318 — PREREGISTRATION. Per-target magnitude redistribution, isolated from the weighting.

Written and hashed **before any candidate was built and before any cloud score existed**.
`work/PREREGISTRATION.sha256` carries the digest of this file. Nothing here is
competition-facing: no submission is built, no slot is implicated, `vcc submit` is never run.

## 0. What step 1 already measured (this is INPUT, not a result to be selected on)

`profile318.py` recovered the shipped arm's **delivered** per-target magnitude profile from the
two board-scored builds on disk — `Q519I2pLWRCW7PAEx5ul` (shipped, 0.2177189) and
`2K5ppxkowMSj2ugoZqYb` (parity, 0.2169733) — as

    m(t,c) = || log-effect delivered for target t in context c, shipped ||
           / || the same, parity ||

on the real 300-target panel × 3 contexts, against the challenge's own control profiles.

| | measured |
|---|---|
| median / mean `m` | 1.01135 / 1.01443 |
| range | **0.9321 – 1.2349** |
| sd(log m) | **0.03289** |
| whole-matrix norm ratio | **1.008550** |
| reproducibility across the 3 contexts | Pearson(log m) **0.9992–0.9995**, ICC ≈ 1.000 |
| row cosine shipped↔parity | median 0.99838, min 0.9485 |

**The mandate's "0.806–2.155, whole matrix +6.0 %" is a fixture number** (316's `test316.py`
on 313's 20-target × 60-gene synthetic fixture). On the real panel the distortion is 6× smaller
in spread and 7× smaller in mean. Everything below is priced on the measured profile.

## 1. The hypothesis, stated so it can fail

**H.** Redistributing delivered magnitude ACROSS targets — holding the direction of every row
and the emission policy fixed — moves the six scored members by more than the same multiset of
multipliers assigned at random, and in the direction that would explain the shipped arm's
+0.00075 over parity.

## 2. Arms (declared before the build; `build318.py` `CONFIGS`)

`m(t)` is applied as `profile_from_effect(..., effect_scale=m[t])` — the exact point where
317's amplitude knob acts, per target instead of globally. Effect matrix, seeds, emitter,
spread policy and call set are byte-identical to `base` in every arm.

`q(t)` = the measured 300-value `log m` profile transferred to h1-test's 85 targets by matched
quantile of the target's own effect row norm ρ(t), **descending**. Justification: over the 300
real panel targets `log m` is monotone in ρ at Spearman **−0.414** (a target with a smaller
predicted response is amplified) and |log m| tracks the shipped↔parity row rotation at
Spearman **+0.650**, which is the analytic mechanism (§3). The transfer takes that relation at
its **strong form** (Spearman = −1 by construction) — a choice that can only help H.

**Every arm is renormalised to ONE whole-matrix effect norm — 1.0085497, the real shipped
arm's own measured value.** A common scalar preserves sd(log m) exactly, so the shuffled
controls stay entropy-matched to machine precision (gate 6), and every contrast among the
redistribution arms is then **pure redistribution** with the uniform component held identical.

| arm | `m(t)` | sd(log m) | what it tests |
|---|---|---|---|
| `base` | 1 | 0 | the baseline; gate 5 asserts it is 221's `rc_callset` byte for byte |
| `ptunif` | 1.008550 flat | 0 | **317's uniform axis** — the whole-matrix component, alone |
| `ptrule` | `exp(q(t))`, renormalised | 0.030964 | **the isolated arm** — the measured profile, replayed |
| `ptshuf1/2/3` | the same multiset × one scalar, permuted (seeds 1/2/3) | 0.030964 | the entropy-matched control |
| `ptamp6` | `exp(6·q(t))`, renormalised | 0.185783 | **power probe** — the same profile, 6× in log |

## 3. Frozen predictions

1. **`pds` cannot pay.** `pds` is a per-target rank statistic and a positive per-target scalar
   cannot reorder a target's own genes; 317 measured `Δpds = +0.000000` for a UNIFORM ×1.20 on
   this scorer. **Prediction: |Δ raw `pds`| < 0.0005 for `ptrule` and every shuffle, and
   < 0.004 (317's whole-axis range over a 2.9× span) for `ptamp6`.**
   *Consequence if it holds:* the shipped arm's board `pds` +0.00225 — **half** of the
   +0.00075 the mandate attributes to "the distortion" — is NOT a magnitude effect. It is the
   direction difference between applying the weighting inside vs outside the norm match
   (median row cosine 0.99838, min 0.9485), i.e. a SHAPE change that 316's parity arm removed
   along with the norms.
2. **The `mse` cost is second-order and positive.** A zero-mean spread on a convex member
   raises it. From 317's own curve (`mse` 0.930874 / 1.015978 / 1.234248 at s = 0.85/1.00/1.20)
   the local curvature is ≈ 3.0 raw per unit s², so **predicted Δ raw `mse` ≈ +0.0016 for
   `ptrule`** (≈ −0.00028 overall) **and ≈ +0.058 for `ptamp6`** (≈ −0.0102 overall).
3. **Separability.** `ptrule` is separable from its own shuffles only if
   |Δmember(`ptrule`) − mean Δmember(`ptshuf*`)| > 2 × sd(Δmember over the 3 shuffles) on at
   least one member. Declared before the numbers exist.

## 4. Decision rule

- **REAL AND PRICED** — `ptrule` beats every shuffle on the priced overall composite by more
  than 2 shuffle-sd, AND the sign is a gain, AND `ptamp6` extends it monotonically.
- **NOT SEPARABLE FROM ITS SHUFFLED CONTROL** — `ptrule` sits inside the shuffle spread.
- **NO-GO** — the composite is negative, or `ptamp6` (6× the measured dose) still cannot reach
  a fraction of the +0.0089 gate an amplitude-changing arm must clear (AGENTS.md §8.4, §3).

`ptamp6` is the power clause: **if 6× the measured redistribution cannot move the composite by
a meaningful fraction of +0.0089, the axis is closed with power, not for want of it.**

## 5. Instrument and its stated limits

The real `cell-eval2` 0.16.0 `--preset vcc2026` on `data/harness/h1-test.h5ad`, one disposable
EC2 worker per arm, RAW six members from `agg_results.csv`, converted to board overall points
through 317 §1.4's fitted scaling map. **Reused known-answer rows: 317 §1.2 (the recorded
`rc_callset` six-member row reproduced bit for bit on a fresh worker, max |diff| 4.44e-16) and
317 §1.3 (`base` sha256 = `048b9286…`, 10/10 workers).** `base` is not re-scored: `cell-eval2
run` reads only (prediction, reference, version, preset) and 317 proved the row reproduces from
the same object; my `base` is gated to the same sha256, so 317's row IS `base`'s row. That
saves one $1.72 worker inside a $12 cap.

**h1-test's limits, carried forward from 317 §2.2 and not repeated as errors:** it is scored
with the H1 union OFF (oracle-leaked 85/85), its raw `mse` at the shipped point is 1.016
against the board's 0.81671, and **it cannot locate the board's amplitude optimum.** No
absolute `mse` level here is read as a board opportunity. What it can answer is the question
actually asked: whether the ASSIGNMENT of a magnitude profile across targets moves the members
more than a random assignment of the same multiset.

## 6. Budget and safety

≤ 6 workers × ≈ $1.72 = **≤ $10.5**, cap $12. No local `cell-eval2`, ever. `uv run python`.
Writes only inside `experiments/318-per-target-magnitude/`. No submission, no `submissions/`,
no edit to `scripts/make_submission_pm_cd4.py` or `scripts/def_pipeline.py`.
