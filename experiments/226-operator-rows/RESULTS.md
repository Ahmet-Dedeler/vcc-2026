# 226 — Is the whitening operator sample-limited? **No. It is SATURATED in `n`, and the knee is at `n/p ≈ 1.4` — production sits at 1.638, past it.**

**Date:** 2026-09-04. Measurement only. **No `vcc`, no `cell-eval2`, no build, no submission, no
production change, nothing written outside `experiments/226-operator-rows/`.** 19
eigendecompositions of a 17,942² covariance, run strictly one at a time.

**Memory, reported honestly.** The scripts carry a hard guard that aborts above **11.0 GiB** of
*current* RSS. **It fired once**, on the first attempt at the curve, at **12.30 GiB** — fancy
indexing a memmap pulls every touched page into the process's RSS and a blocked pass over the
4.2 GiB basis added all of it on top of the 2.6 GiB covariance. The reader was changed to plain
`read()` (page cache, not process RSS) and no run has exceeded the guard since. macOS `ru_maxrss`
high-water still reports up to **12.82 GiB** because it counts the file-backed pages of the basis
*write*; `scripts/build_union_whitener.py`'s own docstring documents why that figure overstates.

---

## Decision, in one line

**More fitting rows do not help.** Lifting production's row count by **+36%** (29,394 → 39,921,
same three donors, same prune rule, composition held fixed) moves the PLOLO transfer cosine
**+0.57% macro / +2.13% aggregate — +0.49 σ / +1.01 σ of the operator's own subsample noise —
and costs `pds` −0.001616 raw (−0.000602 overall).** The last untested dimension of the operator's
geometry is closed: after 193 (regulariser) and 212 (exponent), **the fitting sample is closed
too.** The row axis has exactly one real feature and it is a **cliff, not a slope**: below
`n = p = 17,942` the operator collapses (cosine −35% to −53%, `pds` to 0.586); above `n/p ≈ 1.4`
it is flat to the noise floor.

**Consequence for 225's ladder: the two do not compete and they do not compound.** Every source
added to the stack is *also* fitting data for the operator, and this experiment prices that
by-product at **≈ 0.5 σ**. 222's falsifier 1 found that H1 joined the stack while the operator's
sha was unchanged — H1 entered as a residual source and contributed zero operator rows. **That
was the right outcome. Nothing was left on the floor.**

---

## 0. Phase 0 — known-answer rows. **Seventeen gates, all PASS.**

Three rebuild gates (§0.1) plus fourteen in `gamma226.py` (§0.2 — twelve scalars against 212's
published rows, plus the two shape claims). An eighteenth fell out of falsifier 5 unplanned (§3.5).

Rule 11 in force: every operator is loaded by path and asserted by sha256.

### 0.1 The rebuild path reproduces the shipped operators exactly

| # | gate | expected | measured | |
|---|---|---|---|:--|
| KA-1 | re-fit for the A/B/C bundle, `--basis-prune essential` | sha `e2f65ecf69b7e0b0…` | **identical** | **PASS** |
| KA-2 | re-fit, `--basis-prune none` | sha `5608e807230dea04…` | **identical** | **PASS** |
| KA-3 | this experiment's own `f = 1.0` eigensystem vs the shipped `.npz` | bit-identity | `evals` **True**, `evecs` **True** | **PASS** |

`DEF-RUNBOOK` §4 lists KA-1 as the Oct-22 known-answer row; it passes today on today's inputs, so
the rebuild path is sound and everything below counts. KA-3 is the stronger one — the curve's own
code path, not the builder's, and it reproduces the shipped operator to the last bit.

### 🔴 Correction to the mandate (and to how the operator is usually described)

The mandate names the shipped operator as `data/donors/whitener_union.npz`, sha `e2f65ecf…`, basis
**(33,787 × 17,942)**. Those are **two different files**, and both were re-fitted here to settle it:

| file | basis rows | sha256 | what it is |
|---|---:|---|---|
| **`whitener_union_pruned.npz`** | **29,394** | **`e2f65ecf…`** | **ships (a15)**; `n/p = 1.638` |
| `whitener_union.npz` | 33,787 | `5608e807…` | v11's operator, superseded |

So production's `n/p` is **1.638, not 1.883**, and the mandate's "roughly n ≈ 2p" is closer to
**n ≈ 1.6p**. This is 193's §0 correction confirmed by rebuild rather than by inspection. It does
not change the question; it makes production *nearer* the cliff than the mandate assumed, which is
the direction that would have mattered.

### 0.2 193's vertex and 212's exponent result, reproduced on this experiment's instrument

`gamma226.py`, shipped operator asserted, `kolf2`, both renders, 212's own 15-point γ grid.

| # | gate | expected (source) | measured | \|d\| | |
|---|---|---:|---:|---:|:--|
| KA-4 | γ 0.65, a15 `pds` | 0.611201 (193 §4 / 212) | **0.6112009** | 1.1e−07 | **PASS** |
| KA-5 | γ 0.65, a15 `cos` | 0.0078207 (212) | **0.0078207** | 1.1e−09 | **PASS** |
| KA-6 | γ 0.65, a23 `pds` / `cos` | 0.610128 / 0.0075590 | identical | ≤2.6e−08 | **PASS** |
| KA-7 | γ 0.30, a15+a23, four scalars | 212's rows | identical | ≤2.8e−07 | **PASS** |
| KA-8 | γ 0.00, a15+a23, four scalars | 212's rows | identical | ≤4.9e−07 | **PASS** |
| KA-9 | cosine **monotone decreasing** in γ, all 15 points | 212 §5 | holds; 0.0092404 → 0.0075232 | — | **PASS** |
| KA-10 | a23 `pds` **maximised at γ = 0.80** | 212 §5 | argmax = 0.80 | — | **PASS** |

> **One gate had to be restated, and the restatement is a small correction to 212's summary.**
> I first wrote KA-10 as "a23 `pds` is strictly monotone increasing in γ" — 212 §5's wording. On
> `kolf2` alone it is not: there is one inversion at γ = 0.70 (0.610128 → 0.610088). **That
> inversion is in 212's own published `work/plolo212.json`**, and my run reproduces both numbers
> exactly, so the instrument agrees with 212 bit-for-bit and only my restatement of its claim was
> too strong. 212's claim is about its four-fold aggregate; it should not be read per fold.

### 0.3 The yardstick this experiment is judged against

212 closed γ by a **ratio**, and said why: *"a ratio is instrument-independent."* Along γ the
operator buys cosine and `pds` pays. Measured here, same fold, same code:

```
GAMMA axis, a23 render, wide baseline γ 0.50 -> 0.80:   d(pds)/d(cos) = -5.34
                        γ 0.60 -> 0.70:                  d(pds)/d(cos) = -1.87
```

**If more rows buy cosine at γ's exchange rate, "more rows" is γ in disguise and 212 has already
closed it.** That is the test the row axis has to pass, and it is why cosine alone was never going
to be a sufficient readout.

---

## 1. Phase 1 — the row curve. **A cliff at `n = p`, then flat.**

Composition held **exactly** fixed: the subsample is stratified by donor, so every point carries
production's donor mix to within one row; the axis is production's 17,942 columns; the renorm
reference is production's. Only `n` moves. `kolf2`, a15 config, three challenge frames — the only
leakage-free fold for a production-composition operator (162's rule, and 193 §4's).

| frac | n rows | n/p | rank-def? | floored evs | `cos` (macro) | Δ% vs prod | `cos_agg` | Δ% | `pds` | Δ raw |
|---:|---:|---:|:--:|---:|---:|---:|---:|---:|---:|---:|
| 0.250 | 7,349 | 0.410 | **YES** | 10,593 | 0.0051062 | **−34.7** | 0.0051622 | −36.0 | 0.611874 | +0.000673 |
| 0.375 | 11,023 | 0.614 | **YES** | 6,919 | 0.0036658 | **−53.1** | 0.0035355 | −56.2 | 0.586471 | −0.024730 |
| 0.500 | 14,698 | 0.819 | **YES** | 3,244 | 0.0040882 | **−47.7** | 0.0037366 | −53.7 | 0.593092 | −0.018109 |
| 0.625 | 18,371 | 1.024 | no | 71 | 0.0076129 | −2.66 | 0.0077026 | −4.48 | 0.616350 | +0.005149 |
| — | 19,961 | 1.113 | no | 52 | 0.0076520 | −2.16 | 0.0076957 | −4.57 | 0.610999 | −0.000202 |
| 0.750 | 22,045 | 1.229 | no | 55–64 | **0.0076229 ± 0.0000902** | −2.53 | 0.0077450 | −3.96 | **0.611135 ± 0.002227** | −0.000066 |
| 0.875 | 25,720 | 1.434 | no | 53 | 0.0078650 | +0.57 | 0.0080803 | +0.20 | 0.612830 | +0.001629 |
| **1.000** | **29,394** | **1.638** | no | **48** | **0.0078207** | — | **0.0080641** | — | **0.611201** | — |
| — | 29,941 | 1.669 | no | 47 | 0.0077890 | −0.41 | 0.0076858 | −4.69 | 0.608431 | −0.002770 |
| **cap lifted** | **39,921** | **2.225** | no | **36** | **0.0078653** | **+0.57** | **0.0082356** | **+2.13** | **0.609585** | **−0.001616** |

The 22,045 row is the mean ± sd of **five independent draws** (seeds 1–5). 19,961 / 29,941 / 39,921
come from the **cap-lifted** pool (all donor rows, 12,000-cap removed, same essential prune) — the
only way to extend the axis **above** production without changing composition.

### 1.1 The noise floor, measured two independent ways

| read | how | `cos` sd | `cos_agg` sd | `pds` sd |
|---|---|---:|---:|---:|
| 5 draws at n = 22,045 | seeds 1–5, stratified | **0.0000902 (1.18%)** | 0.0001693 (2.19%) | **0.002227** |
| 2 draws at n ≈ 29.4–29.9 k | production vs the cap-lifted 75% draw | spread **0.41%** | — | spread 0.002770 |

The two agree. `f = 1.0` consumes the whole pool at any seed, so its two "draws" are the same
operator and its sd is exactly 0 — reproducibility, not noise; the band comes from the subsampled
points.

### 1.2 The curve in units of its own σ

| n | n/p | z(`cos`) | z(`pds`) |
|---:|---:|---:|---:|
| 18,371 | 1.024 | −0.11 | +2.34 |
| 19,961 | 1.113 | +0.32 | −0.06 |
| 22,045 | 1.229 | −1.15 … +1.26 (5 draws) | −1.16 … +1.14 |
| 25,720 | 1.434 | **+2.68** | +0.76 |
| 29,394 | 1.638 | **+2.19** | +0.03 |
| 29,941 | 1.669 | +1.84 | −1.21 |
| 39,921 | 2.225 | **+2.69** | −0.70 |

**Read the last four rows.** From `n/p` 1.434 to 2.225 — a **55% increase in fitting rows** — z(`cos`)
goes +2.68 → +2.19 → +1.84 → +2.69. That is one number with noise on it, not a rising curve.

### 1.3 The two things the curve actually says

1. **The cliff is at `n = p`, and it is enormous.** Below 17,942 rows the covariance is singular,
   `Whitener.fit`'s `1e-12` relative floor becomes the operator (10,593 floored eigenvalues at
   n = 7,349 against **48** at production), and both channels collapse. This is the only large
   effect on the row axis and it is a rank artefact, not a sampling one.
2. **The plateau starts at `n/p ≈ 1.4` and production is on it.** The rise from the 22,045 band to
   production is **+2.60% cosine (+4.90 σ of the band mean)** with `pds` **dead flat (+0.03 σ)** —
   real, but **it is below production and cannot be harvested; production is already there.** It
   locates the knee, nothing more. Everything above 25,720 is flat.

### 1.4 The mandate's stopping rule, answered

> *"If the curve is flat from 50% to 100%, the operator is already saturated and this experiment
> ends here."*

The curve is **not** flat from 50% to 100% — but only because 50% (n = 14,698) is **rank-deficient**
(`n/p = 0.819`), which is a different failure from being sample-limited. **On the full-rank part of
the axis, which is the part the question is about, it is flat**, and it stays flat 36% *above*
production. That is the saturation result, reached by a route the mandate's 25/50/75/100 grid would
have misread.

---

## 2. Phase 2 — not triggered, and the strongest available above-production point run anyway

The mandate gates Phase 2 on "Phase 1's curve is still rising at 100%." It is not (+0.49 σ). So no
H1, CD4, 223 or 225 rows were ever put into a fitting basis — which is also the cleanest possible
answer to falsifier 1 (§3.1).

Instead I ran the one point that extends the axis **above** production with **zero** composition
change: lift the 12,000-row-per-donor cap, keep the same three donors and the same essential prune.

```
n 29,394 -> 39,921   (+10,527 rows, +35.8%)
  cos      0.0078207 -> 0.0078653   +0.57%   +0.49 sigma
  cos_agg  0.0080641 -> 0.0082356   +2.13%   +1.01 sigma
  pds      0.611201  -> 0.609585    -0.001616  -0.73 sigma   = -0.000602 overall
```

**Independent corroboration, on operators this experiment did not build.** 193's Phase 2 compared
`whitener_union` (33,787) with `whitener_union_full` (45,531) — also a pure cap-lift at fixed
composition, built by a different code path, and it stored the per-frame cosines without reporting
them. Recomputed from `193/work/phase2-basis.json`:

| comparison | Δ rows | Δ `cos` | in this experiment's σ | Δ `pds` |
|---|---:|---:|---:|---:|
| 193: `union` → `union_full` | +11,744 (+35%) | **+0.96%** | **+0.79 σ** | −0.000305 |
| **226: production → cap-lifted** | **+10,527 (+36%)** | **+0.57%** | **+0.49 σ** | −0.001616 |

Two cap-lifts, two code paths, two operator families, ~35% more rows each: **+0.5–0.8 σ on cosine
and a small negative on `pds`, both times.**

---

## 3. Phase 3 — falsifiers. **All five run, all five pass.**

### 3.1 Falsifier 1 — leak audit. *The dangerous one, and it is clean by construction.*

| object | what it is | in any fitting basis here? | scored on? |
|---|---|:--:|:--:|
| `k562_gwps` (K562) | production donor | **yes** (7,784 rows) | no |
| `orion_hct116` (HCT116) | production donor | **yes** (10,807 / 15,967) | no |
| `orion_hek293t` (HEK293T) | production donor | **yes** (10,803 / 16,170) | no |
| **`kolf2` (KOLF2.1J)** | **the evaluation surface**, 273 panel targets | **NEVER** | **yes** |
| `h1_2025`, `gse295214_h1` | H1 hESC | **never** (Phase 2 not triggered) | never |
| CD4 banks | post-whitening residual only | never | — (enters both arms identically) |
| PerturbMatch | blended into `pre` at `PM_BLEND = 0.5` | never | — (identical at every n) |

- The exclusion is **enforced in code**, not by discipline: `rows226.py` raises
  `LEAK ABORT` if the scored fold or `kolf2` appears in the donor list.
- **h1-test is not used anywhere in this experiment.** `data/donors/h1_2025` covers 85/85 of
  h1-test's targets, so an H1-fitted operator scored on h1-test would leak; the way to be safe was
  to not build one, and Phase 1's verdict meant none was needed.
- **`debench.Surface.pds()` is never imported** by any script in this directory (188 Phase 3's
  constraint). Verified by grep over `experiments/226-operator-rows/*.py`.
- The fold's prediction `pre` is byte-identical at every point on the curve — only the operator
  changes — so no n-dependent artefact can enter through the prediction side.

**One-fold limitation, stated plainly.** 162's leakage rule leaves `kolf2` as the *only* fold on
which a production-composition operator can be swapped at all, so the whole curve rests on one
surface. The mitigations are the five-draw band, the second independent draw at production's own
`n`, and 193's corroboration on a different code path — not a second fold.

### 3.2 Falsifier 2 — row-shuffle null. **Passes, and it is the most instructive result here.**

Refit on the *same* number of rows with each row's gene labels permuted independently — sample size
preserved exactly, covariance structure destroyed.

| n | | `cos` | `cos_agg` | `pds` | floored evs |
|---:|---|---:|---:|---:|---:|
| 22,045 | real | 0.0075728 | 0.0078278 | 0.608548 | 55 |
| 22,045 | **permuted** | **0.0089214** | 0.0095968 | **0.594592** | **0** |
| 29,394 | real | 0.0078207 | 0.0080641 | 0.611201 | 48 |
| 29,394 | **permuted** | **0.0090644** | 0.0098762 | **0.592464** | **0** |
| — | *212's γ = 0, no whitening at all* | *0.0092404* | — | *0.5885361* | — |

**The permuted null has HIGHER cosine than any real operator** — +15.9% at production's row count —
**and much lower `pds`**, landing almost exactly on 212's γ = 0 row. Permuting gene labels drives
`C` toward isotropy, so the permuted operator **is** approximately γ = 0.

Two things follow, and the second is the important one:

1. The null does not reproduce the real operator's row behaviour, so nothing measured here is
   "just sample size". **PASS.**
2. **Cosine alone was never a safe readout for this question.** The single largest cosine gain
   available anywhere in this experiment comes from *destroying* the operator. Any row result
   that raised cosine while `pds` fell would have been γ in disguise — exactly what 212 closed.
   The real row axis does the opposite: over the full-rank range cosine moves with `pds` flat or
   mildly negative, never trading against it at γ's −1.87 to −5.34.

### 3.3 Falsifier 3 — wrong-recipient null. **Passes with a 10–12× separation.**

184 died here (arm +0.0665, control +0.0667, null +0.0663). The control is a fixed derangement of
the target index — the *same* permutation at every row count — so `cos_wrong[i] = cos(pred_i,
truth_{π(i)})` is the same functional on a deliberately mis-paired recipient.

| n | own-pair `cos` | wrong-pair `cos` | ratio |
|---:|---:|---:|---:|
| 7,349 | +0.0051062 | −0.0004713 | 10.8× |
| 14,698 | +0.0040882 | −0.0004407 | 9.3× |
| 18,371 | +0.0076129 | −0.0008154 | 9.3× |
| 22,045 | +0.0075728 | −0.0006745 | 11.2× |
| 25,720 | +0.0078650 | −0.0006707 | 11.7× |
| **29,394** | **+0.0078207** | **−0.0007250** | **10.8×** |
| 39,921 | +0.0078653 | −0.0005423 | 14.5× |

The wrong-recipient cosine is **negative at every row count and becomes more negative** as `n`
grows across the cliff (−0.00047 → −0.00073), i.e. mis-paired predictions get *less* similar while
own pairs get *more* similar. Over the full range: own-pair **+53.2%**, wrong-pair moving the
**opposite way**. No global rescale can produce that. **PASS.**

### 3.4 Falsifier 4 — leave-top-3-out retention, on `cos` and `pds` separately

222 §9 found 14.6% vs 72.7% on the same arm — the two channels have opposite fragility. Reproduced
here, and it is what condemns the residual full-rank cosine wiggle.

| segment | channel | full Δ | minus top 3 of 273 | **retention** |
|---|---|---:|---:|---:|
| cliff, 7,349 → 29,394 | `cos` | +0.0027145 | +0.0020478 | **75.4%** |
| cliff, 7,349 → 29,394 | `pds` | −0.0006733 | −0.0029775 | 442% (sign-unstable, near zero) |
| **full-rank, 18,371 → 29,394** | **`cos`** | **+0.0002078** | **+0.0000715** | **34.4%** |
| **full-rank, 18,371 → 29,394** | **`pds`** | −0.0051489 | −0.0042847 | **83.2%** |

**The cliff effect is robust (75.4% retention); the full-rank cosine effect is not (34.4%).** Two
thirds of the small full-rank cosine difference lives in 3 of 273 targets — the signature of noise,
and independent confirmation of §1.2's σ reading. The channels' fragilities are opposite, exactly
as 222 §9 reports.

### 3.5 Falsifier 5 — panel-identity check, with the failure demonstrated on purpose

`docs/DEF-RUNBOOK.md` §0.1: the operator's axis is `union(donor coverage) − panel`, the `.npz`
stores `evals` and `evecs` and **no gene names at all**, and a new panel fails silently.
`panel226.py` recomputes the axis for three bundles and runs the guard on both a mismatched and a
matched pair:

| bundle | axis length | mispaired vs A/B/C | identical axis? |
|---|---:|---:|:--:|
| `data/challenge` (A/B/C) | 17,942 | — | — |
| 199's synthetic `def-bundle` | **17,942** | **17,729 / 17,942 = 98.81%** | **No** |
| 199's synthetic `def-bundle-r` | **17,942** | 17,511 / 17,942 = 97.60% | **No** |

- **The length is 17,942 all three times**, so length is not a discriminator — which is precisely
  why the failure is silent.
- **98.81% independently reproduces AGENTS.md §8.3's 98.8%**, on a different code path from 199's.
  That is an unplanned known-answer row and it passes.
- The guard **fired** on the mismatched pair and **passed** on the matched pair; `rows226.py`
  carries the same guard and aborts before any scoring, and the `kolf2` fold's own `wcols` was
  asserted equal to the operator axis at every run. **PASS.**

**A finding worth recording for Oct 22:** the only reason this trap exists is that the operator file
carries no axis. Nothing in this experiment's mandate permits changing the artifact, but *an
operator `.npz` that stored its own gene symbols would make the trap unrepresentable*, and the cost
is ~350 KB. Flagged for whoever owns `scripts/build_union_whitener.py`; **not done here.**

---

## 4. Pricing the arm that is on the table

The only actionable move is production → cap-lifted, because production already sits at 29,394.

```
PLOLO, kolf2 fold, a15 render, macro per-target cosine (212's functional):   +0.57%   (+0.49 sigma)
PLOLO, kolf2 fold, a15 render, energy-weighted aggregate cosine (197's):     +2.13%   (+1.01 sigma)
PLOLO pds:                                                                  -0.001616 raw
  -> -0.003575 scaled -> -0.000602 overall  (193 §5's scaled/6 x 1.011)
```

**Quoted in PLOLO units with the fold set and functional named, and NOT converted to the board by a
scalar** — AGENTS.md's 2026-09-04 warning. The two functionals disagree by 3.7× on this arm
(+0.57% vs +2.13%), which is the same macro-vs-aggregate gap that cost 216 four Phase-0 gates.

**Against §8.4's bar.** §8.4's corrected direction-arm bars are **+11.6–13.6%** cosine against the
amplitude floor and **+2.5–2.9%** shape-only, both computed at the direction-arm multiplier
**2.51×**. §8.4 also states that multiplier **does not apply to operator-geometry arms**, where the
measured range is **0.59–1.23×** because `pds` cancels 96% of the `mse` leg. This is an
operator-geometry arm, so the bar is **2.0–4.3× higher**: roughly **+5.1% to +12.3% cosine
shape-only.**

| reading | arm delivers | bar (operator-geometry, shape-only) | shortfall |
|---|---:|---:|---:|
| macro `cos` | +0.57% | +5.1 … +12.3% | **9–22×** |
| aggregate `cos_agg` | +2.13% | +5.1 … +12.3% | **2.4–5.8×** |

**It fails on the most generous reading by 2.4×, and it carries a negative `pds` leg
(−0.000602 overall) that the bar does not even charge yet. NO-GO. No slot.**

### Does it compound with 225's source ladder, or duplicate it?

**Neither — there is nothing to compound.** The by-product value of routing a new source's rows
into the operator's fitting basis is bounded by this experiment at **≈ 0.5 σ / +0.57% cosine per
+36% rows**, and the largest addition 225 could plausibly make is smaller than that. The two
channels are **independent**: a source's value as a *donor* (225's territory) is unaffected by
this, and its value as *operator fitting data* is ≈ 0. 222's falsifier 1 observed that H1 joined
the stack with the operator's sha unchanged; **that was correct, not an omission.**

---

## 5. What is now closed, and the one thing that is not

**Closed.** The operator's three parameters are now all measured and all at their optimum:

| axis | closed by | verdict |
|---|---|---|
| regulariser (`lam`, `f(ev)`) | 193 | `lam = 0`, production's map, interior optimum |
| exponent (`gamma`) | 212 (×3) | γ = 0.65, knife-edge, ratio argument |
| **fitting sample (`n`)** | **226** | **saturated above `n/p ≈ 1.4`; production at 1.638** |

**The one operational corollary that is new, and it belongs in the Oct-22 checklist.** The row
axis has a cliff at `n = p`, and nothing currently checks for it — the builder prints
`fitting (29394, 17942)` and would print `fitting (12000, 17942)` just as happily. Every operator
in service is safely above it, verified:

| operator | n | p | n/p | on the plateau? |
|---|---:|---:|---:|:--:|
| `whitener_union_pruned` (ships) | 29,394 | 17,942 | **1.638** | yes |
| `held_k562_gwps` (062) | 35,744 | 17,806 | 2.007 | yes |
| `held_orion_hct116` (062) | 27,764 | 17,942 | 1.548 | yes |
| `held_orion_hek293t` (062) | 27,554 | 17,942 | 1.536 | yes |

But a D/E/F panel that shifts the axis, or a future leakage-safe operator built from fewer donors,
could land below 1.4 **with every existing check passing** — the same failure shape as §8.3's panel
trap. **Suggested one-line guard for `scripts/build_union_whitener.py`: abort if
`total < 1.4 * len(wcols)`.** Not applied here (production is not this experiment's territory);
recorded for the runbook's owner.

---

## 6. Provenance and operating notes

**`INBOX-BUILDER-WARNING.md` (229 is editing `scripts/make_submission_pm_cd4.py`).** This
experiment **never reads, imports or copies that file**, so the warning does not bite. What it does
depend on, with shas as read (rule 11 — a path and a sha, never a description):

| dependency | sha256 (first 16) | used for |
|---|---|---|
| `scripts/build_union_whitener.py` | `2337092687eaa973` | `DonorView`, `drop_pages`; the KA-1/KA-2 rebuilds |
| `src/vcc2026/whiten.py` | `5dbddefde17192cd` | `Whitener.fit`'s floor rule (inlined, not called) |
| `experiments/162-retrodicting-validator/plolo.py` | `756903d1f40c3ec3` | the PLOLO fold, renderer, endpoints |
| `experiments/162-retrodicting-validator/run_arms.py` | `2f69489a4a18fb62` | a15 / a23 configs |
| `data/donors/whitener_union_pruned.npz` | `e2f65ecf69b7e0b0` | the shipped operator, asserted at every load |

None of these is under edit by another agent. **Nothing under `data/` was written**; the two
rebuilt operators went to `work/` and were deleted after their hashes were checked.

**`INBOX-MEMORY.md` (memory throttle).** Complied with. Every heavy job ran **strictly one at a
time** — the γ sweep and the cap-lift were chained to start only after the preceding
eigendecomposition's process had exited, never backgrounded in pairs. No sample size was reduced
and no falsifier was dropped to fit in memory; the one memory failure was fixed by changing *how*
the basis is read, not by shrinking it. Peak figures in the header are `ru_maxrss`, read in-process
via `resource.getrusage(RUSAGE_SELF)` — the same counter `/usr/bin/time -l` prints as "maximum
resident set size". At the close of this experiment `vm.swapusage` reported **791 MB free**, so no
further jobs were launched.

---

## 7. Reproduce

```bash
uv run python scripts/build_union_whitener.py --basis-prune essential --max-rss-gib 11 \
    --scratch experiments/226-operator-rows/work/scratch \
    --out experiments/226-operator-rows/work/refit_pruned.npz \
    --expect-sha256 e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123   # KA-1

uv run python experiments/226-operator-rows/panel226.py                                # falsifier 5
uv run python experiments/226-operator-rows/gamma226.py                                # KA-4..KA-10
uv run python experiments/226-operator-rows/rows226.py --fracs 1.0,0.5,0.75,0.25,0.625,0.375,0.875 \
    --seeds 1 --out work/rows226.json                                                  # the curve
uv run python experiments/226-operator-rows/rows226.py --fracs 1.0,0.75 --seeds 2,3,4,5 \
    --extra-seed-fracs 0.75 --permute-null --out work/rows226_band.json                # band + falsifier 2
uv run python experiments/226-operator-rows/rows226.py --cap 18000 --fracs 1.0,0.75,0.5 \
    --seeds 1 --out work/rows226_caplift.json                                          # above production
uv run python experiments/226-operator-rows/analyse226.py                              # falsifiers 3, 4
uv run python experiments/226-operator-rows/band226.py                                 # the merged verdict
```

Run them **one at a time** — each holds a 17,942² float64 covariance and its eigensystem. Wall
clock ≈ 2 h 20 min total on this machine with four other agents' jobs resident.
