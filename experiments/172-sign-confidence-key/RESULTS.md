# 172 — a learned per-(target, gene) sign-confidence key: **NO-GO**

Date: **2026-09-02**. Measurement only. No `cell-eval2`, no build, no `vcc prep`, no
submission, nothing written outside this directory.

**Decision in one line: the ordering gap is NOT learnable from AgingLab's feature family.**
Four model classes spanning 22 to 2,900 parameters, trained leave-one-donor-out on
783,000 labelled (target, gene) pairs from three donor lines, produce keys that read
**0.2212 – 0.2429** raw `reach` on h1-test, straddling **0.23533 for `|H|`**, the single
simple key we already have. The LODO-selected class reads **−0.0109
below** `|H|` (paired bootstrap SE 0.0079, P(better) = 0.074); the best class on h1-test reads
**+0.0075 above** it (SE 0.0097, P = 0.778). **Every model is inside the noise of `|H|`, in
both directions.** The brief pre-authorised this outcome: "A key that only matches `|H|` is a
real result too — it would mean the gap is not learnable from these features."

The learned keys *do* beat the information they were built from — +0.0091 (M1) to +0.0275
(gbt) over the un-whitened consensus `|cons|`, and clearly better front precision (pooled
top-1% 0.767 vs 0.708, macro top-1% 0.887 vs 0.805). **Learning works; it just recovers what
spectral whitening already provided for free.**

Against the **learnable** ceiling that `INBOX-calibration.md` established after launch (the
H-oracle, 0.6186, not the full oracle 0.6455 — its numbers reproduce here, §4.2), `|H|` closes
**12.0%** of the gap over the production order and the best learned key closes **13.8%**. The
INBOX's own pre-registered falsifier — *"if the top-5% accuracy does not exceed ~0.75 the
information is not in these features"* — comes out negative: 0.651 → 0.674 pooled.

**Recommendation to 171: deliver `|H|` as the `rank_score`.** It is the best key measured
here, it needs no model, and it has no transfer risk. The learned keys are shipped in the same
file for a successor who wants them.

Five things came out of the work that are worth more than the headline and are flagged at the
end: a reach error bar nobody in this project had, the fact that 94.4% of h1-test's `reach`
lives in 32 targets, a leak in the permutation null the brief specified, a measurement that
**whitening costs 0.011 pooled sign accuracy on the DE pool while buying `pds`**, and a
challenge-power table in which the whole channel shrinks 4×.

---

## 0. What was measured, on what, and what the instrument reproduces

| | |
|---|---|
| surface | h1-test, a15 emission — the only surface validated against the board for `reach` on emission-law changes (157 §1, INBOX-4 §6) |
| DE table | `experiments/164-callset-artefact/de164.a15.npz` (85 × 10,779) |
| effect | `experiments/169-sharpening-six-member/work/h1test-effect.npz` (`H`, `b`, `cpm`, `gate`) |
| alignment | 10,779/10,779 genes, 85/85 targets, by name |
| endpoint | `experiments/153-de-bench/debench.py:Surface.reach(key=...)`, unmodified |
| peak RSS | 6.89 GB (`donors172.py`, memory-mapped Orion pages counted resident); 3.87 GB for the fits. Cap 10 GB, two sibling agents running. |

**Known-answer rows (AGENTS.md §11.1 rule 10). Every script in this directory calls
`harness172.check()` before it prints anything and aborts if any row moves.**

| row | here | expected |
|---|---:|---:|
| `reach(production order)` | **0.18283** | 0.18283 ✓ |
| `reach(\|H\|)` | **0.23533** | 0.23533 ✓ |
| `reach(sig-first, then \|H\|)` | **0.21914** | 0.21914 ✓ |
| `reach(ORACLE)` | **0.64553** | 0.64553 ✓ |
| INBOX-4 §1 decile table, `\|H\|` | 0.621 0.552 0.544 0.537 0.531 0.522 0.521 0.516 0.501 0.499 | 0.621 0.550 0.543 0.536 0.531 0.521 0.521 0.515 0.502 0.499 ✓ |
| INBOX-4 §1 decile table, expression | 0.540 0.551 0.539 0.533 0.535 0.532 0.519 0.537 0.529 0.529 | 0.540 0.550 0.539 0.532 0.534 0.533 0.519 0.536 0.528 0.528 ✓ |
| 157 §3b power matching at `f = 1` | 0.18283 / 0.23533 / 0.64553 | reproduces the native BH surface exactly ✓ |

### 0.1 A correction to the brief's top-X% curve

The brief quotes `|H|`'s top-X% sign accuracy as **1% .779 · 2% .729 · 5% .724 · 10% .691 ·
25% .646 · 50% .582 · all .536**. **That curve does not reproduce here under any convention I
could construct, and it contradicts INBOX-4 §1's own decile table**, which puts `|H|`'s first
decile — the top 10% — at **0.621**, not 0.691.

What does reproduce, on two independent ten-cell rows, is INBOX-4 §1 itself: pool = the
reference's BH-significant *and adjudicable* genes, ordered per target, accuracy **pooled**
over (target, gene) pairs. Under that convention `|H|` reads:

| | 1% | 2% | 5% | 10% | 25% | 50% | all |
|---|---:|---:|---:|---:|---:|---:|---:|
| **pooled** (this document's convention) | 0.707 | 0.666 | 0.651 | 0.623 | 0.579 | 0.557 | **0.534** |
| macro (mean of per-target accuracies) | 0.826 | 0.784 | 0.795 | 0.769 | 0.721 | 0.664 | 0.609 |
| the brief | 0.779 | 0.729 | 0.724 | 0.691 | 0.646 | 0.582 | 0.536 |

The brief's row is most likely a macro curve whose `all` cell came from the pooled statistic —
the macro *decile* d1 is 0.689 against the brief's 0.691. It is not chased further, because
**every key in this document is measured against `|H|` on the same lines of code**, which is
what the comparison needs. Both conventions are reported for every model.

*(This also identifies INBOX-4 §1's "`|predicted effect|`" row as `|H|`, not `|pred_lfc|`,
which reads 0.589 in its first decile.)*

---

## 1. Donor audit — the leak that had to be excluded, and the ones that did not

`donors172.py`, printed in `work/donors172.log`:

| donor | perturbations | genes | h1-test targets covered | |
|---|---:|---:|---:|---|
| **`h1_2025`** | 255 | 18,080 | **85/85 (100.0%)** | **BANNED — pure oracle (AGENTS.md trap #1; 157 §7 records this exact leak firing at reach 0.404)** |
| `k562_gwps` | 9,787 | 8,248 | 82/85 (96.5%) | training + deployment donor |
| `orion_hct116` | 17,767 | 38,606 | 84/85 (98.8%) | training donor |
| `orion_hek293t` | 17,977 | 38,606 | 84/85 (98.8%) | training + deployment donor |
| `k562_essential` / `nadig_jurkat` / `rpe1` | — | — | 18–21/85 | unused |

`h1_2025` is not opened anywhere in this experiment. **The two Orion lines cover 84/85 of
h1-test's targets and ARE used for training — that is not leakage, and the reason is
structural, not a judgement call:** the label in a fold is whether a transferred sign matched
the *held-out donor line's own* response, and the model carries **no gene identity and no
target identity**, so there is no channel through which "target X, gene Y goes up in H1" could
be memorised. What a fold can teach is only the shape of the map from measurable scalars to
sign reliability. The falsifier for that claim is section 5's strict null, which reads the
constant-key floor.

**Folds.** Each holds out one production donor and transfers from the other two exactly as
production does (`MultiDonorModel(..., gene_norm=False, eb=False)`, uniform weights,
`beta = 0.0` per `submissions/a15-cd4f010.vcc.meta.json`):

| held out | sources | note |
|---|---|---|
| `orion_hct116` | `k562_gwps` + `orion_hek293t` | **this is exactly 095's h1-test donor pair** |
| `orion_hek293t` | `k562_gwps` + `orion_hct116` | one Orion line swapped |
| `k562_gwps` | `orion_hct116` + `orion_hek293t` | both sources from the same atlas — the odd fold |

The third fold is weaker by construction: two X-Atlas/Orion lines are far less independent
than K562 + Orion, so "cross-source agreement" means less there. Its numbers are reported
separately in every table and it is consistently the worst fold.

### 1.1 The training pool, and the one approximation in it

`reach` walks the reference's DE pool, so the rows the model is fitted on should be the rows
it will be asked to order. A donor screen has no p-values, only `delta`, so the pool is proxied
by the same Poisson-noise Wald statistic `eb_shrink` already uses in production:
`z ~ |delta| · sqrt(n_cells · depth · p_g)`, which inside a target is `|delta| · sqrt(cpm)`.

**Validated on h1-test, where the true BH set IS known** (`lodo172.py:pool_check`), by taking
each target's top-`n_real` genes under the proxy and asking what fraction are really in the
pool:

| proxy | median recall of the real DE pool | mean |
|---|---:|---:|
| **`\|lfc\| · sqrt(cpm)`** (used) | **0.828** | 0.814 |
| `\|lfc\|` alone | 0.364 | 0.408 |
| `cpm` alone | 0.207 | 0.214 |
| random | 0.043 | 0.104 |

| fold | recipient axis (gated at 5 CPM) | targets | pool/target | rows | base sign accuracy |
|---|---:|---:|---:|---:|---:|
| `k562_gwps` | 8,248 | 3,000 | 68 | 204,000 | 0.5311 |
| `orion_hct116` | 11,105 | 3,000 | 92 | 276,000 | 0.5668 |
| `orion_hek293t` | 12,265 | 3,000 | 101 | 303,000 | 0.5707 |

Targets are sampled 95% carried-by-both / 5% carried-by-one to match h1-test's 81/85, and the
pool fraction (0.826%) matches h1-test's median `n_real` 89 / 10,771 gated. h1-test's own
un-whitened-consensus sign accuracy over the pool is **0.5449**, so the folds bracket it.

---

## 2. The feature map — 21 scalars, no identities

`feats172.py`. AgingLab's published four (expression, saturation, consensus magnitude,
cross-source agreement) plus the extensions the brief asked for. **`build()` is called by the
same code for a LODO fold and for h1-test**, which is the only way to rule out training/serving
skew. Every magnitude feature is a within-target midrank, a within-target log-ratio, or
centred on the recipient's own median, so the build's arbitrary `alpha · global` amplitude
cannot enter the model.

| group | features |
|---|---|
| consensus magnitude | `rank_cons` `log_cons_rel` `rank_geo` `log_geo_rel` |
| cross-source agreement | `concord` `rel_dis` `ratio_mm` `rank_agree` `n_don` |
| donor-panel responsiveness | `rank_z` `log_z_rel` `log_marg` `log_sd` `resp` |
| expression | `log_cpm_rec` `log_cpm_don` |
| saturation | `sat` = `log1p(\|cons\| / (b_ctrl + 0.5))` |
| target-level | `t_cos` `t_lognorm` `t_concfrac` `t_ndon` |

**Every rank is a midrank** (097 §10). The agreement vector is exactly zero wherever a donor is
silent or the donors disagree — a majority of most rows — and ordinal ranks would hand that
block an arbitrary ramp ordered by position in the gene axis, which is the bug 097 §10 caught
in its own per-gene key and in production's `agreement_shrink_vector`.

**Two features are deliberately absent and both omissions are load-bearing:**

1. **The whitened effect `H` itself.** Production whitens with `union_pruned`, an operator
   fitted on donor perturbations; there is no leak-free per-fold analogue (it would have to be
   refitted without the held-out line, on that line's own gene axis), and AGENTS.md §4.2's
   operational corollary says the DE ordering signal is the **un-whitened** cross-donor
   agreement anyway. `|H|` is therefore the *benchmark*, not an input. Section 6 measures what
   that costs.
2. **Recipient control-cell detection rate and dispersion.** Buildable on h1-test, not
   symmetrically buildable for the donor lines without re-deriving three raw h5ads — and
   experiment 161 measured that **exactly these two scalars make a cross-line learned model
   worse** (worst-line Δpds −0.00097 → +0.00382 when dropped). Unmeasured, with that reason.

---

## 3. LODO fold table — model selection, and it does not go where capacity says

`run172.py` §B. Trained on two donor lines, scored on the third with `debench`'s own `reach`
rule applied to the fold's DE-pool proxy. **AGENTS.md §11.1 rule 9 applies to the magnitudes
in this table** — the pool denominator (68–101) is not h1-test's and is not the challenge's, so
only the *ordering* of rows is readable. Best setting per class, chosen inside this table;
the full 24-config grid is in `work/run172.log`.

| model / setting | `k562_gwps` | `orion_hct116` | `orion_hek293t` | **mean** | front@5 | params |
|---|---:|---:|---:|---:|---:|---:|
| `\|cons\|` (un-whitened) | 0.12524 | 0.13279 | 0.12667 | 0.12823 | 0.6474 | — |
| `geo × concord` (013/095) | 0.11218 | 0.13157 | 0.14005 | 0.12793 | 0.6563 | — |
| expression (log cpm) | 0.08293 | 0.08503 | 0.09452 | 0.08749 | 0.5625 | — |
| **M1 logistic, C=0.1, w=top25** | 0.13351 | 0.15088 | 0.16255 | **0.14898** | 0.6699 | **22** |
| M4 gbt, it=100, w=1/rank | 0.13163 | 0.15378 | 0.16081 | 0.14874 | 0.6689 | 2,900 |
| M5 mlp, hidden=(16,) | 0.12958 | 0.14582 | 0.15400 | 0.14313 | 0.6680 | 369 |
| M2 step-additive, C=0.01, w=1/rank | 0.11364 | 0.14375 | 0.12439 | 0.12726 | 0.6645 | 406 |

Four things to read off it, and all four survived to h1-test:

- **Learning beats the raw consensus on every fold**: 0.12823 → 0.14898, +16%, 3/3 folds.
- **Capacity buys nothing.** The 22-parameter logistic and the 2,900-parameter gbt are
  **0.00024 apart** — noise. The AgingLab-scale additive model (406 params) is *worse* than
  the 22-parameter one at every regularisation strength tried, because the step basis spends
  its parameters on shape in a region where the response is close to linear in the
  standardised features. **M1 vs M4 is a coin flip in selection; naming M1 "the winner" is a
  formality and section 4 reports both.**
- **The top-heavy weightings mostly hurt** (constraint 4's answer). Three were tried: a hard
  top-25%-of-each-target cut, `1/rank`, and `exp(−rank/10)`, all derived from a prior model's
  ordering and never from labels. Only `top25` on the logistic helped (0.14652 → 0.14898);
  `1/rank` cost the logistic −0.011 and `exp(−r/10)` cost the additive model −0.004. Sharper
  concentration is worse: it throws away the rows that calibrate the front.
- **The `k562_gwps` fold is the worst in every row**, as predicted from its two same-atlas
  sources.

---

## 4. h1-test — every class refitted on all three folds, evaluated once

`run172.py` §C. Nothing in this table was selected here.

| key | **reach** | d1 | d2 | top1% | top5% | top10% | macro top1% | ρ(pool) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| constant key (gene-name order) | 0.18021 | 0.532 | 0.526 | 0.543 | 0.536 | 0.535 | 0.618 | — |
| production p-values | 0.18283 | 0.546 | 0.555 | 0.558 | 0.550 | 0.549 | 0.649 | 0.0264 |
| `geo × concord` (013/095) | 0.21007 | 0.641 | 0.568 | 0.754 | 0.666 | 0.642 | 0.842 | 0.0696 |
| `\|cons\|` (un-whitened) | 0.21539 | 0.630 | 0.593 | 0.708 | 0.652 | 0.632 | 0.805 | 0.0795 |
| M5 mlp | 0.22121 | 0.653 | 0.578 | 0.727 | 0.672 | 0.654 | 0.845 | 0.0732 |
| **M1 logistic — LODO-selected** | **0.22448** | 0.645 | 0.591 | 0.731 | 0.667 | 0.647 | 0.864 | **0.0852** |
| 50/50 rank blend M1 × `\|H\|` (untuned) | 0.22906 | 0.640 | 0.574 | 0.740 | 0.668 | 0.642 | 0.855 | 0.0785 |
| M2 step-additive | 0.23163 | 0.632 | 0.576 | 0.766 | 0.663 | 0.634 | 0.876 | 0.0645 |
| **`\|H\|` — THE BENCHMARK** | **0.23533** | 0.621 | 0.552 | 0.707 | 0.651 | 0.623 | 0.826 | 0.0604 |
| M4 gbt | **0.24287** | 0.641 | 0.586 | 0.767 | 0.674 | 0.644 | 0.887 | 0.0597 |
| **ORACLE** | **0.64553** | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.0000 |

Full decile and macro curves for every row: `work/run172.json`.

**The learned keys win precision at the top and lose it in the middle.** M4's pooled top-1% is
**0.767 against `|H|`'s 0.707** and its macro top-1% is **0.887 against 0.826** — the models
did learn what the brief asked them to learn. But by d4 they are *below* `|H|` (0.530 vs
0.537), and `reach`'s `k*` walk keeps going past the front, so a sharper front with a softer
shoulder is roughly a wash. That is the mechanism behind the headline.

**The ρ calibration in the brief reproduces and then some.** M1 has the highest pool Spearman
of any key here (0.0852, above `|H|`'s 0.0604) and is *worse* on `reach`. M4 has the *lowest*
ρ of the learned keys (0.0597) and the best `reach`. **Rank correlation is not merely a weak
proxy for this metric — on this table it is anti-correlated with it.**

### 4.1 With error bars, none of it is a result

`err172.py`. Paired bootstrap over the 82 scored targets, 20,000 resamples. Paired matters: the
jackknife SE of a single key's `reach` is ~0.036, while the SE of a *difference* is 0.005–0.016,
because the same 32 small-pool targets dominate both keys (§7.2).

| key | reach | Δ vs `\|H\|` | SE | 95% CI | P(better) |
|---|---:|---:|---:|---|---:|
| production order | 0.18283 | −0.05250 | 0.01550 | [−0.0838, −0.0233] | 0.000 |
| `geo × concord` | 0.21007 | −0.02527 | 0.01432 | [−0.0541, +0.0025] | 0.037 |
| `\|cons\|` | 0.21539 | −0.01994 | 0.00822 | [−0.0371, −0.0051] | 0.003 |
| M5 mlp | 0.22121 | −0.01412 | 0.00887 | [−0.0322, +0.0024] | 0.049 |
| **M1 logistic (selected)** | 0.22448 | **−0.01086** | 0.00785 | [−0.0270, +0.0035] | **0.074** |
| 50/50 blend | 0.22906 | −0.00627 | 0.00598 | [−0.0192, +0.0040] | 0.145 |
| M2 step-additive | 0.23163 | −0.00370 | 0.01205 | [−0.0271, +0.0203] | 0.375 |
| **M4 gbt** | 0.24287 | **+0.00753** | 0.00965 | [−0.0105, +0.0276] | **0.778** |
| ORACLE | 0.64553 | +0.41019 | 0.02617 | [+0.3584, +0.4605] | 1.000 |

And against `|cons|` — the same information the models were given, with all coefficients
collapsed onto one feature. **This is the question "did learning add anything":**

| key | Δ vs `\|cons\|` | SE | 95% CI | P(better) |
|---|---:|---:|---|---:|
| M5 mlp | +0.00582 | 0.00849 | [−0.0109, +0.0226] | 0.756 |
| M1 logistic | **+0.00908** | 0.00513 | [−0.0007, +0.0196] | **0.967** |
| 50/50 blend | +0.01367 | 0.00696 | [+0.0013, +0.0284] | 0.985 |
| M2 step-additive | +0.01624 | 0.01225 | [−0.0071, +0.0412] | 0.911 |
| `\|H\|` | +0.01994 | 0.00822 | [+0.0051, +0.0371] | 0.997 |
| **M4 gbt** | **+0.02747** | 0.01133 | [+0.0068, +0.0510] | **0.997** |

**The reading.** Learning adds a real +0.009 to +0.027 over the un-whitened consensus, and
`|H|` adds +0.020 over the same baseline for free. The learned keys and `|H|` are
statistically indistinguishable from each other (M1 P = 0.074, M4 P = 0.778, and M4 − M1 =
+0.01839 ± 0.00995, P = 0.971 — itself not clean). **Whatever the model learns from
un-whitened donor scalars, spectral whitening had already extracted.**

### 4.2 Answering `INBOX-calibration.md` directly

That INBOX landed after launch. Its three numbers reproduce here (`work/inbox_check.log`):

| INBOX claim | here | theirs |
|---|---:|---:|
| P(emitted sign == design sign) over the pool | **0.94521** | 0.945 ✓ |
| design-sign accuracy over the pool | **0.53733** | 0.537 ✓ |
| full oracle (emitted signs) | **0.64553** | 0.64553 ✓ |
| **H-oracle** (design signs, perfectly ordered) | **0.61862** | 0.62119 — 0.0026 apart, a tie convention in the 0/1 key |

So the **learnable** ceiling is 0.6186, not 0.6455: emission corrupts only 5.5% of pool signs,
and 96% of the headroom is a property of `H` itself. Re-scored against that ceiling, over the
production baseline 0.18283:

| key | reach | **% of the H-oracle gap closed** | % of the full-oracle gap |
|---|---:|---:|---:|
| `\|cons\|` | 0.21539 | 7.5% | 7.0% |
| M5 mlp | 0.22121 | 8.8% | 8.3% |
| M1 logistic (selected) | 0.22448 | 9.6% | 9.0% |
| 50/50 blend | 0.22906 | 10.6% | 10.0% |
| M2 step-additive | 0.23163 | 11.2% | 10.5% |
| **`\|H\|`** | **0.23533** | **12.0%** | **11.3%** |
| M4 gbt | 0.24287 | 13.8% | 13.0% |

**The best learned key closes 13.8% of the learnable gap against `|H|`'s 12.0%** — a 1.8-point
improvement that §4.1's paired bootstrap cannot separate from zero, and it belongs to a class
that was not LODO-selected.

**The INBOX's own falsification criterion is met in the negative.** It pre-registered: *"Whatever
features you add must lift the top 5% well above 0.72 … If after honest LODO the top-5% accuracy
does not exceed ~0.75, the information is not in these features and the experiment should say
so — that is a real result."*

| convention | `\|H\|` top-5% | best learned top-5% | lift | above 0.75? |
|---|---:|---:|---:|---|
| pooled (reproduces INBOX-4 §1) | 0.651 | 0.674 (M4) | +0.023 | **no** |
| macro | 0.795 | 0.831 (M4) | +0.036 | yes — **but the baseline already was** |

Under the pooled convention nothing comes near 0.75; under the macro convention `|H|` clears it
without any model. **Under either reading the learned lift is a fraction of what the criterion
asked for, so: the information is not in these features, and this experiment says so.**

The INBOX's other two suggestions were followed and are worth recording as confirmations:

* **"Your target label should be `sign(H) == sign(real_lfc)`, not the emitted sign."** This is
  what the folds already do — a fold has no emitter, so its label is design-sign correctness of
  the transferred consensus. The convergence is independent, and §6.1 shows why it matters:
  design and emitted signs agree 94.5% of the time on the pool, so the two labels are the same
  problem.
* **"Consider modelling P(wrong) explicitly … a key that pushes confidently-wrong genes to the
  back is worth as much."** A calibrated P(correct) already does both — it is the same ordering
  read from the other end. What it does *not* do is separate "confidently wrong" from
  "unknown", and the d1–d2 rows of the INBOX's table say those are the same genes here
  (`|H|` < 0.013, P(correct) 0.51, P(emission flips) 0.36). Not pursued further; the
  falsification criterion above already answers the parent question.

---

## 5. Permutation null — the briefed one leaks, and that is a finding

`run172.py` §D. Two nulls, because the one the brief specified does not do what it looks like
it does. Before trusting a control, ask what it would read if the thing being tested were
absent — and this one reads almost the same either way.

| null | mean h1 `reach` | sd | range (6 seeds) |
|---|---:|---:|---|
| **N1 — labels shuffled WITHIN (fold, target)**, as briefed | **0.20904** | 0.00501 | 0.2008 – 0.2156 |
| **N2 — labels shuffled across ALL rows** (strict) | **0.17684** | 0.01377 | 0.1631 – 0.1971 |
| for reference: fitted M1 | 0.22448 | | |
| for reference: `\|cons\|`, no learning at all | 0.21539 | | |
| for reference: constant key (gene-name order) | 0.18021 | | |

**Why N1 leaks, measured.** `corr(y, rank_cons)` pooled over all 783,000 fold rows:

| | real labels | after N1 (within-target shuffle) | after N2 (global shuffle) |
|---|---:|---:|---:|
| `corr(y, rank_cons)` | **+0.12357** | **+0.03665** | −0.00022 |

A within-target shuffle preserves each target's own base sign accuracy. `rank_cons` is *not*
uniform inside a target's DE pool — the pool is a top-`K` subset of the live genes, so targets
whose consensus agrees with the recipient's DE genes have both a higher mean `rank_cons` and a
higher base accuracy. That **between-target** association survives the shuffle at 30% of its
original strength, and a monotone fit converts it straight back into a **within-target**
ordering. So N1 reads `|cons|` (0.209 vs 0.215), not the floor, and would have looked like a
"the model is 60% noise" verdict if taken at face value.

**N2 is the honest null**, and it is also the honest measure of how much `reach` a key can
collect from nothing: **0.17684 ± 0.01377**, straddling the constant-key floor 0.18021, with a
seed range of 0.034 raw. The fitted M1 sits 3.5 sd above it. The models are real; they are just
not better than `|H|`.

**Generalisable form of the leak, for the next agent:** *any* label permutation that leaves
coarse-grained structure intact — target, fold, batch — can be converted into fine-grained
ordering by a monotone fit whenever the fine-grained feature is non-uniform within the coarse
unit. The fix is to permute at the level the metric reads, or to permute everything.

---

## 6. Ablation — which features carry it, and what whitening is doing

`run172.py` §E, on the LODO winner. Both surfaces, because they disagree and the disagreement
is the informative part.

| dropped group | h1 `reach` | Δ | mean fold `reach` | Δ | h1 top1% |
|---|---:|---:|---:|---:|---:|
| (nothing: full model) | 0.22448 | — | 0.14898 | — | 0.731 |
| expression | 0.22381 | −0.00066 | 0.14209 | **−0.00689** | 0.725 |
| donor-panel responsiveness | 0.22623 | +0.00175 | 0.14230 | **−0.00668** | 0.724 |
| saturation | 0.22622 | +0.00174 | 0.14729 | −0.00169 | **0.771** |
| target-level | 0.23042 | **+0.00595** | 0.14889 | −0.00009 | 0.734 |
| consensus magnitude | 0.22691 | +0.00244 | 0.14963 | +0.00065 | 0.733 |
| cross-source agreement | 0.22389 | −0.00058 | 0.15120 | +0.00222 | 0.730 |

**The two groups that carry the fold signal are expression and donor-panel responsiveness**
(−0.0069, −0.0067), and neither is measurable on h1-test at this resolution. **Consensus
magnitude and cross-source agreement — the two the whole prior literature points at — can be
dropped with no loss on either surface**, because `rank_z` and `log_z_rel` (consensus divided
by the gene's own across-panel sd) already contain them in a better-normalised form. The
learned coefficients say the same thing (`work/coefs.log`):

| feature | coefficient (log-odds / sd) | univariate corr(y, ·) |
|---|---:|---:|
| `log_z_rel` — consensus ÷ the gene's own panel sd | **+1.163** | +0.1656 |
| `t_concfrac` — per-target sign-agreement fraction | −1.133 | +0.0304 |
| `t_ndon` — donors carrying the target | +1.105 | +0.0129 |
| `log_marg` — the gene's mean \|effect\| across the panel | **−0.908** | +0.1153 |
| `log_cons_rel` — consensus vs the target's median | +0.834 | +0.1610 |
| `rank_z` | −0.737 | +0.1363 |
| `rank_cons` | +0.479 | +0.1236 |
| `log_cpm_rec` — recipient expression | −0.002 | +0.0292 |

**The model's single strongest statement is `+log_z_rel` with `−log_marg`: trust a gene when
this target moved it far more than the donor panel usually moves it, and distrust genes the
panel moves a lot in general.** That is a standardisation, not a magnitude — and it is exactly
the thing `|effect|` ordering cannot express. Raw recipient expression gets a coefficient of
**−0.002**: after the standardised features are in, expression carries nothing on its own,
which independently reproduces INBOX-4 §1's flat expression decile row.

The target-level block is the mirror image: it has three of the seven largest coefficients and
**dropping it raises h1 `reach` by +0.006 and changes fold `reach` by −0.0001**. Those features
are constant within a target, so they contribute exactly zero to a within-target ordering; they
are soaking up between-target base-rate variance, which lets the gene-level coefficients fit
cleanly, and then they cost a little at deployment through the intercept. **A successor should
fit them as per-target offsets and drop them from the delivered key.**

### 6.1 Whitening: it costs sign accuracy on the DE pool and buys the reach back anyway

`probe172.py`. Measured on the h1-test DE pool, 64,230 pairs with a non-zero consensus:

| | |
|---|---:|
| `sign(H) == sign(cons)` — how often whitening flips a called sign | **0.8405** |
| within-target Spearman(`\|H\|`, `\|cons\|`), median | 0.7247 |
| pooled sign accuracy, **un-whitened consensus** | **0.5449** |
| pooled sign accuracy, recovered `H` (whitened) | 0.5373 |
| pooled sign accuracy, emitted `pred_lfc` (what `reach` scores) | 0.5343 |

**Whitening flips 16% of the signs on the pool and costs 0.0076 pooled sign accuracy, and
still produces the better `reach` key (0.23533 vs 0.21539).** It orders better while signing
worse. This is the same shape as §4.2's headline — whitening *lowers* transfer cosine while
*raising* `pds` — one level down, on a different metric, and it is a new measurement.

It also puts a number on the ORACLE. The oracle is 0.64553 against our *emitted* signs; against
the un-whitened consensus signs the same construction would start from a 0.0106-higher per-gene
accuracy. **A route that recovered the un-whitened sign while keeping the whitened profile is
worth roughly 0.011/0.9 ≈ 0.012 raw `reach` at the oracle end and is not in anyone's plan.**
Flagged, not pursued: it changes the emitted profile, which is `pds`'s and 165's territory.

---

## 7. Two instrument findings that outlive this experiment

### 7.1 `reach` on h1-test has an error bar of ±0.036, and nobody had measured it

157 §7 quotes seed spreads of 0.008–0.011 for `reach`; those hold the target set fixed and
re-draw the emitter. **The jackknife SE over the 82 scored targets is 0.036** — 3–4× larger,
and it is the relevant one whenever two different keys or two different builds are compared on
different data. Paired differences are much tighter (0.005–0.016) because the same targets
dominate both sides, which is why §4.1 uses the paired bootstrap and why unpaired reach
comparisons across experiments should be treated with suspicion.

### 7.2 94.4% of h1-test's `reach` comes from 32 targets

`anatomy172.py`. `R_p = k*/n_real`, and at a front accuracy near 0.62 a ≥90%-pure prefix is a
handful of genes:

| key | reach | `k*`=0 | 1–2 | 3–9 | ≥10 | median `k*` | max `k*` |
|---|---:|---:|---:|---:|---:|---:|---:|
| production | 0.18283 | 22 | 38 | 20 | 2 | 1 | 40 |
| `\|cons\|` | 0.21539 | 11 | 38 | 22 | 11 | 2 | 193 |
| `\|H\|` | 0.23533 | 10 | 30 | 30 | 12 | 3 | 163 |
| ORACLE | 0.64553 | 0 | 9 | 20 | 53 | 66 | 3,391 |

| `n_real` bucket | targets | mean `k*` | mean `R_p` | **share of reach** |
|---|---:|---:|---:|---:|
| **1–20** | **32** | 3.16 | **0.5690** | **94.4%** |
| 21–50 | 3 | 2.00 | 0.0800 | 1.2% |
| 51–100 | 6 | 5.50 | 0.0661 | 2.1% |
| 101–300 | 5 | 3.80 | 0.0243 | 0.6% |
| 301–1,000 | 14 | 7.29 | 0.0149 | 1.1% |
| >1,000 | 22 | 14.45 | 0.0056 | 0.6% |

Split the same way (`err172.py`), the whole model table lives or dies on 32 targets:

| key | `n_real` ≤ 20 (94.4% of the mass) | `n_real` > 20 |
|---|---:|---:|
| production | 0.45662 | 0.00760 |
| `\|cons\|` | 0.51751 | 0.02204 |
| M1 logistic | 0.52896 | **0.02961** |
| `\|H\|` | 0.56899 | 0.02179 |
| M4 gbt | **0.57964** | 0.02733 |
| ORACLE | 0.73155 | 0.59047 |

Note the inversion: **M1 is the *best* key on the 50 large-pool targets and only fourth on the
32 small-pool ones**, and the small-pool block decides the headline. This is 157 §3a's finding
(98% from the small-`n_real` third) reproduced independently on 164's DE table with a different
bucketing, so it is a property of the surface and not of one table.

### 7.3 At challenge power the whole channel shrinks 4×

`power172.py`, 157 §3b's proportional matching (`pool_p` = the reference's top-`round(f·n_real_p)`
by its own p-value). `f = 1` reproduces the native BH surface exactly for production, `|H|` and
the oracle — the known-answer row. *(157's own f > 1 rows are 0.01–0.02 away from these; both
tables agree exactly at f = 1, and they diverge above it because which extra genes enter a
scaled pool depends on the reference's p-value ordering below the BH cut, which is a property
of the individual DE table — 157 built its own, this uses 164's. Cross-check, not assertion.)*

| key | f=0.5 | **f=1.0** | f=1.5 | **f=2.0** | f=3.5 | f=5.0 |
|---|---:|---:|---:|---:|---:|---:|
| median `n_real` | 44 | **89** | 134 | **178** | 312 | 445 |
| production | 0.22927 | 0.18283 | 0.13165 | 0.10953 | 0.05061 | 0.01944 |
| `\|cons\|` | 0.25580 | 0.21539 | 0.17332 | 0.15226 | 0.06665 | 0.04407 |
| M1 logistic | 0.27539 | 0.22448 | 0.16261 | 0.15117 | 0.05942 | 0.04718 |
| `\|H\|` | 0.26128 | 0.23533 | 0.18290 | **0.17072** | 0.06544 | 0.04210 |
| M4 gbt | 0.28518 | 0.24287 | 0.19369 | 0.16876 | 0.06950 | 0.04432 |
| ORACLE | 0.66068 | 0.64553 | 0.62658 | 0.62047 | 0.58398 | 0.57847 |

Capture % of the oracle headroom over the production order:

| key | f=0.5 | f=1.0 | f=1.5 | **f=2.0** | **f=3.5** | f=5.0 |
|---|---:|---:|---:|---:|---:|---:|
| `\|H\|` | 7.4 | 11.3 | 10.4 | **12.0** | **2.8** | 4.1 |
| M4 gbt | 13.0 | 13.0 | 12.5 | 11.6 | 3.5 | 4.5 |
| M1 logistic | 10.7 | 9.0 | 6.3 | 8.1 | 1.7 | 5.0 |
| `\|cons\|` | 6.1 | 7.0 | 8.4 | 8.4 | 3.0 | 4.4 |

**The challenge's `n_real` sits between f = 2 and f = 3.5.** `|H|`'s raw gain over production
holds up to f = 2 (+0.0525 at f = 1, +0.0612 at f = 2) and then **collapses to +0.0148 at
f = 3.5**, with the capture percentage going 11.3% -> 12.0% -> **2.8%**. The learned keys track it down and never separate.
This independently reproduces 157 §3b's central claim — *higher power makes the ordering
channel worse, not better* — with a new key family it never screened.

---

## 8. All six members (AGENTS.md §11.1 rule 12)

**This experiment builds nothing and emits nothing.** No prediction was rendered, no counts
were emitted, no call set was recomputed. So five of six are pinned by construction, not by
argument:

| member | Δ from this experiment | why |
|---|---|---|
| `pds` | **exactly 0** | reads the pseudobulk profile; no profile was produced |
| `mse` | **exactly 0** | reads the pseudobulk profile; no counts were emitted |
| `nmae` | **exactly 0** | reads our lfc over the reference's DE genes; the lfc is untouched |
| `fid` | **exactly 0** | reads our own call set; no call set was recomputed |
| `jac` | **exactly 0** | reads our own call set; no call set was recomputed |
| **`reach`** | **counterfactual only**, table in §4 | the scorer reads our Wilcoxon p-values, not a key |

**The exemption is narrow and must not be inherited.** The moment a key is *delivered* — which
is 171's expression-compensating spread schedule — the emission law changes, and rule 12 then
binds on the composition: `mse` must be priced through 096's `claim` route with its
known-answer rows, `fid`/`jac` move because who clears BH moves, and only `pds`/`nmae` stay
bit-identical via the exact-total repair. **None of that is priced here and none of it should
be inferred from this document.** The §4 numbers are an upper bound on what the composition can
deliver, before any realisation loss.

Applying the composed slot gate (rule 13) to the best number in §4: M4's +0.00753 raw over
`|H|` scales to **+0.0085, i.e. +0.0014 overall** — an order of magnitude below the +0.0089
amplitude-changing threshold and below even the +0.0019 shape-only one. **The learned key is
not shippable on its own terms even at its most favourable reading.**

---

## 9. GO / NO-GO

**NO-GO on learning the sign-confidence key from this feature family.** The specific claims:

1. **`|H|` is not beaten.** The LODO-selected class reads −0.0109 (P(better) = 0.074); the best
   class on h1-test reads +0.0075 (P = 0.778). The range spanned by all five learned/blended
   keys, 0.2212 – 0.2429, contains 0.23533 in its interior.
2. **Learning is real but redundant.** +0.009 to +0.027 over the un-whitened consensus is
   solid (P = 0.967, 0.997) — and `|H|` gets +0.020 of that for free, with no model. Whitening
   had already extracted what the donor scalars carry.
3. **Capacity is not the limit; features are.** 22 → 406 → 2,900 parameters move mean fold
   `reach` by 0.0002 (M1 vs M4) and the additive model is *worse* than the linear one. Extra
   capacity buys front precision (macro top-1% 0.826 → 0.887) that does not convert.
4. **The features that carry the fold signal are the standardised ones** — consensus divided by
   the gene's own across-panel variability, and the panel's marginal responsiveness — not
   magnitude, not agreement, not expression. Both are already partly inside `|H|`.
5. **The ordering gap that remains is a sign-quality gap.** `oracle − |H|` = 0.410, and per
   157 §3b's identity `oracle reach ≈ accuracy/0.9`, closing it needs the per-gene sign
   accuracy to move, not the order. At 0.5343 pooled we are 0.0343 above a coin flip, and `|H|`'s own front decile is 0.621.

**What should happen next, in priority order:**

1. **171 should deliver `|H|` as its `rank_score`**, not a learned key. It is the best measured
   key, has no transfer risk, and is already what INBOX-4 recommended. `key172.npz` carries it
   as `key_absH` alongside the learned alternatives.
2. **Re-price the target before spending more on ordering.** §7.3 says the challenge sits where
   capture falls to 2.8–12%, and §8 says the best learned increment is +0.0014 overall. The
   board gap on `reach` is +0.0162. **Ordering cannot close it; sign accuracy can.**
3. **The one live lead is §6.1**: the un-whitened consensus signs the pool 0.0106 better than
   the emitted whitened ones. Recovering the un-whitened sign while keeping the whitened
   profile would be worth ~+0.012 raw `reach` at the oracle end. It touches the profile, so it
   belongs to whoever owns `pds`, and it is handed over rather than recommended.

---

## 10. The key file

| | |
|---|---|
| path | **`experiments/172-sign-confidence-key/key172.npz`** |
| sha256 | **`3bee1878143fbabf7ad7ff4f7ebe6b210fff4b29a12a52a66f39c3bdc5791492`** |
| axes | `genes` (10,779, h1-test DE axis), `targets` (85) — both aligned to `de164.a15.npz` |
| convention | higher = more confident; every array is (85, 10,779) float32 |

| array | `reach` | |
|---|---:|---|
| `key_absH` | **0.23533** | **the recommendation** — `\|H\|` from 169's recovered effect |
| `key_M4_gbt` | 0.24287 | best on h1-test, NOT LODO-selected, inside `\|H\|`'s noise |
| `key_M2_step` | 0.23163 | AgingLab-scale additive model |
| `key_blend_M1_absH` | 0.22906 | 50/50 rank blend, untuned |
| `key_M1_lodo_selected` | 0.22448 | the pre-registered deliverable |
| `key_cons` | 0.21539 | un-whitened consensus, for reference |

`meta` is a JSON string with the donor lists, the banned donor and the reach values, so the
file is self-describing. `deliver172.py` recomputes every `reach` from the saved arrays before
writing, so the file cannot ship a key whose reported number does not belong to it.

---

## 11. What was NOT measured, and why

* **`fid`, `jac`, `mse`, `pds`, `nmae` for a delivered key.** Nothing was built (§8). Pricing
  the composition is 171's, and it is not optional there.
* **A whitened feature in the LODO folds.** No leak-free per-fold whitening operator exists;
  building one means refitting `union_pruned` without the held-out line on that line's own gene
  axis, which is a different experiment. §6.1 measures what the omission costs.
* **Recipient control-cell detection rate and dispersion.** 161 measured that these two make a
  cross-line learned model worse; not rebuilt (§2).
* **The unseen D/E/F panel.** By construction untestable. The design constraint — no gene
  identity, no target identity, 21 measurable scalars, one code path for training and
  deployment — is the only assurance available, and it is the same assurance AgingLab's
  published head has.
* **The brief's top-X% curve.** Does not reproduce and contradicts INBOX-4 §1; §0.1 records the
  conflict rather than resolving it.
* **Whether M1 or M4 is "right".** They are 0.00024 apart in LODO selection and 0.0184 apart on
  h1-test (P = 0.971 paired). The selection is a coin flip and is reported as one.
* **More than 3,000 targets per fold, or more than 3 seeds per null.** 783,000 rows already
  saturate a 22-parameter model; the binding constraint is the 82-target evaluation surface,
  not the training set.

---

## 12. Files

| file | what it produces |
|---|---|
| `harness172.py` | the instrument: `H1Test`, `reach`, `deciles`, `topx`, `check()` — the known-answer gate every other script calls |
| `donors172.py` | streaming donor access, cached per-gene panel statistics, the donor↔h1-test overlap audit (§1) |
| `feats172.py` | the 21-feature map; one `build()` for folds and for deployment (§2) |
| `lodo172.py` | the three LODO folds and the DE-pool proxy validation (§1.1) |
| `fit172.py` | model classes, the step basis, fold `reach`/`front@5`, top-heavy weightings |
| `run172.py` | §A–F: LODO grid, h1-test table, both nulls, ablation |
| `anatomy172.py` | §7.2, where `reach` physically comes from |
| `err172.py` | §4.1 and §7.1, jackknife and paired bootstrap |
| `power172.py` | §7.3, 157 §3b proportional power matching |
| `deliver172.py` | §10, packages `key172.npz` and re-verifies every saved key |
| `work/inbox_check.py` | §4.2, reproduces `INBOX-calibration.md`'s four numbers |
| `work/coefs.py` | §6, the learned coefficients |
| `work/*.log`, `work/*.json` | every number above, verbatim |

No submission was built, no `cell-eval2` was run, no file outside this directory was modified.
