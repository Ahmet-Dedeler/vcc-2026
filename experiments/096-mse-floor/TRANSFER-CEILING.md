# 096 — what actually transfers between cell lines, and what the mse baseline already knows

Two measurements, both on real data, that together bound what any donor-transfer model can
score on `expr_mse_unbiased_capped_norm`. Neither was in the repo before.

## 1. Cross-donor transferability: the on-target knockdown transfers, the response does not

`cross_donor2.py`, K562-GWPS vs Orion-HEK293T, **9,237 shared perturbations, 7,976 shared
genes**, correlating the *perturbation-specific* part (`delta - mean`, because the mse baseline
already is the mean response — see §2).

| what | Pearson rho |
| --- | ---: |
| **the target gene's own knockdown** | **+0.647** (97.5% negative in both donors) |
| top-25 responders (union of the two donors' top-25 by \|effect\|) | +0.185 |
| top-100 responders | +0.105 |
| top-400 responders | +0.064 |
| magnitude-weighted over all genes | +0.118 |
| unweighted over all 7,976 genes | +0.035 |
| does "this perturbation is strong" transfer? (correlation of \|\|s\|\|) | +0.297 |

**The first version of this measurement reported only the last row (rho = 0.035) and would
have concluded that donor transfer is hopeless.** That reading is wrong, and the tell was
available immediately: we hold the best `pds_cosine` on the entire leaderboard, which is
impossible if nothing transfers. Correlating across ~8,000 genes when a few dozen respond
measures dilution, not biology (see `docs/RESEARCH-METHOD.md` rule 3).

**The finding is a split, not a null.** CRISPRi knocking down its own target reproduces
across cell lines at `rho = 0.65`. The downstream transcriptional programme it triggers
reproduces at `rho ~ 0.10-0.19`. Cell-line context governs the response; the knockdown itself
is context-free.

This is exactly the signature of our board profile: **best-on-board `pds`** (each perturbation
is identifiable by its own gene dropping — that part transfers) and **worst-quartile `mse`**
(the diffuse response, which is most of the expression distance, does not).

## 2. The mse baseline already explains 81.5% of the true signal

`ontarget_denoised.py`, on `orion-hct116-thinned` (132 real perturbations). The first pass
computed shares of `||d_p - mean||^2` from finite samples; the metric explicitly subtracts a
jackknife estimate of exactly that sampling inflation, so the naive shares are the wrong
denominator. Redone with a **split-half cross-product estimator** — cells split in two, `d1`
and `d2` formed independently, `d1*d2` unbiased for the true squared signal because the two
halves' noise is independent.

The correction is not cosmetic:

| quantity | naive | de-noised | real fraction |
| --- | ---: | ---: | ---: |
| total signal `sum \|\|d\|\|^2` | 2461.9 | 1082.9 | 44.0% |
| baseline residual `sum \|\|d - mean\|\|^2` | 1558.4 | **200.8** | **12.9%** |

**Only 12.9% of the apparent baseline residual is real.** The rest is sampling noise, and the
metric removes it.

| | naive | **de-noised** |
| --- | ---: | ---: |
| share of the TRUE signal the mean response explains | 36.7% | **81.5%** |
| on-target gene's share of the baseline residual | 1.09% | **8.39%** |

### What this bounds

- **The mean perturbation response explains 81.5% of what a perturbation actually does.**
  Only 18.5% of the true signal is perturbation-specific. This is why `b = 0.993` is a hard
  baseline, and why 376 of 385 teams score exactly 0 on this member.
- **An oracle predicting the mean response plus the exact on-target knockdown, and nothing
  else, would score `mse_raw ~ 0.910`** — comfortably below the baseline, worth about
  `+0.085` scaled (`+0.014` overall). The on-target gene is 8.4% of the residual, a hundred
  times its share of the gene axis.
- Combining the two measurements: we transfer on-target at `rho = 0.65` (capturing
  `rho^2 = 0.42` of its 8.4%) and downstream at `rho ~ 0.1-0.19` (capturing 1-3.6% of the
  remaining 91.6%). That puts a **realistic donor-transfer ceiling at `mse_raw ~ 0.94`**,
  i.e. `+0.054` scaled, **`+0.009` overall**.
- **The leader's 0.827 captures 16.7% of the residual.** On-target can supply at most 8.4%,
  so they are getting ~9% from the downstream response — implying `rho ~ 0.31` there,
  against the `0.10-0.19` that cross-donor transfer delivers. **Their advantage is a better
  model of the response, not better calibration.**

## 3. Correction to `docs/BOARD-ANALYSIS-2026-08-28.md`

That document said `mse` is 50% of our gap to first place and implied it was addressable.
The first half is right; the second is not, and the distinction matters for where we spend
the next month.

**The addressable part is about `+0.009` overall, not `+0.028`.** The rest of the mse gap is
the leader having a genuinely better model of the perturbation-specific response — the exact
quantity that cross-donor transfer caps at `rho ~ 0.19`. No amount of amplitude tuning
reaches it.

That is still worth having: `+0.009` is larger than every emission-side gain this project has
booked, including a8's `+0.0041`. But it should be pursued knowing it is a bounded prize.

## 4. The strategic consequence

The residual is **diffuse** — de-noised, the top-10 genes carry 7.7%, top-200 carry 32.7%,
top-1000 carry 70.5%. There is no small set of genes to get right. Combined with §1, the
implication for model design is:

**Predict the mean response, add the on-target knockdown at its measured transfer strength,
and shrink the downstream response toward zero in proportion to how little it transfers.**
That is a per-component shrinkage, not the single global amplitude the current pipeline
applies — and §1 says the correct shrinkage differs by a factor of ~3.4 between the on-target
coordinate (`rho` 0.65) and everything else (`rho` ~0.19).

**The current build applies one amplitude to both.** That is the specific, testable defect
this experiment identifies.

---

## 5. RETRACTION of §2's on-target conclusion — the board refutes it

§2 concluded that "an oracle predicting the mean response plus the exact on-target knockdown,
and nothing else, would score `mse_raw ~ 0.910`", from a split-half de-noised estimate that the
on-target gene is 8.39% of the baseline residual.

**Two leaderboard entries measure this directly, and it is wrong.** Mining the board's
`description` field (131 of 385 teams populate it — a source this project had never read):

| team | self-described method | `mse_raw` |
| --- | --- | ---: |
| **Kun Lab** | *"Zero-effect floor: each perturbation is 400 control cells resampled with replacement from its own context. **No perturbation signal by construction**; used to calibrate the 0-point of the metric scale."* | **1.0424** |
| **Lumos** | *"Freeze NTC cells; scale 271 GWPS-covered target genes by fold_expr (clip 0-1); 29 genes NTC-only. **On-target only, no off-panel trans.**"* | **1.0425** |

**On-target-only and predict-nothing score identically, to four decimal places.** The
on-target knockdown is worth `+0.0001` on this member, not the `0.08` §2 projected — a factor
of about 800.

### Why the estimate failed

The split-half estimator correctly de-noised the *reference*, but the metric's own correction
is not the same operation. `N_p` subtracts `C_p` and `rho*min(C_hat_p, C_p)` from the
numerator while the denominator is the measured distance from control — so the numerator is
de-noised and the denominator is not. Taking shares of a fully de-noised residual, as §2 did,
asks a different question from the one the scorer asks, and the answer differs by two orders
of magnitude. **The error was reasoning about a ratio by de-noising both halves of it when the
scorer de-noises only one.**

Kun Lab's 1.0424 is the number that anchors this member: **predicting no change at all scores
`raw ~ 1.042`**, and the scored baseline `b ~ 0.993` sits only slightly below it. Everything
between 1.042 and the leader's 0.827 is the whole playable range.

### What survives from §1 and §2

The **cross-donor transfer measurements in §1 stand** — they are direct correlations between
two donor matrices and involve no metric modelling. On-target transfers at `rho = 0.647`;
downstream at `rho ~ 0.10-0.19`. What does *not* survive is the inference from those
correlations to an mse score, because that inference went through the broken ratio.

**The lesson, for `docs/RESEARCH-METHOD.md` rule 6:** knowing a metric is clamped is not the
same as knowing what its denominator is. Before estimating a member's headroom, find a board
entry that has *measured* the relevant anchor. Two teams published theirs in a field we were
not reading.
