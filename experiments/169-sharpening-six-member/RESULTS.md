# 169 — magnitude sharpening, priced on all six members. **NO-GO**, and `nmae` is why

**Date:** 2026-09-02 · measurement only. No `cell-eval2`, no build, no `vcc prep`, no
submission, nothing written to `submissions/`, no file outside this directory modified.
Peak RSS 7.3 GB.

---

## Decision, in one line

**NO-GO on the sharpening arm, at every exponent, in both transform spaces.** The lever is as
large as 030 and 168 said — **+0.019452 scaled `pds`** on the instrument rule 13 actually
names, positive on all four folds, **2.9× the largest value in this project's history and
1.54× 162's own "nobody has ever cleared this" budget bar.** Five of the six members are now
measured and clean, including the `reach` toll that 168 could not price and that this
experiment **falsifies** (measured `+0.00061` overall against a charged `−0.00243`). It fails
on the sixth: `nmae` prices at **−0.00023** overall on HCT116 and **−0.01289** on h1-test — a
56× disagreement that neither gate-depth nor amplitude matching closes — and the surface that
says no is the one the board has validated on that member. The nearest historical arm, **a13**
(the other per-gene magnitude reweighting), was gated positive on both harnesses and lost
**−0.0203** on the board with raw `nmae` at 1.0185. Same family, same member, same shape.

**a15 stays exactly as it is. Zero flags change.**

---

## 0. What 168 left open, and what happened to each

| 168's open item | status |
| --- | --- |
| **`reach`/`fid`/`jac` unmeasured; rule 13 charges a blanket `−0.0146` scaled `reach` toll** | **CLOSED, and the toll is falsified.** Five MATCHED seeds on h1-test: Δ`reach` = **+0.00818 ± 0.00373** raw (t = 2.19) = **+0.000605 overall**, 95 % CI `[−0.000161, +0.001372]`. The toll is `−0.002428` overall; **the entire interval lies above it** and the point estimate has the opposite sign. §4 |
| **`nmae` has ONE surface** | **CLOSED as a measurement, and the answer is the experiment's result.** The second surface disagrees with the first by 56×. §5–§6 |
| **`mse` proxied by `sse*`** | **CLOSED.** 096's `claim` route, four known-answer gates bit-for-bit. `claim` **invariant** (Δ = −0.00001 against 5.3 % of headroom to 164's ρ cliff), `rA = 0.99903`, `mse` = **+0.000037** overall at the null. §3 |
| **the cliff's mechanism, and is `q` derivable?** | **CLOSED as an exact identity.** The optimum is set by a threshold whose value is **1/2**, not by a fit. §5 |
| **refine the q grid in [1.0, 1.3]** | Done; the knee moves to 1.25–1.30 once the renormalisation is fixed. §1 |
| **does sharpening double-count a23's mass fix?** | **No.** The mass fix's own `pds` gain does not erode with `q` on either surface. §1 |
| *(not asked, and it decides the result)* | 098 §4's **gate-depth artefact on `nmae`** — a rule-9-class problem I was initially blind to. §6.1 |

---

## 1. The arm, and one correction to 168's specification

    h_g  ->  sign(h_g) * |h_g|^q       on the CPM-gated coordinates, row renormalised

168 (following 030) renormalised to the panel's **median** row norm, which leaves the
across-target norm spread free to move — measured, the delivered row-norm CV goes
`0.153 → 0.180` at q = 1.15. `pds_cosine` is invariant to per-row scaling, so the two forms are
**identical on `pds`** and 168's `pds` table is a statement about both. Everywhere else the
per-row form strictly dominates:

| q | Δ`pds` worst surface | Δ`nmae` HCT116, **median** (168) | Δ`nmae` HCT116, **row** |
| ---: | ---: | ---: | ---: |
| 1.15 | +0.01145 | +0.0036 | **+0.0028** |
| 1.25 | +0.01799 | +0.0327 | **+0.0144** |
| 1.30 | +0.02082 | +0.0491 | **+0.0403** |

Same `pds`, `nmae` 2.3× cheaper at q = 1.25, and — the reason it matters for `mse` — the row
form leaves `A` untouched: **`rA = 0.99990`** (Orion pair) and **`0.99903`** (challenge panel,
panel coordinates zeroed). Everything below uses the **row** form.

**Known-answer gate (`sharpen169.py`, aborts on failure):** the median form with the mass fix
reproduces 168's `work/sharpen.json` at q = 1.15 / 1.25 / 1.50 on both surfaces and on exact
`nmae` — **9 of 9 assertions PASS**.

**Refined grid, row renorm + mass fix, worst of the two Orion lines:**

| q | 1.05 | 1.10 | **1.15** | 1.20 | 1.25 | 1.30 | 1.50 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Δ`pds` | +0.00405 | +0.00740 | **+0.01145** | +0.01440 | +0.01799 | +0.02082 | +0.02920 |
| Δ`nmae` HCT116 | −0.0001 | +0.0008 | **+0.0028** | +0.0094 | +0.0144 | +0.0403 | +0.1864 |

The knee is real and it sits between 1.25 and 1.30, not between 1.15 and 1.25 as the median
form suggested.

### Composition with a23's mass fix: they do **not** double-count

If sharpening were re-doing the mass fix's work, the mass fix's own gain would shrink as `q`
rises. Δ`pds` attributable to `VCC_MASS_FIX=minnorm`:

| q | 1.00 | 1.10 | 1.15 | 1.25 | 1.30 | 1.50 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| HCT116 | +0.00035 | −0.00023 | −0.00052 | +0.00093 | +0.00145 | +0.00197 |
| HEK293T | +0.00278 | +0.00254 | +0.00318 | +0.00353 | +0.00330 | +0.00376 |

No systematic erosion on either surface (board a23 realised `+0.001737` scaled, between the
two). They reweight different structure and compose.

### The monotone-ordering claim, as a measurement

A power transform preserves the `|effect|` ordering, so every `|effect|`-keyed policy should be
untouched. On the challenge panel, where `emission_spread_vector` keys on `|effect|` exactly,
the top-8,000 keep-block overlap with the incumbent is **1.00000 at every q**. On h1-test,
where 157's board-validated convention keys on `|log2(pcpm/cellcpm)|` — a *different* key — the
overlap is 0.984. Both reported: the first is the claim, the second is the proxy's own slack.

---

## 2. `pds` on the instrument rule 13 names — the largest signal this project has produced

Rule 13's gate is written in terms of **PLOLO** (162: four donor folds over the real 300-target
panel, rendered into the challenge control frames, 6/6 board signs, LOO RMS 0.001494 scaled).
The Orion-pair number is *not* that number — `pds` is a retrieval statistic and the decoy pool
there is 132 targets against PLOLO's 269–294. So 162's folds were re-run with the transform
inserted where production would put it.

**Two known-answer gates, both exact:**

| gate | expected | measured |
| --- | ---: | ---: |
| q = 1.00 reproduces 162's own a23 rows, 4 folds × 4 frames | — | **worst \|Δpds\| = 0.00e+00**, 16/16 bit-identical |
| a15 → a23, the mass fix's own PLOLO price | +0.004014 scaled (162) | **+0.004014** |

**The arm, challenge frame, 162's pre-committed predictor (unweighted 4-fold mean):**

| q | k562_gwps | orion_hct116 | orion_hek293t | kolf2 | mean raw | fold sd | **scaled** |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1.05 | +0.003787 | +0.003727 | +0.003672 | +0.000893 | +0.003020 | 0.001228 | +0.006680 |
| 1.10 | +0.007421 | +0.006940 | +0.007139 | +0.002236 | +0.005934 | 0.002142 | +0.013126 |
| **1.15** | **+0.010639** | **+0.010408** | **+0.010680** | **+0.003448** | **+0.008794** | 0.003088 | **+0.019452** |
| 1.20 | +0.013705 | +0.013582 | +0.013880 | +0.004731 | +0.011475 | 0.003895 | +0.025382 |
| 1.25 | +0.016881 | +0.016434 | +0.017196 | +0.006060 | +0.014143 | 0.004674 | +0.031284 |

Positive on **all four folds at every q**. For scale: the largest PLOLO-predicted Δ`pds` in this
project's history is a15's **+0.006656** scaled, and 162 Gate 4 computed that a shape-only arm
needs **> +0.0126** to clear its own unpriced noise — "nearly 2× the largest value in our entire
history … under this gate none of the six historical geometry arms earns a slot." **q = 1.15
predicts +0.019452: 1.54× that bar, 2.9× a15. This is the first arm this project has produced
that clears Gate 4 at all.**

Same run, corroborating: PLOLO's weighted sign credit `sgn_w` rises **+0.00148** and `sse*`
improves **−0.000092** at q = 1.15.

### 2.1 Adjudicating the two inboxes — the `pds` gain is real, and it is not an Orion blob

Two messages arrived in this directory mid-run and both bear on §2.

**`INBOX-from-168.md` UPDATE 3** closes the family as a "TERMINAL NO-GO" because h1-test read
GLOBAL q = 1.15 at `+0.00112` against Orion's `+0.01145`, and lays down "never size an arm on
Orion; use h1-test". **The rule is right; the replacement instrument is not.** The number that
must size a `pds` arm is `pds_cosine` in the *scorer's* space — the delivered log-pseudobulk
effect over all 18,533 genes minus the 300 panel coordinates, ranked against the panel's own
decoys. `debench.Surface.pds()` on an h1-test DE table is a **different statistic**: log2 fold
change, a CPM-5 gate, and **85** decoys against the challenge's 300. It has never been validated
against a board `pds` point. PLOLO has — 6/6 signs, LOO RMS 0.001494 — and PLOLO reads this arm
at **+0.019452** scaled, positive on all four folds. I reproduce the h1-test DE-table figure
(+0.00443 for this arm) and report it in §7's inputs, but it is not the quantity rule 13's gate
is written in.

**`INBOX-from-161.md`** asks whether 161's `+0.0206` mean raw `pds` for `|C|²` across five lines
is the same Orion overstatement, and proposes the discriminator: each surface's
target-independent energy share (published: h1-test 3.6 %, orion-hct116 26.4 %,
orion-hek293t 44.1 %). Run on PLOLO's own four folds (`blob169.py`):

| fold | targets | scored genes | blob share | Δ`pds` q = 1.15 | scaled |
| --- | ---: | ---: | ---: | ---: | ---: |
| k562_gwps | 269 | 7,448 | **0.6 %** | +0.010639 | +0.023534 |
| orion_hct116 | 290 | 17,806 | **0.3 %** | +0.010408 | +0.023023 |
| orion_hek293t | 294 | 17,806 | **0.3 %** | +0.010680 | +0.023624 |
| kolf2 (leak-free) | 273 | 17,806 | **0.0 %** | +0.003448 | +0.007626 |

**PLOLO's truths carry essentially no target-independent component** — 0.0–0.6 %, an order of
magnitude below h1-test's 3.6 % and two below the Orion *harnesses*' 26.4 / 44.1 %, because
PLOLO's truth is the held donor library's own panel effect matrix rather than a harness
pseudobulk with its depth structure. **The gain is positive on all four, including two folds
that are not Orion at all and one that is leak-free. It is therefore not blob alignment**, and
161's `|C|²` result is corroboration rather than a re-measurement of the same artefact. (The
`+0.854` correlation between blob share and Δ`pds` over these four folds is carried entirely by
`kolf2`, n = 4, and `kolf2` is more plausibly smaller because it is the fold using a
non-production donor with the shipped operator — 162 flags it as the leak-free one. I would not
read a trend into four points, three of which are at 0.3–0.6 %.)

**So the family's `pds` term is real, larger than either inbox concluded, and survives the
discriminator 161 proposed.** It dies on `nmae` (§6), not on `pds` — and that is a materially
different closure from 168d's, because it means a *fold-change-space* variant of the same lever
is still live (§6.4, §11).

168 UPDATE 3 also asked, as the highest-value remaining item, for **a measured `reach` toll for
a profile-reshaping arm**, since rule 13's `−0.0146` is calibrated on the emission arms a21+a22
only. §4.2 is that measurement.

---

## 3. `mse` — priced on the real challenge panel, and it is neutral

162 §3 measured that a cross-line validator cannot price `mse` (LOO RMS 0.112 against a
predict-zero null of 0.044 — 2.5× **worse** than predicting nothing). The only instrument that
has retrodicted a board `mse` move is 096's decomposition on `data/challenge/context_{A,B,C}`:

    mse − 1 = ( s²A − 2sX + Kc − min(s²B, claim) ) / ΣD

**Known-answer gates (`claim169.py`), all passing:**

| gate | expected | measured |
| --- | ---: | ---: |
| flat 0.25/φ0 | 0.11860 (096) | **0.11860** (−0.00 %) |
| flat 1.00/φ0 | 0.19216 (096) | **0.19438** (+1.16 %) |
| flat 0.25/φ4 = a8 | 0.16478 (096) | **0.16652** (+1.06 %) |
| a15 as shipped | 0.24581 (164) | **0.24581** (+0.00 %) |
| `A` / `B` round trip at s = 1 | 0.283614 / 0.281236 | **exact, +0.000 %** |
| a15 raw `mse` / `score_mse` | 0.88228 / 0.11112 (164) | **0.88228 / 0.11112** |

**The arm at a15's `VCC_ALPHA = 0.982`:**

| q | `A` | `rA` | `B` | `claim` | Δ`claim` | `Ĉ_p` | ρ branch |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1.00 | 0.28274 | 1.00000 | 0.28106 | 0.24621 | — | 16.65 | `s²B` |
| 1.15 | 0.27671 | 0.97869 | 0.27509 | 0.24620 | **−0.00001** | 16.64 | `s²B` |
| 1.30 | 0.26658 | 0.94284 | 0.26505 | 0.24620 | −0.00001 | 16.64 | `s²B` |

* **`claim` is invariant to 1 part in 20,000.** 164 established a15 at 99.7 % of the dispersion
  cap with `s²B = 0.23268` against `claim = 0.24581` — 5.3 % of headroom, beyond which the raw
  `mse` slope is `−4.375` with no amplitude rescue. **Sharpening does not spend it.** The
  a19/164 failure mode is excluded by measurement, not by an exemption claim.
* The 2.1 % drop in `A` is **not the arm**: `‖h‖²` is preserved to 1.00000 by construction, and
  the drop is the ON-TARGET coordinate absorbing more of the fixed budget (its share goes
  118.8 → 233.8 as q goes 1.00 → 1.15). It is an artefact of pricing on **v11's** delivered
  effect, which has no mass fix and leaks render residue onto the zeroed panel coordinates.
  Zero those 300 coordinates, as a15's mass fix does, and **`rA = 0.99903`**.
* What enters (\*) is `A − B`, the common mode, and it moves −1.1 % relative.

**`X` is the one term needing the withheld truth**, and 156/162 warn against pricing a profile
arm by scaling a harness cosine — so it is priced at three settings and the gate uses the null:

| q | rX hct116 | rX hek293t | rX PLOLO 4-fold | Δ`mse` @ rX = 1 | @ worst measured rX |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1.15 | 1.10994 | 1.17702 | **1.10349** | **+0.000037** | +0.002760 |
| 1.25 | 1.18702 | 1.29775 | 1.17357 | +0.000074 | +0.004641 |

Three independent surfaces — two held Orion lines and PLOLO's four folds in the challenge frame
— agree the alignment rises 10–18 % at q = 1.15. **The gate uses `+0.000037`;** the upside is a
transferred harness cosine and rule 13 does not credit it for a shape-only arm.

---

## 4. `reach`, `fid`, `jac` — measured, and rule 13's blanket toll is falsified

Priced on **h1-test** (157 §1: it reproduces both board raw `reach` points to within 0.014,
while repaired HCT116 misses a15's by 2.6×). Donors `k562_gwps` + `orion_hek293t`;
**`data/donors/h1_2025` covers 85/85 of h1-test's targets and never enters**, here or in any
ordering key.

### 4.1 Getting a sharpened profile onto h1-test, and proving the route

Every existing h1-test path (093/095/157/164) re-emits from the built prediction's per-gene
panel totals — i.e. holds the **profile** fixed, which is the axis sharpening has to move. So
the artifact's log-space effect is recovered from its emitted counts:
`h = log1p(TS·Z·p) − b`, with `p = (totals + 0.5)/Σ(totals + 0.5)` because
`emit_counts(mode="quasi")` sets `totals = floor(n·L·profile)` (`donor_transfer.py:380`).

**Checks that could have failed, and what each would read if the premise were wrong:**

| check | if wrong | result |
| --- | --- | --- |
| `Σ expm1(b) == TS` | any wrong control | **50000.0000** |
| off-gate `h` is exactly 0 for every target, so its ACROSS-TARGET sd must sit at the quantisation floor | a wrong gate leaves target-dependent effect in the mis-assigned band; a wrong `b` leaves a shared offset the size of a log-CPM difference | sd **0.00151** at `min_cpm 5.0` vs **0.01357** at 10.0, against an on-gate **0.06555** |
| the same, run as a threshold scan over {2, 5, 10} | — | **selects `--min-cpm 5.0`**, independently reproducing 164 §3.1, which established the same fact from a completely different statistic |
| common mode of the recovered effect | — | 12.1 % of median \|h\| against 096's independently measured 9.2 % of the row norm |
| median SNR on gated genes | — | **30×** |

**Known-answer rows on the realisation:** the 164 code path (`KA_a15_raw`, no recovery)
reproduces 157/164's median `n_pred` **6035** and `reach` **0.18283** — bit-identical. The
recovery round trip at q = 1.00 lands **0.00034** from it against a re-seed spread of 0.00547.

### 4.2 `reach` with an error bar, which is the only way it can be read

157 recorded a same-policy re-seed spread of 0.008–0.011 raw here, **the same size as every arm
delta ever measured on it**. So the incumbent and the arm were emitted at **five matched seeds**
and compared pairwise (identical reference table, control block, gene axis and seed; only `q`
differs).

| arm | seed-paired Δ`reach` raw | sd | sem | t | 95 % CI raw | **95 % CI overall** |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| **h-space q = 1.15** | **+0.00818** | 0.00834 | 0.00373 | 2.19 | `[−0.00218, +0.01853]` | **`[−0.000161, +0.001372]`** |
| h-space q = 1.25 | +0.01415 | 0.00421 | 0.00243 | 5.83 | `[+0.00370, +0.02460]` | `[+0.000274, +0.001821]` |
| lfc-space q = 1.15 | +0.00825 | 0.00683 | 0.00305 | 2.70 | `[−0.00022, +0.01673]` | `[−0.000016, +0.001238]` |

Perturbation-paired (85 pairs, the higher-powered estimator): +0.00818, t = 1.67, CI overall
`[−0.000114, +0.001325]`. Under 157's **proportional** power matching (`f = 2`, `f = 3.5`, which
brackets the challenge regime — `Surface(n_real_target=N)` is **not** used, because 157 §3a
measured that flattening deletes the small-`n_real` tercile carrying 98 % of our reach) the
point estimate wanders between `−0.0052` and `+0.0034` raw, i.e. also consistent with zero.

> **Rule 13 charges `−0.014569` scaled = `−0.002428` overall for any arm that changes the
> emitted profile. The measured value is `+0.000605` overall and its entire 95 % interval sits
> above the toll. The toll is falsified for this arm.**

Mechanism: sharpening *shrinks* the call set (median `n_pred` 6059 → 5827, −232) and the
survivors are the higher-`|effect|` calls, whose sign purity is higher (154: top decile 0.675
against 0.531). 168 flagged this as open with a plausible near-zero answer; it is near zero with
a slightly positive tilt.

### 4.3 `fid` and `jac`

Rule 9 is in force. `fid` is read as **precision** `mean_p(k/n_pred)` (no `n_real`; on the
challenge `n_pred ≈ 4,424` against `n_real ≈ 240`, so the yield cap is saturated and precision
*is* `fid` there), and `jac` is rebuilt on the board's own counts from the measured recall.

| arm | Δ precision | sem | Δ recall | sem | Δ`fid` overall | Δ`jac` overall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| h-space q = 1.15 | −0.00006 | 0.00004 | −0.01675 | 0.00805 | **−0.000032** | **−0.000419** |
| lfc-space q = 1.15 | +0.00000 | 0.00015 | −0.01080 | 0.00451 | **+0.000002** | **−0.000271** |

Both negligible; precision is inside its own seed noise.

---

## 5. The `nmae` cliff — the mechanism, exactly, and why `q` is derivable

168 asked for `Σ|h−r|` split by sign-correct/sign-wrong and overshoot/undershoot. The split is
an **identity**, not a fit. For reals `h, r`, `|h − r| = |h| + |r| − 2m` with
`m = min(|h|,|r|)·1[sign h == sign r]`, so over the reference's DE set:

> **`nmae_p = 1 + ( Σ_g |h_g| − 2 Σ_g m_g ) / Σ_g |r_g|`**  (I)

`cliff169.py` gates (I) against `debench.Surface.nmae()` and aborts below 1e−9 agreement: it
agrees to **0.0e+00** on HCT116 and **1.1e−16** on h1-test.

**(I) says the metric charges the L1 BUDGET we spend and refunds TWICE the sign-correct OVERLAP,
capped at `|r|` per gene. There are no other terms.** One unit of L1 moved onto a gene is worth
`1 − 2·(dm/d|h|)`, and `dm/d|h|` is exactly `1[sign correct AND |h| < |r|]`. Therefore:

> **An L1 reallocation improves `nmae` iff it lands where the sign-correct-AND-unsaturated rate
> exceeds 1/2, and worsens it below. The threshold is 1/2 and it is not fitted.**

(I) also reproduces from first principles the law `docs/CURRENT_STATE.md` records
(`NMAE = 1 − ⟨sign(lfc_ref), lfc_hat⟩`): that is (I) in the limit where nothing saturates.

Per decile of `|h|` at q = 1, and the factorisation is the whole story:

| decile | L1 share | `P(correct)` | `P(unsat)` | `e = P(correct & unsat)` |
| ---: | ---: | ---: | ---: | ---: |
| **HCT116** 1 | 0.0018 | 0.5133 | 1.0000 | 0.5133 |
| 5 | 0.0574 | 0.4911 | 1.0000 | 0.4911 |
| 9 | 0.1634 | 0.5838 | 0.9979 | 0.5828 |
| **10** | **0.4025** | **0.8287** | **0.7823** | **0.6295** |
| **h1-test** 1 | 0.0060 | 0.5164 | 1.0000 | 0.5164 |
| 5 | 0.0621 | 0.5336 | 0.9886 | 0.5276 |
| 9 | 0.1746 | 0.5675 | 0.7489 | 0.4250 |
| **10** | **0.3488** | **0.6172** | **0.5680** | **0.3540** |

**`P(sign correct)` rises with `|h|` on both surfaces** — 0.51 → 0.83 and 0.52 → 0.62,
independently confirming 154's top-decile purity result. What differs is `P(unsat)`: the head
runs out of headroom on h1-test and does not on HCT116. **That is the entire cliff**, and it
corrects 168's refutation rather than reversing it: 168 was right that the *count* of
overshooting genes barely moves (1.12 % → 3.40 %) and wrong to conclude the overshoot story was
dead — the operative variable is the **L1 mass** on saturated genes, 8.8 % → 16.7 % on HCT116
between q = 1.00 and 1.15, and 22.6 % → 29.3 % on h1-test.

**Derived vs fitted.** The frozen-membership linearisation of (I) predicts Δexcess from the
per-decile `e` alone and under-predicts the damage by 2–5× (HCT116 q = 1.15: predicts −11.7,
actual +8.5); the gap is exactly the membership term, i.e. the head crossing `|r|` as it grows.
So the *direction* and the *threshold* are closed-form; the *optimum* needs one measurement of
where the head saturates, which is a property of the amplitude, not of `q`. **A fitted `q` would
not transfer to D/E/F. `q* = where the L1-weighted head efficiency crosses 1/2` would, and that
is measurable on any panel.**

---

## 6. Why this is a NO-GO

### 6.1 The gate-depth artefact: a rule-9 problem on `nmae`, which I was initially blind to

Experiment 098 §4:

> "every `nmae` number this project has ever taken was read on a reference gate 2–20× shallower
> than the board's, and the **sign of the amplitude response flips** across that range … the
> challenge's reference calls ~2,000 genes per row; our harnesses call 89 (orion) and 1,120 (h1)."

Re-read with the depth swept (`nmaegate169.py`), Δ`nmae` at q = 1.15, raw, buildable `h` space:

| surface | BH 0.05 | BH 0.15 | BH 0.35 | BH 0.60 | BH 0.85 |
| --- | ---: | ---: | ---: | ---: | ---: |
| mean \|S_p\| (HCT116) | 90 | 247 | 767 | 2,149 | 5,411 |
| **HCT116** | +0.00261 | +0.00140 | **+0.00082** | +0.00079 | +0.00053 |
| mean \|S_p\| (h1-test) | 1,102 | 1,486 | 2,072 | 2,936 | 4,575 |
| **h1-test** | +0.08613 | +0.06390 | **+0.04158** | +0.05447 | +0.02239 |

The depth matters — it cuts HCT116's cost 3× — and it does **not** reconcile the surfaces. They
differ by 30–50× at every depth. (Independently, the h1-test **emitted** tables — the scorer's
own path, five matched seeds, sem 0.00006 — read Δ`nmae` = **+0.04640** at q = 1.15, so the
disagreement is not an artefact of the analytic evaluator either. §10 flags where the analytic
and emitted forms diverge.)

### 6.2 The two surfaces, and why neither can be waved away

| | HCT116 | h1-test | challenge |
| --- | ---: | ---: | ---: |
| raw `nmae` at q = 1 (BH 0.05) | 0.94786 | 0.95838 | **0.944251** (board, a15) |
| \|h\| ≤ \|r\| | 98.88 % | 91.01 % | not measurable |
| L1 on saturated genes | 8.80 % | 22.62 % | not measurable |
| median \|h\|/\|r\| | 0.1266 | 0.2825 | not measurable |
| **Δ`nmae` @ q = 1.15, overall units** | **−0.00023** | **−0.01289** | — |

* **For HCT116**: its `nmae` level is 4× closer to the board's and it sits in the
  low-saturation regime the challenge is usually assumed to be in.
* **For h1-test, and this is the stronger argument**: 098 §4 anchored a **board-measured
  derivative** (`d(raw nmae)/d(ln s)` = +0.0062 and +0.0097 from a10→a11 and a11→a12) against
  both harnesses and concluded, verbatim: *"h1-test at BH 0.35–0.45 reproduces the board's
  value. **orion never does, at any depth** … orion's `nmae` derivative should not be trusted at
  all — its 89-gene gate is 22× shallower than the board's."*
* **AGENTS.md rule 1 says take the worse surface.** Here the board-validated surface *is* the
  worse surface. Both considerations point the same way.

**I could not establish the challenge's own saturation regime, and one inherited number that
would have settled it does not survive checking.** `docs/CURRENT_STATE.md` states "`|h| ≤ |r|`
on 99.3 % of gated genes"; tracing it to `docs/LOG.md` and 098 (d1) shows it is a **harness**
measurement (99.36 % / 99.27 % on h1 / orion), not a challenge one — and on h1-test I measure
**91.01 %**, from the artifact 098 itself names, reproduced two independent ways (my recovery,
`Σ|h|/Σ|r| = 0.3108`; and 098's own `bulk_lognorm_means` convention, 0.3217; against 098's
reported 0.138), on a reference block that matches 098's **exactly** (58 rows, 64,930 pairs,
mean `|r|` 0.3639). **That number should not be used as a challenge-regime anchor until someone
re-derives it.** Flagged rather than laundered.

### 6.3 The board's own known-answer row for this family is a13, and both harnesses failed it

**a13 is this arm's nearest historical neighbour**: the other per-gene magnitude reweighting of
the effect matrix. It was gated positive on both harnesses, shipped, and lost **−0.0203
overall** — the largest single-arm loss this project has recorded — because raw `nmae` went to
**1.0185**, i.e. **Δ +0.074 raw**. 098 §4's post-mortem: at challenge depth its reconstruction
of a13's key still reads `−0.0005 / −0.0010` on the two harnesses. **Both harnesses were off by
0.075 raw and in the wrong direction on the one board test of an `nmae` shape arm.** 162
independently classifies `nmae` as not locally priceable (3/6 signs, p = 0.892) and excludes
a13 from every fit because it is not reconstructible — **so a13 is not in the 0.0019 threshold's
RMS either**, and it is the arm most like this one.

h1-test reads this arm at **+0.046 to +0.086 raw `nmae`** — the same order as a13's realised
board catastrophe. That is not proof and I do not claim it is. It is why a `+0.0032`-overall
expected value with a `−0.0094`-overall worse-surface reading does not earn a slot.

### 6.4 The alternative transform space — measured, better, still short

§5's identity says the damage is L1 landing on saturated genes. The builder's effect space is
`log1p(TS·p)`, which is **linear** below ~20 CPM and logarithmic above, so for the same fold
change `|h|` is ~2.4× larger on a 100-CPM gene than on a 10-CPM one: **a power transform of `h`
concentrates on high-EXPRESSION genes, not on large-FOLD-CHANGE genes.** That is exactly the
six-member split — `pds` and `mse` read the effect in `h` space, `nmae`/`reach`/`fid`/`jac` read
log2 fold changes.

So the arm was rebuilt in fold-change space (`sharpen_lfc`: `h → lfc → power → back to h`, then
the **builder-space** row norm restored — without that restore it silently bundles a 9 %
amplitude cut, `rA = 0.833`, i.e. two variables in one flag, which is the v3 mistake).

| q = 1.15 | `h` space | `lfc` space |
| --- | ---: | ---: |
| PLOLO Δ`pds` scaled | **+0.019452** | +0.017631 (91 %) |
| challenge `rA` | 0.99903 | 0.97810 |
| Δ`claim` | −0.00001 | +0.00010 |
| Δ`nmae` h1-test, emitted, 5 seeds | +0.04640 | **+0.00860** (5.4× cheaper) |
| Δ`nmae` HCT116 @ BH 0.35 | +0.00082 | **+0.00054** |
| Δ`reach` overall | +0.000605 | +0.000611 |
| Δ`fid` / Δ`jac` overall | −0.000032 / −0.000419 | +0.000002 / −0.000271 |
| six-member sum, `nmae` from h1-test | **−0.009422** | **+0.000930** |
| six-member sum, `nmae` from HCT116 | +0.003240 | +0.003170 |

**Better on every member and it still does not clear +0.0019 on the worse surface.** One caution
that belongs on the record: the h1-test DE-table `pds` (a log2-FC-space statistic over 85
decoys, *not* the scorer's `h`-space `pds`) reads **−0.0228** for this arm where it reads
**+0.0044** for the `h`-space one. PLOLO — which computes `pds` in the scorer's own space with
269–294 decoys and is the board-validated instrument — says +0.017631. I trust PLOLO, but a
sign disagreement between two `pds`-shaped statistics is exactly the kind of thing that should
be resolved before this arm is revisited.

---

## 7. The six-member table

Overall board units versus a15 + a23. `pds` from PLOLO × 1.011/6; `mse` at the rX = 1 null (what
rule 13 credits for a shape-only arm); `reach`/`fid`/`jac`/`nmae h1` from five matched seeds on
h1-test; `nmae hct` from HCT116 at the board-validated gate depth BH 0.35. `nmae` uses the
conservative span 0.60 (AGENTS.md §4 range 0.57–0.63); our own two board points solve to 0.6674,
which would shrink every `nmae` entry by 10 %.

### `h`-space arm (`sharpen`)

| q | pds | mse | reach | fid | jac | **nmae h1** | nmae hct | **SUM (h1)** | SUM (hct) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1.05 | +0.001126 | +0.000010 | — | — | — | — | −0.000033 | — | +0.001103 |
| 1.10 | +0.002212 | +0.000023 | +0.000823 | −0.000071 | −0.000235 | −0.011625 | −0.000112 | **−0.008873** | +0.002640 |
| **1.15** | **+0.003278** | +0.000037 | **+0.000605** | −0.000032 | −0.000419 | **−0.012890** | −0.000229 | **−0.009422** | +0.003240 |
| 1.20 | +0.004277 | +0.000055 | +0.000534 | −0.000106 | −0.000609 | −0.016340 | −0.000547 | −0.012190 | +0.003604 |
| 1.25 | +0.005271 | +0.000074 | +0.001047 | −0.000143 | −0.000682 | −0.016946 | −0.001026 | −0.011378 | +0.004542 |

### `lfc`-space arm (`sharpen_lfc`, builder-norm restored)

| q | pds | mse | reach | fid | jac | **nmae h1** | nmae hct | **SUM (h1)** | SUM (hct) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1.05 | +0.001056 | −0.000001 | — | — | — | — | −0.000019 | — | +0.001036 |
| 1.10 | +0.002024 | +0.000001 | — | — | — | — | −0.000072 | — | +0.001953 |
| **1.15** | **+0.002971** | +0.000006 | **+0.000611** | +0.000002 | −0.000271 | **−0.002388** | −0.000149 | **+0.000930** | +0.003170 |
| 1.25 | +0.004684 | +0.000025 | +0.001103 | −0.000202 | −0.000492 | −0.017735 | −0.000410 | −0.012617 | +0.004708 |

Rows at q ∈ {1.10, 1.20} carry `reach`/`fid`/`jac`/`nmae h1` from a **single** matched pair;
the matched-seed ladders were run at q = 1.15 (five pairs) and q = 1.25 (three), and the
seed-paired sem on `reach` is 0.00028 overall, so those two rows are indicative only. Rows at
q ∈ {1.05, 1.30} have no h1-test emission at all and their `SUM (h1)` is therefore not shown.

**The `nmae` column is the entire decision.** Every other member is measured, small, and either
neutral or favourable.

---

## 8. Rule 13's composed slot gate, run

**Gate 0 — eligibility.** Effect matrix only: no emitter change, no count law, no per-cell
library change. PLOLO is eligible. **PASS.**

**Gate 2 — amplitude exposure.** Row-wise L2 preserved exactly by construction; `rA = 0.99990`
(Orion) and `0.99903` (challenge, panel coordinates zeroed); `VCC_ALPHA` untouched.
**SHAPE-ONLY → threshold +0.0019.**

**Gate 1 — `pds` magnitude.** `|0.899 · Δpds| = 0.01749` against a 0.00246 bar — **7.1×. PASS.**

**Gate 4 — the budget.** `+0.019452` against the `+0.0126` shape-only bar — **1.54×, and the
first arm in this project's history to clear it. PASS.**

**The composite, `h`-space arm at q = 1.15:**

| accounting | Δoverall | vs +0.0019 |
| --- | ---: | --- |
| rule 13 literal (k = 1.011, blanket `reach` toll) | +0.000850 | FAIL |
| k = 0.899 (162's all-arm fitted scale), blanket toll | +0.000486 | FAIL |
| **k = 1.011, `reach` MEASURED (this experiment)** | **+0.003883** | PASS |
| k = 0.899, `reach` measured | +0.003519 | PASS |
| k = 0.899, `reach` at its 95 % lower bound | +0.002753 | PASS |
| k = 0.899, all six, `mse` at the null, **`nmae` from HCT116** | +0.002876 | PASS |
| k = 0.899, all six, `mse` at the worst measured rX, `nmae` from HCT116 | +0.005599 | PASS |
| **k = 0.899, all six, `mse` at the null, `nmae` from h1-test (board-validated)** | **−0.009785** | **FAIL** |
| **`lfc`-space arm, same accounting, `nmae` from h1-test** | **+0.000601** | **FAIL** |
| `lfc`-space arm, `nmae` from HCT116 | +0.002840 | PASS |

**Five of the six members clear comfortably. The sixth flips the answer depending on which of
two harnesses prices it, and the one the board has validated on that member says no.**

---

## 9. GO / NO-GO

### NO-GO on both arms, at every q. Do not spend a slot. `a15` stays exactly as it is:

```
VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.10 VCC_MASS_FIX=minnorm \
  --min-cpm 10.0 --overdispersion 4.0 --tail-spread 1.0 --keep-top 8000 \
  --keep-spread 0.25 --whitener-id union_pruned \
  --whitener data/donors/whitener_union_pruned.npz \
  --whitener-metadata data/donors/whitener_union_pruned.json
```

**Zero flags change. There is no recommended arm.**

### The invocation that would ship it, and the prediction I would preregister

Recorded so it is not re-derived, and so the prediction is on the record if a later runtime
resolves the `nmae` surface question. The builder patch is **`sharpen-builder.diff` — written,
not applied**, per the mandate; `git apply --check` passes and the patched module parses.
Diffed against a15's meta it is **exactly one** difference, `config.sharpen_q` — the patch adds
that field precisely so rule 10's diff-the-artifact check can see the arm, and it carries a
comment that `effects.whitening.after_all_median_norm` is computed **upstream** of
`profile_from_effect` and is therefore **not** a falsifier for this flag.

```
VCC_SHARPEN_Q=1.15 VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.10 VCC_MASS_FIX=minnorm \
  --min-cpm 10.0 --overdispersion 4.0 --tail-spread 1.0 --keep-top 8000 \
  --keep-spread 0.25 --whitener-id union_pruned \
  --whitener data/donors/whitener_union_pruned.npz \
  --whitener-metadata data/donors/whitener_union_pruned.json
```

**Preregistered overall prediction: `+0.0029` if `nmae` behaves as HCT116 reads it, `−0.0098`
if it behaves as h1-test reads it.** Per member at the HCT116 reading: `pds +0.00291`,
`reach +0.00061`, `mse +0.00004`, `fid −0.00003`, `jac −0.00042`, `nmae −0.00023`.

**Falsifiers, in the order the board can test them:**

1. **`pds` must rise by `+0.0175 ± 0.0035` scaled.** This is the prediction I would stake the
   instrument on: 6/6 board signs, LOO RMS 0.0015, and this is 2.9× the largest signal PLOLO has
   ever produced. If `pds` does not move, PLOLO is broken and 162 must be reopened.
2. **`reach` must not fall by more than `0.0015` scaled.** Measured `+0.0045` scaled, 95 % CI
   `[−0.0012, +0.0102]`. A drop near the blanket `−0.0146` voids 157 §1's h1-test calibration.
3. **`mse` must move by less than `±0.002` scaled**, and the built artifact's `claim` must land
   within 3 % of `0.24620` (`mech_chal_stats` + `jackknife_correction`). Anything larger means
   the profile change reached the dispersion, which §3 measures at 1 part in 20,000.
4. **`nmae` IS the experiment.** HCT116 predicts `−0.0038` scaled, h1-test `−0.192` scaled — **a
   single board observation separates them by 50×**, and it would settle a question that has now
   cost two experiments (098, this one) and one board arm (a13).
5. The artifact's config must differ from a15's in **exactly one field**, `sharpen_q`.

---

## 10. What I did NOT measure, and why

| not measured | why | can it move? |
| --- | --- | --- |
| **the challenge's own `nmae` saturation regime** | it needs the withheld reference; the one document stating it traces to a harness measurement I could not reproduce on the artifact it names (§6.2) | **This is the open question, and it is the whole result.** |
| `nmae` reconciliation between the ANALYTIC and EMITTED evaluators on h1-test | they differ ~2× (q = 1.15: analytic +0.0861 at BH 0.05, emitted +0.0464). The analytic form uses the **unfloored** profile while `emit_counts` floors `n·L·profile`, and this gate is dominated by low-count genes. 098's `validate.py` matched them to +0.0003 on 40 rows of a12; it does not hold here. §7 uses the **emitted** value (the scorer's own path); the analytic ladder is used only for the gate-DEPTH comparison, where the bias is common to both arms of every ratio. | Both say "large cost on h1-test"; the reconciliation would change the size, not the sign. |
| `reach`/`fid`/`jac` at q ∈ {1.05, 1.30} | matched-seed ladders were run at 1.15 and 1.25 only; single-seed readings there are inside the 0.00028-overall sem and are omitted rather than reported as if resolved | No — bounded by the sem |
| **`cell-eval2`** | hard constraint. Never run. | — |
| a **build** of any kind | hard constraint; the patch is deliberately unapplied | — |
| the `median`-renormalisation arm **on h1-test** | a real defect, caught and reported rather than dropped: `sharpen(mode="median")` on a **single row** degenerates to `mode="row"` (the median of a one-element list is that element), so `de169.q1.15_mf_med.npz` is **bit-identical** to `de169.q1.15_mf.npz` on every field. The renormalisation comparison was done correctly where the whole panel is transformed at once (§1), and that is what the decision rests on. | It would only tighten §1 |
| the `lfc`-NORM variant's full pricing | it carries a 9–17 % amplitude cut (`rA = 0.833`), i.e. two variables; measured and reported as a control (its h1-test `nmae` is **−0.00057**, better than either arm, and that gain is the amplitude cut, not the space) | Not shippable as one variable |
| per-arm `VCC_ALPHA` re-optimisation | second variable, and `rA = 0.999` means there is nothing to re-optimise | No |
| `q < 1` | 030 established the family is monotone in `q` on `pds` and `q < 1` loses there | Untested |

### Integrity

* **Eleven known-answer gates, each of which could have failed, all passing:** 9 assertions
  reproducing 168's `sharpen.json`; 16/16 bit-identical PLOLO rows against 162's `arms.json`
  plus its a23 delta to `+0.004014`; four `claim` gates plus the `A`/`B` round trip and a15's raw
  `mse` against 164; 157/164's `n_pred` 6035 and `reach` 0.18283 bit-identical; and identity (I)
  against `debench` to 0.0e+00 / 1.1e−16.
* **Three of my own controls failed and are reported as failures, not quietly dropped:** the
  `median` h1-test arm (above); my first `nmae` reading, taken at native BH 0.05 before 098 §4's
  gate-depth artefact was found; and my first `lfc`-space definition, which preserved the wrong
  row norm and smuggled in a 9 % amplitude cut.
* **Before trusting any control I asked what it would read if the thing were absent.** The
  off-gate recovery test would read a shared log-CPM-sized offset if `b` were wrong and
  target-dependent spread in the mis-assigned band if the gate were wrong; it reads 0.00151
  against an on-gate 0.06555, and it *selects* the right gate out of three candidates. The
  `reach` seed pairing would have read the arm delta as significant from one seed; it was not
  believed until five matched pairs put a CI on it.
* One inherited number is **flagged as not reproducible** rather than repeated: 098 (d1)'s
  `|lfc_hat|/|lfc_ref| = 0.138` and 99.36 % unsaturated on h1-test, against 0.311 / 0.322 and
  91.0 % measured here two independent ways on the artifact 098 names, over a reference block
  that matches 098's exactly (§6.2).

---

## 11. The one thing worth doing next

The mechanism in §5 is exact and it names the defect: **the builder's effect space is
expression-weighted, so a power transform there is not a fold-change sharpening.** `pds` and
`mse` are paid in that space and gain; `nmae`, `reach`, `fid` and `jac` are paid in fold-change
space and do not. §6.4's fold-change arm keeps 91 % of the `pds` gain and cuts the `nmae` cost
5.4× on h1-test — the right direction, and still short.

But the higher-value target is not a better exponent or a better space. It is **the surface
question**: no instrument this project holds can price an `nmae` shape arm, both harnesses
missed a13 by 0.075 raw *in the wrong direction*, and that single gap is now blocking the
largest `pds` lever on the board — the only arm ever built here that clears 162's own Gate 4.
**A harness that retrodicts a13's `nmae` unlocks this arm and every arm like it**, and it is
worth more than another exponent sweep.

---

## 12. Files

| file | what it produces |
| --- | --- |
| `sharp.py` | the two transforms — `sharpen` (builder space) and `sharpen_lfc` (fold-change space), with the renormalisation choice as an explicit argument |
| `totals169.py` / `cellcpm169.py` | streamed caches of the h1-test prediction totals and the mean-per-cell control CPM (keeps peak RSS off 164's 8.8 GB) |
| `probe_effect.py` | §4.1 — recovers the h1-test effect, proves the recovery, re-derives the artifact's `--min-cpm 5.0` |
| `sharpen169.py` | §1 — refined `q` grid, four transform variants, mass fix both ways, both Orion lines; the 168 known-answer gate |
| `plolo169.py` | §2, §6.4 — `pds` on 162's PLOLO folds, both spaces, two known-answer gates |
| `claim169.py` | §3 — `mse` via 096's `claim` route on the real challenge panel, four known-answer gates (`SPACE=lfc`, `MASSFIX=0` switch the variant) |
| `xratio169.py` | §3 — `rA`, `rX` on the two Orion lines, all variants |
| `realize169.py` | §4 — the h1-test emissions, 27 arms including the matched-seed ladders |
| `six169.py` | §4 — the four DE members per arm, rule 9 handled explicitly |
| `reach169.py` | §4.2 — the matched-seed paired `reach` estimate and its CI, all three arm families |
| `cliff169.py` | §5 — identity (I), its gate, the per-decile efficiency table |
| `regime169.py` | §6.2 — each surface's saturation regime, factorised |
| `nmaechal169.py` | §6.2 — saturation-free limit and the amplitude-matched sweep |
| `nmaegate169.py` | §6.1 — the gate-depth ladder, both spaces, both surfaces |
| `blob169.py` | §2.1 — the target-independent energy share of each PLOLO fold, adjudicating both inboxes |
| `gate169.py` / `final169.py` | §7, §8 — the six-member tables and rule 13's gate, assembled from the JSONs |
| `sharpen-builder.diff` | the builder patch. **Written, not applied.** `git apply --check` passes. |
| `de169.*.npz` | the h1-test DE tables (real + predicted `p`/`padj`/`lfc`, plus the per-gene spread vector) |
| `work/*.log`, `work/*.json` | every table above, reproducible by re-running the script named beside it |
