# 153 — DE-bench: a deterministic, seed-free instrument for the three DE members

Date: **2026-08-31**. Measurement only. No submission, no builder, no `cell-eval2`, no
challenge truth. Surface: `experiments/095-clean-harness/de_thinned.npz`
(orion-hct116-thinned, 132 perturbations x 10,261 gated genes).

## Why

`reach`, `fid` and `jac` have only ever been measured here by emitting 400 count vectors and
running the scorer. That path carries seed noise larger than the effects being chased —
095 §12 has two gamma-Poisson draws at reach 0.1675 vs 0.1535, and `realize_confirm.log` has
a same-policy re-seed pair **0.01222 apart**. But the definitions do not need the emitter:
`nmae`'s gate is reference-only, `reach` is a pure ordering statistic, `fid` and `jac` are
exact given a call set. All four are closed-form in (predicted lfc, call set, key).

## Validate the validator

Reproduces `cell-eval2` on the same prediction and surface (095 §11) to printed precision:

| member | DE-bench | cell-eval2 |
| --- | ---: | ---: |
| `reach` | 0.07702 | 0.0770 |
| `fid` | 0.54084 | 0.5408 |
| `jac` | 0.01203 | 0.0120 |
| `nmae` | 0.93092 | 0.9309 |

And independently reproduces 095 §6's ordering table: `|predicted lfc|` 0.10082 (+0.00448),
oracle 0.62498 (+0.10318). Runtime: under a second, zero seed variance.

**One bug caught by the check:** `nmae`'s context value is the **mean of per-perturbation
ratios**, not the ratio of sums (that is `mse`'s rule — BRIEF §7 vs §2). Using the wrong one
reads 0.0078 high. Also implemented: the >=10-gene gate, and non-finite predicted lfc read
as zero.

## 1. `reach` is capped by sign accuracy, not by the ordering key

Sign accuracy on the reference's DE pool, by `|predicted lfc|` percentile within perturbation:

| bucket | pairs | sign accuracy |
| --- | ---: | ---: |
| p90-100 | 1,185 | **0.7359** |
| p75-90 | 1,678 | 0.6383 |
| p50-75 | 2,810 | 0.5623 |
| p25-50 | 2,772 | 0.5325 |
| p0-25 | 2,834 | 0.5141 |

Monotone, and **it tops out at 0.736 against a 90% purity bar**. No stratification on our own
magnitude reaches the threshold `reach` requires. Consequences, measured: **30/132 = 22.7% of
perturbations have `k* = 0`** — the first gene in the order is wrong and the ranking never
recovers — and the median `k*` is **2** out of a median pool of 68.

`reach` is therefore a head-quality metric, not a ranking-quality metric. The lever is not a
better global order; it is a small, genuinely high-precision head.

## 2. Our UP predictions are anti-informative

Over the 11,309 (perturbation, gene) pairs in the reference's DE pool:

| | n | accuracy |
| --- | ---: | ---: |
| we predict DOWN | 5,794 | **0.6666** |
| we predict UP | 5,515 | **0.4742** |

Reference DOWN-rate **0.5979**; ours **0.5123**. 81.9% of perturbations have a reference
DOWN-rate above 0.5 (median 0.618, p10 0.466, p90 0.749) — so the asymmetry is a broad
property of the panel, not a few targets.

Predicting DOWN everywhere would score **0.598**, above our actual **0.5727**. Our per-gene
sign, used as a hard label, is worth *less than the marginal*. It is still informative
(P(ref down | we say down) 0.667 vs P(ref down | we say up) 0.526, a 14-point spread) — it is
just uncalibrated. The deployable form is a calibrated P(down), not a hard sign.

## 3. Asymmetric damping of the UP half — NO-GO, in three seconds

Damping the anti-informative half (`lfc <- lfc * lam` where `lfc > 0`) is deployable and reads
only our own sign. It loses monotonically on both DE members:

| lam | `nmae` | Δ | `reach` (\|lfc\| key) |
| ---: | ---: | ---: | ---: |
| 1.00 | 0.93092 | — | 0.10082 |
| 0.75 | 0.93323 | +0.00232 | 0.09797 |
| 0.50 | 0.93824 | +0.00733 | 0.09816 |
| 0.00 | 0.95003 | +0.01911 | 0.08903 |

Amplifying the informative DOWN half also loses (`mu` 1.5 -> `nmae` +0.02367).

**Mechanism, and it is the useful part:** damping helps only if the wrong-sign mass exceeds
the right-sign mass *weighted by our own magnitude*. It doesn't, even in a half whose
unweighted accuracy is 0.474 — so `|lfc_hat|` carries sign-correctness information **within**
the anti-informative half. That is why the top decile reads 0.736 while its own half reads
0.474, and it is an argument for a graded magnitude-aware confidence key over any sign-keyed
rule. Cost: one afternoon of scorer time saved.

## 4. Power matching — the regime artefact, reproduced and cured

`fid`, `jac` and `reach`'s normalisation all carry `n_real` in a denominator or union.
`n_real` is 68 here against a challenge value inferred at ~2,000 (094
`A8-BOARD-RESULT.md:46-72`), so every call-set cut looks free here and is charged there. That
artefact has cost three board predictions: v14's `fid`, v15's `fid`, a8's `jac`. Credit for
this correction: the peer session, `experiments/154`.

`Surface(..., n_real_target=N)` takes the reference's top-N by its own p-value instead of its
BH set — same biology, emulated power. Truncating our call set to the top M by predicted
`p_adj`, priced in both regimes:

| regime | M=500 | M=2000 | M=4000 |
| --- | ---: | ---: | ---: |
| native (`n_real` 68) | **+0.0394** | +0.0172 | +0.0008 |
| power-matched (`n_real` 2000) | **−0.2718** | −0.0014 | +0.0001 |

**The sign flips and the magnitude moves 7x.** The native harness says truncation is the
biggest free win available; the power-matched bench says it is catastrophic and that the
optimum is `n_pred >= n_real` — which is what the board twice delivered. The bench
retro-predicts a known live result the native harness gets backwards, which is the standing
requirement for trusting a harness here (AGENTS.md §11.1 rule 3).

**Immediate consequence for the emission work:** `nb_phi0.20` cuts `n_pred` 3,934 -> 630. In
the power-matched regime that is Δ`fid` **−0.3896**, i.e. about **−0.219 overall**. Its
harness-measured +0.02479 `fid` gain is the artefact at full strength.

## 5. `n_real` is NOT identified, and the three board constraints disagree

Enrichment of our call set over chance is stable on the harness across a 44x range of
`n_real`:

| regime | `n_real` | recall | enrichment vs random |
| --- | ---: | ---: | ---: |
| native | 68 | 0.573 | 1.459 |
| top-250 | 249 | 0.523 | **1.339** |
| top-1000 | 999 | 0.476 | 1.219 |
| top-3000 | 2,999 | 0.432 | 1.107 |

Against the board (v11: `n_pred` 5,961, `jac` 0.0321, ~9,900 gated):

| assumed `n_real` | implied recall | implied enrichment |
| ---: | ---: | ---: |
| **240** | 0.804 | **1.335** |
| 1,000 | 0.216 | 0.360 |
| 2,000 | 0.124 | **0.206** |

The board's constraint meets the harness's own enrichment curve at **`n_real` ~ 240**
(1.335 vs 1.339). For `n_real` ~ 2,000 to hold, our challenge enrichment must be **0.206** —
five times *worse than random* — while the same pipeline reads 1.15 at that `n_real` on the
harness. That is a 5.6x break in a quantity otherwise stable to +-25%.

But the truncation evidence runs the other way: the native regime predicts truncation *helps*
(+0.0172 at M=2000) and the board says it slightly hurt, while the power-matched regime
predicts ~flat, which matches. Hard bounds from board algebra alone: `jac <= n_real/n_pred`
gives **`n_real` >= 191**; precision <= 1 on v14 gives **`n_real` <= 4,716**.

**So `jac` level, `jac` change and `fid` change imply ~240, ~2,000 and >2,317 respectively.**
At least one standard assumption is wrong, most likely that harness call-set behaviour
transfers. This is worth resolving because it is the difference between `jac` being a
+0.05 opportunity and being closed: no team among 474 exceeds `jac` 0.037 against a 0.375-0.423
replicate anchor.

**The clean identifier is one board slot.** Hold the profile *exactly* fixed and change only
the call set in a known way (raise `n_pred` — the safe direction, and the one v14->v11
supports). `pds`, `mse` and `nmae` are then exactly invariant, and the `fid` and `jac`
responses identify coverage and enrichment separately.

## Status

Instrument built and validated; four analyses run. No promotion is claimed and no arm is
authorized — any build off §2/§3/§5 needs its own preregistration. §3 is a terminal NO-GO for
the sign-keyed damping family on this surface. Everything here is one surface; h1-test
confirmation is not done.

Files: `debench.py` (instrument), `q1_ordering.py`, `q2_asymmetry.py`, `q3_asym_shrink.py`,
`q4_regime.py`, `q5_nreal.py`.

## 6. Re-screening the PDS-killed families — 104 is reach-positive, and small

`q8_rescreen_104.py`. 104's frozen HCT116 artifact carries `baseline` / `agreement` / `full`
effect matrices on the **same 132 panel targets and the same 10,261 gated genes** as the
repaired harness, so the arms score against a real reference DE table with no alignment loss.

| arm | sign acc | reach | Δreach | Δoverall |
| --- | ---: | ---: | ---: | ---: |
| baseline (incumbent) | 0.5219 | 0.08889 | — | — |
| agreement-only ridge | 0.5219 | 0.09464 | +0.00575 | +0.00108 |
| full STRING/context ridge | 0.5219 | 0.09488 | +0.00599 | **+0.00113** |

Top-decile (by `|effect|`) sign accuracy on the reference DE pool: baseline **0.6675**,
agreement 0.6886, full **0.6945**.

Two readings, and the second matters more than the first.

**The PDS gate did reject a reach-positive signal.** 104 was closed on a top-20% sign-precision
drop and a permutation-null failure; on the members carrying 96% of the board gap it is
positive on both surfaces of the endpoint — reach +0.00599 raw and head purity +0.027.

**But it is worth +0.0011 overall, not +0.02.** Aggregate sign accuracy is *identical* across
all three arms (0.5219) because 104's intervention is positive confidence reweighting with no
sign flips — it moves better genes into the head without changing any gene's direction. That
is exactly the ceiling §1 predicts: a reweighting family can improve which genes sit in the
head, but the head's purity still tops out at 0.69 against a 0.90 bar.

**Do not read the gate criticism as "+0.019 was sitting behind the gate."** The gate was
mis-specified -- that stands, and 104 is the proof -- but what it rejected is worth a
thousandth. The leaderboard panel shows three teams taking +0.10 raw reach (~+0.019 overall
each) at near-zero `pds` cost; they are getting that from **different signs**, not from better
ordering or reweighting of the same ones. Conflating the two would be the same error in a new
direction. (Correction owed to the peer session, 154.)

So the gate-design criticism stands and its cash value is small. Combined with the peer's
154 (no available key isolates a >=90%-pure segment; best buildable key captures 3.4% of the
oracle headroom), the conclusion is that **`reach` cannot be bought by ordering or by
reweighting. It requires different signs, which is a profile problem.**

## 7. Blocked: 108 cannot be re-screened from its frozen artifact

108's `work/frozen.npz` holds `base` / `reconstruction` / `reweight` / `blend10` on a
300-target gallery that is **disjoint from the VCC panel by design** (0/132 overlap with the
harness; the gene axes share 6,778). Its held HCT116 truth was opened, used for the PDS gate,
and not persisted. Re-screening the regulon family on DE endpoints therefore requires
recomputing reference DE for those 300 non-panel targets from the Orion HCT116 source — a
real data job, and one that re-opens a consumed held surface for a new endpoint. Not started;
flagged as the decision it is rather than done quietly.

## 8. Whitening priced against `reach` — the missing half of the trade, and it does not rescue it

`q9_gamma.py`. Production's `WHITEN_GAMMA = 0.65` was chosen on a surface where
`CURRENT_STATE.md` records, in the same paragraph as the exchange rate, "reach was unreadable
again." Powers add in the eigenbasis, so applying `W(delta)` to 104's frozen HCT116 baseline
gives effective gamma `0.65 + delta` exactly. Gene axis reconstructed and **verified**:
`sorted(canonical index of (donor union & canonical) - panel)` = 17,942 columns, matching
`whitener_union.json`'s `n_cols`.

| eff. gamma | sign acc | head purity | reach | Δreach | Δoverall |
| ---: | ---: | ---: | ---: | ---: | ---: |
| **0.00** | **0.5610** | **0.7070** | **0.13360** | +0.04320 | **+0.00814** |
| 0.20 | 0.5491 | 0.7010 | 0.11103 | +0.02063 | +0.00389 |
| 0.35 | 0.5431 | 0.6863 | 0.10115 | +0.01075 | +0.00203 |
| 0.50 | 0.5388 | 0.6681 | 0.09898 | +0.00858 | +0.00162 |
| **0.65 (production)** | 0.5366 | 0.6646 | 0.09040 | — | — |
| 0.80 | 0.5339 | 0.6612 | 0.08602 | −0.00438 | −0.00083 |
| 1.00 | 0.5306 | 0.6448 | 0.07521 | −0.01519 | −0.00286 |

**Monotone, no interior optimum: `reach` always wants less whitening.** The effect is real and
had never been measured — whitening costs **0.0432 raw reach**, and it costs per-gene sign
accuracy (0.5610 -> 0.5366) and head purity (0.7070 -> 0.6646) on the way, which is the
mechanism `AGENTS.md` §4.2 predicts from the 32% transfer-cosine drop.

**And it does not rescue the trade.** 088 §1 measured the other side directly: `g=0.00`
scores `pds` 0.4189 against production's 0.8063, i.e. **−0.387 raw pds = −0.143 overall**,
against the +0.00814 measured here. Even the mild end is unattractive: gamma 0.65 -> 0.45 is
worth roughly +0.0013 overall on reach, against a pds cost `CURRENT_STATE.md`'s sweep already
put at break-even on h1-test and negative on orion.

**Conclusion: production's gamma is correct, and my hypothesis that it was mis-tuned because
reach was unmeasurable is refuted.** The reach cost of whitening is genuine, it is now
quantified for the first time, and it is an order of magnitude smaller than the pds it buys.

The one surviving variant is a **per-gene** blend rather than a global scalar: reach's gain
lives in genes whose sign whitening flips, while pds's gain lives in the high-variance shared
programme directions, so un-whitening only the low-variance coordinates might keep most of
both. That has never been tested -- only the global scalar has -- but on these magnitudes its
ceiling is about +0.008 overall and it would have to give back none of the pds to be worth it.

## Standing summary

Four channels now closed with mechanisms, all bounded at the same order:

| channel | result | value |
| --- | --- | ---: |
| emission ordering (`spread`, φ) | 095 + peer 154: best key gets 3.4% of oracle headroom | ~+0.001 |
| confidence reweighting | 104 re-screen: no sign changes, head purity 0.69 vs 0.90 bar | +0.0011 |
| call-set truncation | power-matched: catastrophic under H_A, closed under H_B | <= 0 |
| whitening gamma | monotone in reach, dwarfed 17x by the pds it buys | +0.0013 usable |

All four are bounded by the same fact: **none of them changes a gene's sign.** The board's own
within-team record shows three teams taking +0.10 raw reach (~+0.019 overall each) at
near-zero `pds` cost, so that gain exists and is not coming from any transform of our
prediction. It is coming from a better directional model.

## 9. The whitener's ridge `lam` — a two-surface, four-member improvement, never swept

`q10_spectral.py`, `q11_ridge.py`, `q12_confirm_hek.py`.

`whiten.py:15` defines `W = (C + lam*I)^(-gamma/2)` and **every call site uses the default
`lam = 0.0`.** With `lam = 0` the weights are `ev^(-gamma/2)` and the union operator's
eigenvalues run down to **8.0e-12**, so ~17,000 near-null directions are amplified
astronomically. §4.2 already names that cost — "whitening amplifies low-variance directions,
so `|predicted effect|` becomes a worse per-gene confidence signal afterwards" — and §8 above
measured it: 0.0432 raw reach and 0.024 of sign accuracy.

So `gamma` was doing two jobs at once: suppressing the dominant shared programmes (which buys
`pds`) and amplifying the null space (which destroys per-gene sign). `lam` separates them.
That is why no `gamma` was ever net positive on both surfaces — the second knob was never
turned.

**HCT116** (`de_thinned.npz`) and **HEK293T** (`de_hek293t_ref.npz`, built for this, reference
side only, scorer-exact via 092's `de_kernel`):

| gamma | lam | surface | pds_proxy | reach | sign | head | nmae |
| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 0.65 | 0.000 | HCT116 | 0.78990 | 0.09040 | 0.5366 | 0.6646 | 0.96895 |
| **0.80** | **0.020** | HCT116 | **+0.0009** | **+0.0102** | +0.0041 | **+0.026** | **−0.0046** |
| 0.65 | 0.000 | HEK293T | 0.72172 | 0.04122 | 0.5122 | 0.5625 | 0.98539 |
| **0.80** | **0.020** | HEK293T | **+0.0101** | **+0.0036** | +0.0028 | **+0.015** | **−0.0012** |

**On HEK293T every `lam > 0` cell tested improves all four members.** On HCT116 the pattern is
monotone in both knobs. `gamma 0.80 / lam 0.020` is non-negative on 4 members x 2 surfaces,
which is the two-surface unanimity this repo requires before a build.

Rough value: reach +0.0019 / +0.0007 overall, `nmae` +0.0013 / +0.0003 overall, `pds`
non-negative on both. **Conservatively +0.002 to +0.004 overall from one parameter.**

**Drop-in:** `Whitener.apply(..., lam=0.020)` at `gamma=0.80`; or, algebraically identical,
write a whitener npz with `evals' = evals + lam` and apply the standard gamma through the
existing `--whitener-id` path that a2-prune already shipped.

**Caveats.** `pds` here is an lfc-space proxy (0.78990 where cell-eval2 reads 0.7953),
validated for ranking by reproducing 088's gamma=0 collapse in sign and ordering — treat its
deltas as directional. Applied to 104's `baseline`, not the exact a15 artifact, which also
carries alpha 0.982, the CD4 residual and PerturbMatch. `renorm` holds the median row norm
fixed but alpha is not re-solved, and `n_pred` must be re-probed before a slot is spent
(a15's fourth falsifier).

## 10. CORRECTION — the lfc-space pds proxy was wrong in sign, and the pick changes

§9's `pds` was an lfc-space proxy. BRIEF §1 defines the scorer's effect as a **difference of
log1p-pseudobulk profiles**, `b = log(1 + TS * P/sum P)` with **TS = 5e4** — not a log2 ratio
of arithmetic-mean CPM. log1p compresses low-count genes and a log2 ratio does not, so the
proxy over-weighted exactly the low-expression coordinates the ridge acts on. **It was biased
on the axis under test**, and on one cell out of ten it had the sign wrong.

`exact_pds.py` rebuilds it from the harness h5ads to the letter (pseudobulk, TS = 5e4, log1p,
panel targets removed, ties sharing the average position). Production reads pds **0.80557** on
HCT116. Real effects cached at `real_delta_{hct116,hek293t}.npz`.

**Under exact pds, `gamma 0.80 / lam 0.020` -- §9's recommendation -- is pds-NEGATIVE on
HCT116 (−0.00434).** Retracted.

Fine grid ranked by the **worst** of the two surfaces (`q14_fine.py`), which is the criterion
that would have caught every transfer failure in this repo's history:

| gamma | lam | HCT Δpds | HCT Δreach | HEK Δpds | HEK Δreach | worst Δoverall |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **0.80** | **0.007** | +0.00133 | +0.00579 | +0.00341 | +0.00475 | **+0.00240** |
| 0.90 | 0.015 | +0.00087 | +0.00692 | +0.00278 | +0.00405 | +0.00206 |
| 1.00 | 0.020 | +0.00191 | +0.00680 | +0.00145 | +0.00408 | +0.00158 |
| 0.80 | 0.020 | **−0.00434** | +0.01021 | +0.00254 | +0.00360 | +0.00161 |

**`gamma 0.80 / lam 0.007` is the pick**: highest worst-surface value, and **no member
regresses on either surface**. There is a clean ridge — as gamma rises the best lam rises with
it — which is the signature of a mechanism rather than a lucky cell.

Expected value **≈ +0.0024 overall**, consistent across both surfaces.

**`nmae` caveat.** My renorm scales every configuration to production's median row norm;
`whiten.py`'s `renorm=True` scales to the un-whitened input's. Both are global scalars, so
`pds` and `reach` are unaffected (scale-invariant); `nmae` is not. Read the nmae deltas as
directional only, and re-measure at build time with the solved alpha.

## 11. Robustness of the `lam` arm against the reach artefact (155's warning)

Experiment 155 found that raw `reach` is **not monotone in direction quality** on this surface:
an arm with *worse* head purity (0.6093 vs 0.6646) posted nearly 2x the reach (0.1747 vs
0.0904) through a heavy tail of long >=90%-pure prefixes on large-`n_real` perturbations —
mean k* 14.8 vs 5.0 but **median k* 1 vs 3**, and 36% of perturbations at zero. A
target-independent consensus direction can buy reach that way. `q15_reach_anatomy.py` tests
whether tonight's ship candidate gains reach through that channel.

| surface | arm | reach | mean k* | med k* | k*=0 | reach small / mid / large `n_real` |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| HCT116 | production | 0.09040 | 5.02 | 3.0 | 22 | 0.1502 / 0.0695 / 0.0516 |
| HCT116 | **0.80/0.007** | 0.09619 | 5.50 | 3.0 | 23 | 0.1532 / 0.0790 / 0.0565 |
| HEK293T | production | 0.04122 | 4.14 | 2.0 | 29 | 0.0818 / 0.0312 / 0.0107 |
| HEK293T | **0.80/0.007** | 0.04598 | 4.72 | 2.0 | 32 | 0.0897 / 0.0373 / 0.0109 |

**Cleared, on the discriminator that matters.** The artefact's signature is head purity DOWN
while reach goes UP. The ship candidate is the opposite: **head purity up on both surfaces**
(HCT 0.6646 -> 0.6845, HEK 0.5625 -> 0.5697) with reach up, and `pds` up as well — a
target-independent consensus direction cannot do that, and 155's artefact arm indeed collapsed
`pds` to 0.5331. The gain is also spread across all three `n_real` terciles rather than
concentrated in the large tail, which is where the artefact lives; on HEK293T the large tercile
barely moves at all.

**Two caveats that belong on the record.**

1. **Median k* is flat on both surfaces** (3.0 and 2.0, unchanged). The entire gain is in the
   mean. `reach` is defined as the mean of `k*/n_real`, so this is what the metric rewards —
   but the typical perturbation is not improved.
2. **`k*=0` gets slightly worse**: 22 -> 23 on HCT116 and 29 -> 32 on HEK293T. These are
   deterministic, not seed noise, so 1 and 3 perturbations genuinely lose their entire reach
   contribution. Small against the aggregate gain, but real, and the direction is wrong.

The arm remains shippable. The honest description is a modest aggregate gain, spread across
strata, with the anti-artefact signature intact and a small tail regression.

## 12. CORRECTION — screened on the wrong operator; the finding survives, the optimum moves

The peer session found that a15 is built on **`union_pruned`**, not `union` — `--whitener-id`
defaults to `union` and silently builds on the wrong basis. Sections 9-11 were screened against
`experiments/095-clean-harness/out/whitener_union.npz` (max eigenvalue **8.036**). Production
ships `data/donors/whitener_union_pruned.npz` (max eigenvalue **5.895**). There are three of
these operators on disk on the same 17,942-column axis; I picked the experiment-local one.

Re-ran the 30-cell grid on the production operator (`q16_pruned.py`):

| gamma | lam | HCT Δpds | HCT Δreach | HEK Δpds | HEK Δreach | worst Δoverall |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **0.90** | **0.020** | +0.00179 | +0.00712 | +0.00359 | +0.00352 | **+0.00226** |
| 0.80 | 0.010 | +0.00202 | +0.00472 | +0.00312 | +0.00296 | +0.00197 |
| 0.90 | 0.015 | +0.00260 | +0.00659 | +0.00312 | +0.00267 | +0.00191 |
| 0.80 | 0.007 | +0.00307 | +0.00423 | +0.00231 | +0.00307 | +0.00166 |

**Survives:** the mechanism, the sign, and the shape. The ridge helps on the production operator
too, there is still a monotone region rather than an isolated cell, and `0.80 / 0.007` remains
positive on all three members on both surfaces.

**Moves:** the optimum is `gamma 0.90 / lam 0.020` at **+0.00226** worst-surface, ~36% above
`0.80 / 0.007`'s +0.00166. The optimal `lam` tracks the eigenvalue spectrum, and the pruned
operator's is smaller — which is exactly what the mechanism predicts, and is a further point in
its favour rather than against it.

**Remaining loose thread.** Whether 104's frozen `baseline` was itself built on the pruned
basis is unverified from the artifacts I hold. If it was plain `union`, the un-whitening step
recovers the raw effect through a slightly wrong operator and every delta here carries that
error. The identity check at `(0.65, 0.0)` reproduces `base` exactly by orthonormality whatever
the operator, so it cannot discriminate. The reason to think it is not fatal: the conclusion is
stable in sign and family across two genuinely different operators (8.036 vs 5.895 max
eigenvalue), which is the same robustness argument as the two-surface agreement.

**Process note worth carrying into AGENTS.md §11.1**, from the peer: before any board arm,
diff the finished build's meta against the incumbent's and count the config differences. It
caught this one — "two differences where there should have been one."

## 14. The emission artefact, priced in closed form on the REAL challenge controls

The scorer's DE is a Mann-Whitney U of our 400 cells against the reference's 18,400 real
controls, p-value from the normal approximation to the rank sum (BRIEF §3). So for **any**
emission law the entire call set is analytic — no emitting, no seeds:

```
Fbar = P(emitted > control) + 0.5 P(emitted == control)
z    = (Fbar - 0.5) * sqrt(12 n m / (N+1)),   n=400, m=18400  ->  factor 68.5
```

Clearing BH-0.05 over ~10k genes needs |z| ~ 2.8, so **|Fbar - 0.5| > 0.041 makes a gene
significant at ZERO predicted effect.** Measured on `data/challenge/context_{A,B,C}.h5ad`:

| emitter | Fbar median | \|z\| median | significant at ZERO effect | called UP |
| --- | ---: | ---: | ---: | ---: |
| minimum-variance (all cells at the mean) | 0.5833 | 5.7 | **95.3%** | 100.0% |
| **production** (quasi, spread .25, phi 4) | 0.5589 | 4.1 | **65.5%** | 97.6% |

(Contexts B and C reproduce A for the min-variance case: 95.9% and 96.8%.)

**The artefact is not a nuisance term. It is our call set.** The mean of a right-skewed count
distribution sits above its median, so placing cells at the predicted mean puts the median gene
nine standard errors from null before any biology is applied.

This explains, in one mechanism, a list of separately-recorded facts:

- `n_pred` ~4,400-5,900 against `n_real` ~300: we over-call because the emitter is not null
  under the null, not because our effects are large.
- 094's `|predicted effect|` vs realised p-order correlation of **+0.071**: the call set is
  selected by per-gene skew, which is unrelated to the model.
- The 76.5% up-bias, and §2's finding that our UP calls are **anti-informative** (0.4742): at
  zero effect the artefact calls 97.6-100% UP. Those are not wrong predictions; they are not
  predictions.
- `fid` ~0.50 and `jac` at baseline for us and most of the board.

**A resampling emitter removes it exactly.** Draw emitted values from the control's own law and
at zero effect emitted and control are exchangeable, so `Fbar = 0.5` identically. Over-dispersion
only damps: phi=4 took the artefact 95.3% -> 65.5%, and that partial move was worth **+0.019
raw reach** on the board (a8). Roughly twice as much artefact remains as phi=4 removed.

It also predicts the Aginglab v9->v11 signature that nothing else explained — reach +0.100,
**jac flat**, mse much better — because fixing the null changes *which* genes are called without
necessarily changing *how many*, while raising emitted dispersion toward the reference's, which
is where the mse jackknife credit lives.

**Caveat.** The quasi row models zero effect everywhere and predicts 65.5% against a board
`n_pred` nearer 44%, so it overstates the residual somewhat — real effects move genes both ways
and the BH threshold sits on a different p-value distribution. The direction and order of
magnitude are not in doubt; the exact residual is.

## 15. CORRECTION to §14 — integer ties shrink the residual from 65% to 19%

§14 (and the peer's independent run) modelled a **float point mass** at the predicted mean. The
emitter does not do that. `donor_transfer.py:299`: `L = int(min(median(lib_sizes), max_counts))`
— **one fixed library size for all 400 cells** — and `totals = floor(n*L*p)` allocated as
integers. So emitted CPM sits on a lattice of `1e6/L ~ 49.7 CPM` while control CPM
(`count * 1e6 / lib_i`, libraries 3,275-52,420) does not. Exact ties concentrate **at zero**,
where every zero-count emitted cell ties every zero-count control cell, and Mann-Whitney scores
ties at 0.5 — pulling Fbar toward 0.5 and returning genes to the null band.

`q20_integer_null.py`, real integer law, exact tie handling:

| context | L | emitter | Fbar med | sig at ZERO effect | up |
| --- | ---: | --- | ---: | ---: | ---: |
| A | 20,109 | integerised min-variance | 0.5538 | **74.5%** | 100.0% |
| A | 20,109 | **production quasi .25 / phi 4** | 0.5204 | **19.0%** | 87.9% |
| B | 19,946 | integerised min-variance | 0.5564 | 71.3% | 100.0% |
| B | 19,946 | production quasi .25 / phi 4 | 0.5218 | 21.6% | 86.1% |
| C | 20,034 | integerised min-variance | 0.5504 | 67.4% | 100.0% |
| C | 20,034 | production quasi .25 / phi 4 | 0.5168 | 16.0% | 80.1% |

**The mechanism in §14 stands; its size under the shipped emitter does not.** phi=4 removed
about **three quarters** of the artefact (74% -> 19%), not a third.

**Sizing, corrected.** phi=4 moved 74 -> 19 points and bought +0.019 raw reach. A resampling
emitter moves 19 -> 0. Scaled linearly that is **~+0.007 raw reach ~ +0.0013 overall**, not the
+0.0153 §14 implied. §14's "twice as much remains" is retracted: about a quarter remains.

**First coherent decomposition of the call set.** a15's `n_pred` ~4,400 of ~9,900 gated (~44%)
splits into ~19 points (~1,900 genes) artefact and ~25 points (~2,500) effect-driven. Neither
earlier estimate was consistent with the observed `n_pred`; this one is.

**And it costs the Aginglab explanation.** If only ~1,900 calls are artefact, removing them
cannot produce their +0.100 reach. The emitter is **no longer a sufficient explanation** for
their v9->v11 jump, and it is withdrawn as the leading hypothesis. What they did is open again.

**Recommendation reversed:** do not do the full resampling rebuild next. It is ~+0.001-scale —
the same bucket as everything else closed tonight — for a much larger build. The residual is
identifiable in closed form (the genes still significant at zero effect), so a targeted per-gene
spread/phi bump on exactly those ~1,900 captures most of it for a flag rather than a rebuild.

Credit: the tie diagnosis is the peer session's. It inferred it from a shape mismatch between
two runs — higher median |z| with *fewer* genes over the bar is impossible without different
low-|z| tails — before either of us knew the mechanism.

## 16. Most of the residual is removable without a rebuild — stop collapsing `lib_sizes`

`emit_counts` **already receives** per-cell `lib_sizes` and throws the variation away:
`donor_transfer.py:299`, `L = int(min(median(lib_sizes), max_counts))`. Control libraries on
context A span **3,275 to 52,420 — 16-fold**. Collapsing them is what puts every emitted cell on
one CPM lattice while real cells are spread across many, which is what concentrates ties at zero.

`q21_varlib.py`, closed form, exact ties, total-preserving allocation:

| emitter | Fbar median | sig at ZERO effect | up |
| --- | ---: | ---: | ---: |
| min-variance, fixed L | 0.5538 | 74.5% | 100.0% |
| **production quasi .25 / phi4, fixed L** | 0.5206 | **19.0%** | 87.8% |
| **quasi .25 / phi4, VARYING libraries** | **0.5047** | **7.1%** | **63.5%** |
| min-variance, VARYING libraries | 0.5484 | 63.5% | 100.0% |

**Not collapsing the library sizes cuts the residual artefact 19.0% -> 7.1%, a 63% reduction**,
and de-biases the calls from 87.8% up to 63.5% up. Dispersion and library variation compose:
either alone leaves most of it (74.5% -> 63.5% for libraries alone), both together nearly close
it. Scaled off the phi=4 benchmark (55 points bought +0.019 raw reach) the 12 points here are
worth roughly **+0.004 raw reach ~ +0.0008 overall** — modest, but it is a small change to an
existing function rather than the resampling rebuild, and the direction de-biasing may be worth
more than the count on `fid`, where §2 measured our UP calls at 0.4742.

**Caveat that must be checked before building.** Varying the emitted per-cell totals changes the
`mse` jackknife: `C_hat` is computed from per-cell totals `l_i` (BRIEF §2), and the noise credit
`rho * min(C_hat, C)` is where a15's entire `mse` member comes from. This is not a free change —
it touches the one member we are currently winning on outright. Price `C_hat` under varying
libraries before shipping it.

**Process note.** The first version of this file floored the expected count per cell
independently instead of allocating the gene total, silently destroying up to one count per cell.
It drove Fbar to 0.39 and inverted the call direction — a confident, plausible, completely wrong
result. It was caught only because it disagreed with §15 on the *identical* fixed-L
configuration. The fixed-L row is retained in the script permanently as a self-check that must
reproduce 19.0%.

## 17. `C_hat` under varying libraries — the blocker, priced

BRIEF §2's delete-one jackknife, computed exactly on context A with 400 emitted cells against
**400 real control cells** as the comparator (that is precisely what `C_hat` would be if we
emitted real cells, so it stands in for the cap's scale):

| configuration | `C_hat` | vs 400 real cells |
| --- | ---: | ---: |
| **400 REAL control cells** | 14.3 | 1.000 |
| min-variance, fixed L | 3.8 | 0.265 |
| **production quasi .25 / phi4, fixed L** | 9.4 | **0.656** |
| **quasi .25 / phi4, VARYING libraries** | 9.0 | **0.629** |
| min-variance, VARYING libraries | 2.6 | 0.183 |

**Not blocking.** Varying the libraries moves `C_hat` by **−4%** (0.656 -> 0.629) and stays well
under the comparator, so it cannot overshoot the cap. Since the credit is
`rho * min(C_hat, C_p)` and `C_hat` is below the cap, the credit falls by about the same 4% —
a small `mse` cost, not a cliff. It is a cost, not a gain: worth trading against the artefact cut
(19.0% -> 7.1%) and the up-bias cut (87.8% -> 63.5%), but it should be measured, not assumed.

**Implementation check:** 0.656 independently reproduces the repo's documented figure that our
emitted cells sit at **"0.67x"** the reference's dispersion (`COMPETITIVE-INTEL` Part 0). That
number was derived by a completely different route, which is good evidence this jackknife is
implemented correctly.

### The larger consequence, which partly rehabilitates the resampling build

Production sits at **0.656** of the dispersion real cells carry. A resampling emitter would sit
at **1.0 by construction** — it emits real-cell-like counts — which is roughly **+52% `C_hat`**.

That is a much bigger argument for resampling than the reach artefact §15 shrank to ~+0.0013,
and it lands on **`mse`**, where a15's `+0.120` scaled member and our board-6th raw 0.8828 live.
§15 downgraded resampling on reach grounds; this re-opens it on `mse` grounds, which is a
different and possibly larger case.

**It is not automatic and I cannot close it from here.** The credit is `rho * min(C_hat, C_p)`:
capped by the reference perturbation's own correction `C_p`, which is unobservable for the
challenge, and scaled by `rho`, itself bounded by the submission's across-perturbation spread `B`.
If the cap is already binding, extra `C_hat` is wasted — and AGENTS.md §4.1b measured that
emitting noise past the compensation point is *strictly harmful* at `rho < 1`. So the honest
statement is: resampling raises the quantity the credit is built from by ~50%, and whether that
converts depends on two terms we cannot measure offline.

## 18. The resampling law, measured — `C_hat` = 1.000 and the null holds exactly

The law: draw 400 real control cells, apply the predicted per-gene fold change multiplicatively,
renormalise each cell to **its own** original library size, integerise preserving that library
exactly. Fold changes taken from 104's frozen `baseline` — our own shipped effect magnitudes
(|eff| median 0.0152, p99 0.1008) — swept over a scale so the result does not rest on one assumed
log convention.

| fold-change scale | `C_hat` | vs 400 real cells | null-gene Fbar median | sig at null |
| ---: | ---: | ---: | ---: | ---: |
| 0.0 (zero effect) | 14.3 | 1.000 | 0.5005 | **0.0%** |
| 0.5 | 14.3 | 1.000 | 0.5005 | 0.0% |
| **1.0 (shipped magnitude)** | **14.3** | **1.000** | **0.5005** | **0.0%** |
| 2.0 | 14.3 | 1.001 | 0.5005 | 0.0% |

Against production: `C_hat` **0.656 -> 1.000** and spurious calls **19.0% -> 0.0%**,
*simultaneously*. No parametric setting of `spread` or `phi` can do both, because both buy
dispersion by scattering genes out of significance — which is why a15's tiered policy tops out
at 0.656.

**The fold-change step costs nothing**, and the reason is worth stating rather than hiding: our
effects are ~1.5% multiplicative at the median gene, which is tiny against cell-to-cell
biological variation, so resampling preserves the dispersion structure essentially untouched.
`C_hat` = 1.000 is close to tautological at our effect sizes — and that *is* the finding, not a
weakness in the test.

**Why the credit converts.** `COMPETITIVE-INTEL-2026-08-28.md:162` measured reference
`C_p` = 14.54 / 14.77 / 16.37 on contexts A/B/C against our `Ĉ_p` = 10.30, with **the cap firing
on 0 of 900 rows**. A non-binding cap means extra `C_hat` is credited rather than wasted, which
was the load-bearing caveat in §17. That document's estimate for closing the gap is **0.349 of
raw `mse`** — ours is 0.8828.

**Treat that magnitude as unverified.** It was estimated without re-emitting, which is precisely
the failure mode AGENTS.md rule 10 exists for and the same shape as the counterfactuals 095 §1
found evaporating. The direction and mechanism are well-founded; the size is not.

**Caveat on the null set:** null-gene Fbar is measured on genes whose own effect is < 1e-4, which
are largely coordinates the donor does not cover. The structural claim — resampling is null where
there is no effect — is confirmed; the null set is not a random sample of genes.

## 19. Stratified null, and the coverage trap that makes H_A/H_B gating

`q24_stratified.py` — both emitters carrying the identical effect-bearing profile, so like-for-like.

| | PRODUCTION quasi fixed-L | RESAMPLING |
| --- | ---: | ---: |
| overall significant | 17.9% | **1.2%** |
| concordance with our own predicted sign | 79.4% | **95.0%** |
| lowest \|eff\| decile: significant | 10.7% | **0.0%** |
| ...its concordance / up-rate | 43.4% / 92.5% | n/a |

**The null does not degrade at small-but-nonzero effects.** Fbar median is 0.4993-0.5004 in
every decile. The bottom decile is the artefact in pure form: production calls 10.7% of
near-zero-effect genes, at 43.4% concordance, 92.5% up.

**But the level is the problem.** 17.9% -> 1.2%. The production model under-calls against
reality (17.9% here vs a15's actual ~44% of gated), so correcting by that factor puts resampling
near ~3%: `n_pred` ~300 against a15's ~4,400. A **ten- to fifteen-fold cut**.

`fid = k / max(n_pred, n_real)`:

- **H_B (`n_real` ~240):** coverage stays 1, `fid` becomes pure precision, and at 95%
  concordance that is potentially a large gain.
- **H_A (`n_real` ~2,500):** coverage 300/2500 = 0.12, `fid` ~0.50 -> ~0.11, raw −0.39,
  **≈ −0.22 overall** — several times any plausible `mse` gain.

**So §18's conclusion that the `mse` case is independent of H_A/H_B was wrong.** The *credit
term* is independent; the *arm* is not. Removing the artefact removes the thing inflating
`n_pred` to satisfy the coverage term — **we have been getting `fid` coverage for free from a
bug.** a20 is gating, not a reading aid.

## 20. Variant B (resample + total repair) strictly dominates variant A

`q25_variantB.py`, at our shipped effect magnitude:

| variant | `C_hat` | vs real | sig% | pseudobulk cos | rel L2 |
| --- | ---: | ---: | ---: | ---: | ---: |
| A resample, no repair | 14.3 | 1.000 | 1.2% | 0.999650 | 0.0264 |
| **B resample + total repair** | **14.7** | **1.032** | 1.0% | **0.999999** | **0.0013** |

1. **`C_hat` survives the repair and is slightly higher** — 1.032x the real-cell value against
   production's 0.656. The repair redistributes counts *between* cells, so the cell-to-cell
   structure carrying the dispersion is untouched. `C_hat` and pseudobulk-exactness are **not**
   in conflict.
2. **The null survives** the repair: 1.0% vs 1.2%.
3. **Pseudobulk fidelity sizes what A risks**: A is rel L2 **0.0264** — a 2.6% profile error in
   the log1p-CPM space `pds` and `nmae` read. B is **0.0013**, twenty times tighter, restoring
   the invariance property (`donor_transfer.py:310-320`) that every prior emission arm rested on.

**Build B. There is no reason to build A.**

## 21. Pre-registering a20 — and a sign error caught in my own reasoning

a15 ships `--keep-top 8000 --keep-spread 0.25 --tail-spread 1.0`; a20 is `--keep-top 14000`,
which exceeds the gated count so every gene moves to spread 0.25. Both the peer and I reasoned
"lower spread -> tighter law -> more genes significant, so a20 raises `n_pred`."

`q26_a20_predict.py` models it and says the opposite: **3,451 -> 2,358, a 32% cut.**

**The mechanism is a non-monotonicity this repo has documented twice.** Lower spread raises
significance from 0 to 0.25 (q20: 74.5% -> 19.0%). Beyond that it reverses: at spread 1.0 the
scattered part is the entire gene total, so the emitted law's realised Fbar wanders far from 0.5
draw to draw *in either direction* and clears |z| > 2.8 by chance, while at spread 0.25 Fbar sits
consistently near 0.5206 (z ~ 1.4), below the bar. The high-dispersion tail **manufactures**
random significance. 095 §12: "past it the scatter is wide enough that genes start being called
significant in the **wrong direction**"; 093 §7 saw it at phi=1.0. I reasoned from the monotone
half of a curve the repo had already mapped.

**The discriminator survives, but only with the sign known:**

| | `n_pred` RISES | `n_pred` FALLS |
| --- | --- | --- |
| H_B (`n_real` ~240) | coverage saturated -> `fid` FALLS (dilution) | less dilution -> `fid` RISES |
| H_A (`n_real` ~2,500) | coverage recovers -> `fid` RISES | coverage worsens -> `fid` FALLS |

The two worlds predict opposite `fid` signs under either convention, so a20 remains a clean
probe — but read backwards it returns exactly the wrong answer, gating a build that swings
between +0.06 and −0.22.

**Resolution does not require trusting either model:** a15's falsifier 4 already probes
challenge-control `n_pred` before submission. Run it on a20's artifact and read the direction off
the artifact. Five minutes against an hours-long build decision.

**Caveat on my number:** one draw per gene is unbiased for the expected count but noisy, and the
two configurations consume independent RNG streams so they are not paired. The 32% gap is larger
than that noise should produce and the mechanism is independently documented, so the sign is
probably right — but a direct probe of the real artifact beats both of our reasonings.

## 22. No parametric law reaches `n_pred` ~1,000 — but a mixed law does

a20's measured `n_pred` (4,420 -> 3,789, −14.3%) confirms §21's sign but leaves `fid` coverage
**saturated under both H_A and H_B**, so a20 cannot separate them. The discriminating arm needs
`n_pred` ~1,000, where H_B leaves coverage at 1 and H_A drops it to 0.40 — a factor-two swing on
`fid` against roughly flat.

`q27_search.py`: 6 `keep_top` x 15 (spread, phi) keep-levels x 15 tail-levels, with per-gene
P(significant) computed once per level over 5 draws and policies composed by selection.

```
PARAMETRIC REACHABLE RANGE, whole grid:  model 1,986-9,107  ->  calibrated 2,860-13,116
target 1,000 actual = 694 model          ->  NOT REACHABLE by any parametric policy
```

**a20 at 3,789 is already near the parametric floor.** Pushing dispersion further re-inflates
`n_pred` through §21's non-monotonicity rather than reducing it.

**Calibration, carried as known-answer anchors rather than assumed:**

| | model | artifact | ratio |
| --- | ---: | ---: | ---: |
| a15 | 3,462 | 4,420 | 1.277 |
| a20 | 2,363 | 3,789 | 1.604 |
| | | mean | **1.440** |

The two ratios disagree by 13%, so any target carries +-13%. A 1,000 target lands ~880-1,150 —
still well inside (240, 2,500), so the discriminator survives the calibration uncertainty.

**The mixed law.** Per-gene P(significant) is 0.2379 parametric at a15's keep level against
**0.0120** resampled (§19). Mixing per gene — structurally what `mode="mixed"` already does:

| target actual `n_pred` | fraction parametric | genes parametric / resampled |
| ---: | ---: | --- |
| 700 | 0.164 | 1,624 / 8,305 |
| **1,000** | **0.256** | **2,546 / 7,383** |
| 1,500 | 0.411 | 4,082 / 5,847 |
| 2,000 | 0.566 | 5,619 / 4,310 |

Put the **top genes by |predicted effect|** on the parametric emitter and resample the rest:
that keeps the confident calls, silences the artefact in the tail, and should raise precision
at the same time — §1 has top-decile sign accuracy at 0.74 against 0.51 in the bottom decile.

**It is not a throwaway probe — it is the resampling arm at partial exposure.** If `fid` holds,
we are a quarter into the build we want; if it collapses, we learn H_A at 25% exposure instead
of discovering it at 100% and eating the −0.22. Either way it puts the first point on the
coverage curve between a15's 4,420 and the full arm's ~170, where we currently have none.

## 23. a19 FAILED on the board — and the failure indicts this instrument

```
member    a15         a19         delta
pds       0.745871    0.744688    -0.00118
mse       0.120057    0.024149    -0.09591   <-- the entire loss
nmae      0.094239    0.102690    +0.00845
fid      -0.030397   -0.026398    +0.00400
reach     0.115368    0.118338    +0.00297
jac      -0.000619   -0.000053    +0.00057
overall   0.174087    0.160569    -0.01352
```

Predicted **+0.002**, delivered **−0.01352** — past the −0.001 falsifier, so by the criterion
registered before the result, **the local instrument is refuted for this lever.** Rank 12 -> 26.

**Every member the screen measured moved as predicted.** reach beat prediction (+0.00297 against
+0.00109/+0.00018); nmae improved; fid and jac improved; pds missed by 0.0017. The arm lost
entirely on **`mse`, which sections 9-13 never measured on any operator at any cell.**

**The lesson, and it is mine.** A two-surface, three-member screen certified an arm on the
strength of the three members it happened to include. For `spread`/`phi` arms that would have
been defensible — the exact-total repair (`donor_transfer.py:310-320`) *guarantees* the profile
is untouched, so `pds`/`mse`/`nmae` cannot move. **A whitener change carries no such guarantee**:
changing the eigenvalue floor changes the effect vector's shape, `mse` reads the profile
directly, and a15's alpha was solved at gamma 0.65 and never re-solved at 0.80. I flagged
"alpha is not re-solved" as a caveat and then recommended the arm anyway.

**Standing rule for this file and for AGENTS.md: price every arm on all six members, or state
explicitly which are unmeasured and why they provably cannot move.**

## 24. The mixed law priced on the members that matter — and one build detail that decides it

`q28_mixed_all6.py`. Top-x genes by `|predicted effect|` on the parametric emitter, the rest
variant-B resampled.

| x | genes parametric | `C_hat` | vs real | pseudobulk rel L2 | model sig% |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.000 | 0 | 14.8 | 1.034 | 0.0015 | 1.0% |
| 0.164 | 1,628 | 14.9 | 1.040 | 0.0015 | 8.2% |
| **0.256** | **2,542** | **14.8** | **1.034** | **0.0015** | **10.1%** |
| 0.411 | 4,081 | 14.4 | 1.012 | 0.0015 | 10.4% |
| 1.000 | 9,929 | 11.5 | 0.809 | 0.0015 | 11.4% |

**`mse` is protected, and should improve.** Emission reaches `mse` through exactly two terms.
`||b_hat - b||^2` is pinned — rel L2 **0.0015**, identical to pure resampling. And `C_hat` goes
**0.656 -> 1.034**, a +58% credit increase, on a cap that fires 0/900 rows. That is the opposite
of a19, which moved the *operator* and therefore the profile's shape with alpha stale.

**THE BUILD DETAIL THAT DECIDES IT: repair gene totals AFTER mixing, not before.** The two blocks
carry different per-cell library structures (parametric at fixed L, resampled at real libraries),
so splicing genes between them breaks the pinning: **rel L2 0.0155 without the post-mix repair
against 0.0015 with it.** A 1.55% profile drift is exactly the channel that cost a19 its `mse`
member. Assert rel L2 < 0.002 on the built artifact.

**Revised `x` mapping — measured, not the additive-rate model.** §22's additive estimate
under-predicts: it gave 7.0% at x=0.256 against 10.1% measured, because mixing changes the
per-cell library structure and therefore the resampled genes' significance too. Use this table:
x=0.256 -> ~1,003 model -> **~1,445 actual**; for ~1,000 actual target **x ~ 0.14**. Both sit
inside the discriminating window (H_A coverage 0.58 and 0.40 respectively against H_B's 1.0).

## 25. `fid`'s base rate is the NON-significant mass, not the reference's DE pool

a21 (x=0.14) probed at `n_pred` **1,141** — on target, discriminator intact — but its call-set
up-rate rose 62% -> 85%, raising the concern that an 85/15 up-split against a "59.7% DOWN"
reference must damage precision.

**That reads the wrong population.** §2's 59.7%-down is over the reference's *DE set*. BRIEF §4
counts `k` over our called genes that are **adjudicable** — any gene with a defined non-zero
reference lfc, significant or not — and with `n_pred` >> `n_real` most of our calls are
reference-non-significant genes, which drift UP:

| reference lfc UP-fraction | HCT116 | HEK293T |
| --- | ---: | ---: |
| over the reference's DE set | 0.4021 | 0.4377 |
| over ALL adjudicable gated genes | 0.5327 | 0.5579 |
| **over reference NON-significant** | **0.5338** | **0.5611** |

This independently reproduces `CURRENT_STATE`'s "reference non-significant genes drift 54.8% up"
from a different route. Expected precision with no skill, `u*q + (1-u)*(1-q)`:

| up-rate `u` | HCT116 | HEK293T |
| ---: | ---: | ---: |
| 0.62 (a15) | 0.5079 | 0.5139 |
| 0.85 (a21) | **0.5229** | **0.5406** |

**Raising the up-rate helps `fid` by +0.015 to +0.027 raw**, not hurts. The no-skill assumption is
supported: a15's raw `fid` is 0.5029, on that line — which is also why `CURRENT_STATE` records
that removing the up-bias wholesale *costs* ~−0.046.

**Where it flips, and by how little.** a21's calls concentrate on high-|effect| genes, likelier to
sit in the reference's DE set (base rate 0.40 up). Under H_A at enrichment ~1.3, ~33% of the
1,141 calls are DE-set genes, blending `q` to ~0.487, where u=0.85 gives 0.491 against u=0.62's
0.497 — a **−0.006** difference. Under H_B only ~3% are, and the up-bias helps outright. Honest
range: **−0.006 to +0.027 raw**, sign set by enrichment.

**Attribution is safe because the signal is 14x larger:**

| | a15 raw `fid` | a21 predicted |
| --- | ---: | ---: |
| H_A (`n_real` 2,500) | 0.503 | ~0.49 x (1141/2500) = **0.224** |
| H_B (`n_real` 240) | 0.503 | ~0.49-0.53, roughly flat |

A 0.28 raw swing against a +-0.02 confound. **This is not the a19 shape** — there the unmodelled
member was seven times larger than the entire predicted effect; here it is a fourteenth of the
signal. Larger `x` and the inverse mask were both rejected: larger `x` trades a measurable signal
for an unnecessary confound, and the inverse mask silences exactly the genes where skill lives
(top-decile sign accuracy 0.74 against 0.51 in the tail).

**Independent confirmation on the built artifact.** a21's up-fraction by confidence decile reads
`0.935, 0.306, 0.150, 0.137, ...` with the discontinuity exactly at the `param_mask` boundary —
direct evidence on a real package that the resampled tail went quiet as §18-19 predicted.

## 26. There is no mse/reach frontier — 15 teams hold both, and we are not one of them

Q30 asked whether Aginglab's reach jump was headroom we still have. The board answers something
sharper. Live snapshot, 514 teams:

**Teams with raw `mse` < 1.00 (actually collecting the credit) AND raw `reach` > 0.20: fifteen.**

| team | avg | raw mse | raw reach | raw nmae | pds |
| --- | ---: | ---: | ---: | ---: | ---: |
| Aginglab.com | 0.21122 | 0.8270 | 0.2621 | 0.8892 | 0.8199 |
| Jurassic Park | 0.19145 | 0.8453 | 0.2324 | 0.9252 | 0.8217 |
| We527 | 0.18942 | 0.8454 | 0.2132 | 0.9186 | 0.8186 |
| Haystack | 0.17987 | 0.8310 | 0.2443 | 0.9366 | 0.7989 |
| PumpkinChicken | 0.17187 | 0.8751 | 0.2351 | 0.9334 | 0.7884 |
| ... 10 more | | | | | |
| **Ahmet Dedeler** | **0.17409** | **0.8738** | **0.1812** | **0.9443** | **0.8360** |

**Among the top 30 by overall: 18 collect `mse`, 26 have reach > 0.20, and 14 have BOTH.**

So the trade we have been reasoning about is not a frontier. Holding `mse` credit *and* reach
above 0.20 is **normal** among the teams ahead of us. Our reach deficit is specific to our
pipeline, not imposed by our `mse` position.

**Q30's supporting observation, which set this up:** among within-team pairs with reach gains
> +0.04, the ones starting at raw `mse` ~1.0-1.2 were not "already tuned" — 1.0 is the *baseline*,
i.e. scoring zero on the member. Every one of those gains came with `mse` getting **worse**. The
board's normal move is to buy reach by over-dispersing past the point where `mse` is collectable.
That is the trade *most* teams make. **The fifteen above did not have to make it.**

**And we are the most pds-heavy team in both groups.** Our `pds` 0.8360 is the highest in the
"collects mse" set and above every member of the "both" club (0.788-0.849), while our reach is
near the bottom and our `nmae` (0.9443) is among the worst. We have over-invested in the one
member we already lead.

**The prize, restated with this constraint removed:** taking Aginglab's reach (0.2621) and `nmae`
(0.8892) while holding our own `mse` and `pds` is **+0.0152 + 0.0149 = ~+0.030 overall -> ~0.204**,
near first. Fifteen teams demonstrate that combination is reachable at our `mse` level. Nothing
about our position forces the trade we have been assuming.

## 27. The transfer ceiling, and why the remaining gap is target-specific or nothing

Four measurements, in sequence. All on HCT116's DE pool, HEK293T as the donor, truth on both
sides — no estimation noise anywhere.

**(a) The cross-context ceiling.** HEK293T's *measured* effect predicts HCT116's per-gene sign at
**0.7051** (0.7974 on the donor's strongest quintile). Our shipped prediction: **0.5724**.

**(b) Almost all of it is NOT target-specific.** Predicting HCT116 target X from HEK293T target
**Y != X** scores **0.6836** (20 permutations, sd 0.0099). So:

| | sign accuracy |
| --- | ---: |
| chance | 0.500 |
| shared, perturbation-independent programme | **0.684** |
| + target-specific transfer | 0.705 |
| **target-specific excess** | **+0.0215** |

**Cross-context transfer carries ~2 points of target-specific sign information.** That is the
ceiling on the entire donor-transfer paradigm from one donor, and ~40 experiments have been
competing for it.

**(c) The shared programme is enormous and we deliberately remove it.** Its mean over targets
scores **0.7880** alone — far above our shipped 0.5724. Whitening and row-centring exist to
strip it (112: centring removes the translation exactly), because it is non-discriminative.
The donor tables are already centred, which is why adding *their* mean back does nothing
(§q36: 0.5724 -> 0.5661).

**(d) Adding it back: the largest DE movement this project has recorded, priced out by `pds`.**

| beta | sign acc | head purity | reach | nmae | pds |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.00 | 0.5724 | 0.7359 | 0.10082 | 0.93092 | 0.74613 |
| 0.10 | 0.6171 | 0.7873 | 0.12368 | 0.90797 | 0.68292 |
| 0.35 | 0.7025 | **0.8447** | 0.19982 | 0.85362 | 0.60381 |
| 0.60 | 0.7480 | **0.8759** | 0.30650 | 0.81206 | 0.56980 |
| 1.00 | 0.7830 | 0.8641 | **0.37599** | 0.78636 | 0.54569 |

Head purity reaches **0.876** against the 0.90 bar `reach` needs — nothing else in this project
has come close. `reach` nearly quadruples. `nmae` improves by 0.145 raw.

**And every beta is net negative.** At 0.35: reach +0.0186, `nmae` +0.0215, `pds` **−0.0525**
overall. The ratio is ~2.5:1 against and constant in beta — there is no interior optimum. This is
129's "adding broad direction makes each query more similar to decoys faster than to its own
truth", now measured on both sides at once instead of argued.

### What this means for the remaining 0.040

**The shared channel is exhausted by construction** — it is worth 18 points of sign accuracy and
`pds` charges 2.5x for every one of them.

**The target-specific channel is worth +0.0215 of sign accuracy, total, at zero estimation
noise.** That is the ceiling every donor-transfer experiment has been competing for, and it
explains the ~40 NO-GOs better than any individual post-mortem did: they were not failing to
find the signal, they were splitting two points of it.

**So Aginglab's advantage cannot be shared direction** — that costs `pds` and their `pds` is
0.8199 against our 0.8360, a 1.9% gap that buys only ~+0.003 on this exchange rate. Their reach
0.2621 and `nmae` 0.8892 must be **target-specific**, from a source or mechanism that is not
one-donor cross-context transfer.

**The strategic consequence:** more donors, better blending, better whitening and better
reweighting are all competing for +0.0215. Closing 0.040 requires target-specific per-gene
direction from somewhere the donor-transfer paradigm does not reach.

### 27b. The shared-direction intervention does NOT transfer — refuted on the second surface

§27(d) found a clean interior optimum on HCT116: at `beta = 1.0`, reach `0.10082 -> 0.37599`,
`nmae` `0.93092 -> 0.78636`, `pds` `0.74613 -> 0.54569`, **net +0.01806 overall**. The relative
profile error also *fell* (0.9988 -> 0.9587), so it was not a metric trick — the prediction got
genuinely closer to truth.

Second surface, donors swapped so neither line uses its own mean (HEK293T evaluated with
HCT116's mean):

| beta | head purity | reach | nmae | pds | NET overall |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.00 | 0.5633 | 0.04044 | 0.98578 | 0.74746 | — |
| 0.60 | 0.7033 | 0.04305 | 0.97288 | 0.58004 | **−0.05766** |
| 1.00 | 0.7587 | 0.03670 | 0.96434 | 0.55291 | **−0.06648** |
| 2.00 | 0.8229 | 0.05148 | 0.94404 | 0.52759 | −0.06740 |

**Negative at every beta, by 0.046 to 0.069.** Head purity rises just as much (0.563 -> 0.823)
and `nmae` improves — but **`reach` barely moves** (0.040 -> 0.051, against HCT116's 0.10 -> 0.40)
while `pds` craters identically. The gain converts on one surface and not the other.

**Terminal NO-GO for the intervention.** The HCT116 result was surface-specific, and shipping on
it would have been a second a19 — a large single-surface gain on members a screen happened to
measure. Two-surface worst-case is the criterion and it says no, unambiguously.

**What survives is the measurement, and it is the durable part:**

- cross-context transfer's **target-specific** sign information is worth **+0.0215**, total, at
  zero estimation noise;
- the **shared** component is worth **+0.18** but is priced out by `pds` at ~2.5:1 *and* does not
  transfer;
- our shipped prediction sits at 0.5724 against a 0.7051 one-donor ceiling, of which 0.6836 is
  shared and unusable.

So the donor-transfer paradigm has roughly **two points of per-gene sign headroom**, and ~40
experiments have been splitting it. Closing 0.040 overall requires target-specific direction from
outside that paradigm.
