# Where the remaining score is, and how confident I am about each route

> 🚫 **Slots buy information, never rank.** Never re-submit an already-scored package and never
> submit to restore a rank. The objective this document exists to serve is **first place on
> 2026-11-05**, not the daily leaderboard row — rank is informational until the Oct 22 test
> release. See `AGENTS.md` §6.

> **START AT ADDENDUM 13** (2026-08-30) if you want the current position. This document is
> chronological and the body below is the 2026-08-21 framing, most of which has been
> superseded — `mse` "structurally unreachable" (§3) is retracted, the amplitude direction
> in Addendum 5 §3 is retracted, and Addendum 7's "the ordering is steerable" is retracted
> by Addendum 8. The STATE route is closed by Addendum 9. Each retraction is recorded in
> place rather than edited away, because the
> reasoning that produced the wrong call is usually the useful part.
>
> The short version: a14 remains the best submitted artifact. Calibration and
> emission tuning are spent; the current gap is mostly per-gene direction.
> Experiment 080 is now a terminal NO-GO. New work must preserve the
> same-target incumbent and add independently validated residual signal.
> See `docs/CURRENT_STATE.md`.

Written 2026-08-21, after v1 (0.1000, rank 3) and v3 (0.0825, rank 4). Leaders: Vivai 0.1341,
IsoCell 0.1207. This is a working hypothesis document, not a plan of record — several sections
are guesses I have labelled as such.

---

## 1. First, a framing correction: there is nothing to reverse-engineer

It is tempting to think of the scorer as a black box to be probed with submissions. It is not.

- **`cell-eval2` is the scorer.** Same code Arc runs, MIT licensed, installed locally.
- **The metric spec is published** — 708 lines of derivations, vendored in `docs/reference/`.
- **The anchors are computable.** `cell-eval2 baseline` + `prep-real-bundle` produce the `b` and
  `r` ends of every metric's scale from any reference dataset.
- Their `-r4` anchor revision is one version past what is published, but that only affects how we
  read *their* numbers; our harnesses compute their own anchors and are internally consistent.

**The one thing we do not have is ground truth for contexts A, B and C.** Everything else is open.
So the entire question is: *how do we build the best possible stand-in for those three lines?*
Submissions are not a probe to be spent on exploration — they are a scarce confirmation channel
(2/day) for hypotheses that already survived local screening.

---

## 2. What our validation stack can and cannot see

| harness | what it scores | status | validated? |
|---|---|---|---|
| **LODO / PDS** (exp 010) | `pds` only | working, 11 s per config | **Yes** — retro-predicted the v1→v3 regression at −0.0104 vs the board's −0.0089 |
| **h1-test count-level** | all six metrics | anchor splits running | not yet |
| h1-fast count-level | all six | working | **No** — lied about `jac` (0.127 predicted, 0.0035 actual) |
| H1 transfer cosine | nothing directly | working | **No** — +20% cosine produced −0.019 `pds` |

Two failure modes burned us, and they are different:

1. **Wrong biology.** H1 hESC is embryonic; A/B/C are almost certainly somatic. Donor weights
   tuned against H1 did not transfer.
2. **Wrong objective.** Mean cosine and 300-way PDS ranking genuinely disagree — averaging donors
   smooths predictions toward each other, raising mean alignment while *reducing* discrimination.
   Measured: equal-weight donor averaging is **worse than a single donor** on PDS (0.7264 vs
   0.7246 gene_norm) while being clearly better on cosine.

LODO fixes both. But it scores only `pds`, and `nmae` was 87% of v3's loss. **Neither tier alone
would have caught v3.** Both are needed.

---

## 3. Where the remaining 0.034 to #1 actually is

| metric | gap to Vivai (scaled) | as overall | my read |
|---|---|---|---|
| `pds` | +0.1044 | **+0.0174** | Hard. Donor tuning nearly exhausted. |
| `nmae` | +0.0491 | **+0.0082** | **Underexplored.** Never tuned against a trustworthy harness. |
| `reach` | +0.0480 | **+0.0080** | **Underexplored.** Never touched at all. |
| `jac` | +0.0041 | +0.0007 | Nearly closed. |
| `fid` | −0.0014 | −0.0002 | Already ahead. |
| `mse` | 0 | 0 | **Structurally unreachable** — see below. |

**`mse` is closed for us, and the reason is exactly quantified.** `min_α mse = 1 − cos²θ`, and
beating the 0.989 baseline needs cos θ > 0.105. Our measured transfer cosine is ~0.08.

> **Corrected 2026-08-23 — see Addendum 4.** The paragraph above was true when written and is now
> false in three places. `mse` is **not** closed "permanently, for everyone": Ibrahim Mansour
> has held raw `mse` below 1 since 2026-08-22T00:20 and now sits at `0.9095`, which forces his
> cos θ ≥ 0.3008. It does **not** need directional accuracy to "roughly double" — experiment 077
> measured our floor at `0.99297` against a `0.989` threshold, so the requirement is
> **cos θ: 0.0838 → 0.1049, a factor of 1.25**. What remains true, and is confirmed twice over by
> exp063 and exp077, is that *no rescaling of our current effect directions can score it* — we are
> already sitting on the floor, and the floor is `0.004` too high.

**So half the gap is `pds` and half is `nmae` + `reach`.** I have spent almost all effort on `pds`
and almost none on the other two — partly because they need the count-level harness, which is only
now becoming available. That asymmetry is the clearest signal about where to look next.

---

## 4. Breakthrough candidates, ranked

Each is scored by my honest guess at *probability it helps* × *size if it does*. These are
guesses; treat the ordering as a hypothesis about where to spend time, not a forecast.

### 4.1 Cross-donor agreement as per-perturbation confidence — **testing now, strongest candidate**

Experiment 011 (done): agreement between two donors predicts how well the blend transfers to a
third. Spearman **0.324** pooled (p = 2e-20), consistent across all three held-out lines. Q5
transfers at 0.0795 vs Q1 at 0.0246 — **3.2× better**.

The property that makes it usable: **agreement is computable without the target.** On A/B/C we
have no ground truth, but we can still measure whether K562, HCT116 and HEK293T agree about a
given perturbation.

Why it should pay:
- `α*_p = cos_p · ‖δ_p‖/‖d_p‖`, so the optimal per-perturbation magnitude scales with how good
  that perturbation's direction is. Agreement estimates exactly that.
- Per-perturbation scaling is **`pds`-neutral** — cosine distance is invariant to scaling one
  perturbation's own effect — so this is a free move on our best metric.
- It should directly serve `nmae` (stop overshooting on perturbations we cannot predict) and
  `mse`, which are precisely the metrics that punish confident wrongness.

This is the first idea that attacks `nmae` with a mechanism rather than a knob. **Highest priority.**

### 4.2 `reach` — completely unexplored

`reach` asks: rank the reference's DE genes by *our* confidence; how deep does directional purity
stay above 90%? Baseline 0.072, replicate 0.968 — **the widest span of any metric**, and we sit at
0.042 scaled.

We currently rank genes implicitly, by whatever p-value ordering falls out of the emitted counts.
That is almost certainly not the best available ordering. Candidates: rank by donor effect
magnitude, by cross-donor agreement *per gene*, or by predicted |lfc| × agreement. This costs one
count-level harness run to test and nobody has looked at it.

### 4.3 Identify what A/B/C actually are — *speculative, and I want a rules check*

The basal expression profile is a fingerprint. CCLE/DepMap publishes RNA-seq for ~1,400 cell lines.
Matching A/B/C against that catalogue would likely identify the tissue, possibly the line.

If it works, it tells us which donor to weight — the single question basal-similarity weighting
failed to answer in experiment 010.

**Two caveats, and the second is important:**
- Arc probably chose lines *without* public Perturb-seq, or the task would be trivial. So this
  likely yields tissue identity, not a matched donor.
- **Rules.** Allowed (Ahmet, 2026-09-26): conditioning on or blending a same-line public screen is fine.

### 4.4 More donor lines

We have 3. scPerturb + Nadig give Jurkat, HepG2, RPE1, K562-essential, plus neurons and iPSC lines.
They cover **0/300** of the VCC panel (the panel excludes common-essential genes), so they cannot
serve as donors — but they **can** add LODO folds, which would tighten every weight estimate we
make. Cheap, low ceiling, improves confidence rather than score.

### 4.5 Per-gene confidence, not just per-perturbation

Same logic as 4.1 one level down: for each (perturbation, gene), do the donors agree on the sign?
Genes where donors agree get full magnitude; genes where they disagree get shrunk. This is a
*directional* refinement, so unlike 4.1 it could move `pds` — but experiment 004 and 006 both
found that per-gene reweighting schemes lose to uniform. Prior: skeptical, worth one cheap test.

### 4.6 Things I have ruled out — do not revisit without a new mechanism

| approach | result |
|---|---|
| Low-rank / SVD truncation | loses monotonically; spectrum is flat (800 comps = 51% var) |
| Supervised per-gene transfer weights | 0.067 vs 0.095 uniform, even fitted on the same line |
| Rank-normalise, binarise, clip, keep-most-variable | all lose 30–60% of signal |
| Empirical-Bayes shrinkage | +10.6% cosine, **−0.017 overall on the board**. Wrong objective. |
| Basal-similarity donor weighting | worse than every fixed weight, monotonically |
| Equal-weight donor averaging | worse than a single donor on PDS |

The pattern across all of these: **the transferable signal is high-dimensional, spread across many
modestly-responding genes, and carried in the magnitudes.** Anything that compresses, reweights or
discretises it loses. This is also why the published deep models underperform — every bottleneck
destroys exactly what transfers.

---

## 5. Open questions I cannot currently answer

- **Is `pds` improvable at all beyond donor choice?** Everything that worked was data acquisition
  (more donors) rather than method. That may mean the ceiling is set by how much public
  perturbation data exists in lines resembling A/B/C — which would be a real, structural limit.
- **What are the leaders doing that we are not?** Vivai and IsoCell have identical raw `pds`
  (0.7835), which strongly suggests a shared direction source. Their `nmae` differs (0.1204 vs
  0.0715), so the difference between them is calibration, not model. We are at 0.7383 raw `pds` —
  behind both. Something in their direction estimate is better than ours and I do not know what.
- **How much does the H1-vs-somatic mismatch cost?** LODO transfers at ratio 0.85 to the board on
  the one pair we have. One calibration point is not a calibration curve.
- **Does the even emitter actually help on the real board?** It beat NB decisively on the fast
  harness, and our `fid`/`jac` beat Jurassic Park's published NB configuration — but we have never
  run the controlled comparison on a trustworthy harness.

---

## 6. What I would do next, in order

1. **Per-perturbation confidence from cross-donor agreement** (4.1). Mechanism, not a knob, and it
   targets `nmae` — half the remaining gap that I have never properly attacked.
2. **`reach` ordering experiment** (4.2). Widest span of any metric, completely untouched, one run
   to test.
3. **Finish the count-level harness** and re-validate the emitter and α against it. Every DE-metric
   conclusion currently rests on a harness that demonstrably lied once.
4. Ask about 4.3 before building it.

Submissions are confirmation, not exploration. Two per day is plenty if each one carries a
hypothesis that already survived both harnesses.

---

# Addendum, 2026-08-21 22:45 UTC — the picture after whitening

Everything above was written when `pds` looked nearly exhausted ("Hard. Donor tuning nearly
exhausted.") and the conclusion of the day was that *uniform, full-rank, fixed-weight donor
transfer is remarkably hard to improve on*. Twelve ideas had been tested and the three that
worked were all more data or better emission, never better inference.

That conclusion was wrong, and it was wrong in an identifiable way.

## What was actually being measured

Every rejected idea — SVD truncation, per-gene weights, EB shrinkage, basal weighting, equal
donor averaging, inverse-variance pooling, MODZ — was an attempt to build a **better estimate
of the effect**. All of them failed, and the pattern was read as "the estimate is as good as
the data allows".

The real constraint was different: **the estimate was fine and the metric's geometry was
wrong for it.** `pds` is a retrieval metric in a space where a few shared transcriptional
programmes carry most of the variance. Every prediction was being compared to every true
effect through a cosine that those programmes dominate. Flattening the space
(`experiments/017`, `019`, `020`, `024`) is worth **+0.070 macro PDS** — four times the entire
multi-donor gain, and larger than everything else in the notebook combined.

The tell was there in the graveyard the whole time. Experiment 004 found that truncating the
donor SVD lost *monotonically* — the low-variance directions carry signal. The conclusion
drawn was "stop compressing". The conclusion available was "so try amplifying them".

## The revised gap analysis

| metric | v3 (board) | after whitening (est.) | Vivai | still available |
|---|---|---|---|---|
| `pds` | 0.504 | 0.60–0.64 | 0.627 | small |
| `nmae` | −0.020 | ~0.02 | 0.120 | **+0.10 scaled ≈ +0.017 overall** |
| `reach` | 0.055 | ~0.07 | 0.090 (Bio_NUDT 0.134) | **+0.06 scaled ≈ +0.011 overall** |
| `fid` | −0.047 | ~−0.05 | −0.041 (best −0.016) | +0.03 scaled ≈ +0.005 |
| `jac` | 0.004 | ~0.005 | 0.008 | ~0 |
| `mse` | 0 | 0 | 0 | 0 until cos > 0.105 |

`pds` is close to converted. **The remaining ~0.03 is `nmae` + `reach` + `fid`, all four of
which are read off a Wilcoxon on emitted counts** — so the next binding constraint is the
emission model, not the effect estimate.

## The one idea that is now clearly next

Experiment 013 proposed a per-gene **emission variance** channel: emit trusted genes at the
variance floor so they are significant and rank first, and untrusted genes with inflated
variance so they fall down the p-value ordering — leaving the mean profile untouched, so
`pds`/`mse`/`nmae` do not move while `reach` and `fid` do.

Experiment 023 supplies the missing half of the argument. After whitening, `|predicted effect|`
is a *worse* per-gene confidence signal, while the un-whitened cross-donor geometric signal
holds up. So the two channels want **different objects**, and there is no conflict between
them:

- **mean profile** ← whitened effects (serves `pds`, `mse`, `nmae`)
- **per-gene emission variance** ← un-whitened cross-donor agreement (serves `reach`, `fid`)

That is the first design in this project where one metric can be improved through a degree of
freedom the others cannot see. It needs the count-level harness, and `h1-fast` is now known to
be usable for `fid`/`reach`/`nmae` (experiment 023) even though it lies about `jac`.

## What did not change

- `mse` is still structurally unreachable and still must not drive a shared parameter.
- The target's basal control profile has now failed **three** times — donor weighting (010),
  per-gene modulation (017), and as the parameter of a principled change of variables (026).
  It defines the gene universe and nothing else.
- The graveyard stays closed: experiment 022 re-ran it under whitening and only donor
  weighting changed rank.

---

# Addendum 2, 2026-08-22 01:40 local — two transforms, and what the metric actually rewards

The first addendum said `pds` was close to converted and the remaining ~0.03 lived in
`nmae` + `reach` + `fid`. Half of that is now wrong: a second `pds` lever turned up that is
almost as large as whitening, and the two compose.

## Where `pds` stands

| | 3-line LODO (VCC panel) | 4-line LODO (independent) |
|---|---|---|
| incumbent (4:1 donors, nothing else) | 0.7378 | 0.7732 |
| equal donors + whitening | 0.8011 | 0.8365 |
| **+ magnitude sharpening** | **0.8415** | **0.8644** |

**+0.104 macro PDS** from the incumbent, replicated on a harness with different folds, a
different predictor count and a mostly different gene population.

## What the two levers have in common, and it is the useful part

Both are **smooth reweightings that beat their hard-truncation counterparts**:

| basis | hard version | smooth version |
|---|---|---|
| directions across perturbations | SVD truncation loses monotonically (exp 004) | spectral flattening **+0.070** |
| magnitudes within a perturbation | top-k loses or barely helps (exp 030c) | power transform **+0.038** |

And neither improves *alignment*. Whitening **lowers** the transfer cosine 32% while raising
`pds` 9.5%. Sharpening makes the prediction far spikier than the truth (matching exponent
1.05, optimal exponent 2.0–2.5) and raises `pds` anyway.

**The through-line: `pds` rewards discrimination, and discrimination is not accuracy.** Every
rejected idea in this notebook — SVD truncation, per-gene supervised weights, EB shrinkage,
inverse-variance pooling, MODZ, `log1p` re-expression, noise-covariance whitening — was an
attempt to make the estimate *better*. Every one failed. The two things that worked make the
estimate measurably *worse* by the natural yardstick and reallocate its similarity so each
prediction resembles its own target more than its 299 decoys.

That is a statement about the competition, not just about our model. `mse` is 0.000 for every
team on the board and `fid` is negative for every team, so the field is being ranked almost
entirely on a retrieval metric while optimising a regression objective.

> **Corrected 2026-08-23.** No longer true of `mse`. One team broke it, and that team leads the
> board precisely because of it. See Addendum 4.

## What is still open, in order

1. **Does sharpening survive `nmae` at count level?** The `pds` gain is exact. The cost is
   estimated by a proxy that puts the `nmae`-optimal alpha at 0.75 where the real metric
   measures 1.08 — wrong on the axis that decides it. Experiment 031 is running the q = 1.0 /
   1.5 / 2.0 arms on the somatic harness. This is the only thing standing between us and a
   second large increment.
2. **The variance channel** (exp 027): +0.028 scaled `reach`, 3.4σ, mean profile provably
   untouched. Also awaiting the somatic harness.
3. **`nmae` still has no lever of its own.** It rewards alignment, and both of our `pds`
   levers cost alignment. If sharpening ships, `nmae` gets worse, and it is already our worst
   metric relative to the field (−0.020 against Vivai's +0.120). At some point the right move
   may be to accept `nmae` as the price and spend nothing on it.
4. **More somatic donor lines** — still the only thing that would raise the ceiling rather
   than the efficiency. Equal weighting winning means each new line now helps directly.

## Plan of record for 2026-08-22 (written before v6's score is known)

Two submission slots per UTC day, and **rank uses the latest one**, so the order matters more
than the count.

1. **Read v6's score first, and extract two things from it beyond the number.**
   `scripts/analyse_submission.py` does both. (a) A **third LODO→board calibration point**:
   the two we have are 0.010 apart and imply slopes between 0.33 and 0.86, and v6 is a change
   six times larger than the gap they span, so it collapses that uncertainty. (b) The **implied
   transfer cosine**, inverted from raw `mse` — the one quantity that says whether `mse` is
   ever reachable and the one no local harness can measure against the real contexts.
2. **Read experiment 031** (somatic HCT116 harness, 132 real VCC-panel targets). Its q = 1.0 /
   1.5 / 2.0 arms decide sharpening; its `spread` arm decides the variance channel. Both are
   currently blocked on an instrument, not on an idea.
3. ~~**Slot 1: v7**~~ **Cancelled — experiment 031 refuted sharpening before the window
   opened.** On somatic VCC-panel data at count level, q = 1.5 takes raw `pds` 0.7633 → 0.8068
   and raw `nmae` 0.9609 → **1.1208**; q = 2 takes `pds` to 0.8319 and `nmae` to **1.2596**.
   The six-metric average falls 0.1042 → 0.0803 → 0.0558. The experiment 018 proxy put `nmae`
   at q = 2 at 1.05 against the real 1.26 — **understated by 4x**, which inverts the decision.
   `submissions/v7.h5ad` stays unbuilt-into-a-`.vcc` and unsubmitted.
4. **Slot 2 (or slot 1): the variance channel** (`--emit quasi`, k = 1000, spread = 0.15), if
   the remaining arms of 031 replicate experiment 027's +0.028 scaled `reach` on somatic data.
   Unlike sharpening, it preserves the pseudobulk by construction, so it cannot move `nmae` —
   which is precisely the failure mode that killed sharpening.

If v6 *regresses* instead, `submissions/v4b.vcc` is the restore path and the diagnosis is that
the 2026 contexts' effect covariance differs from the donor lines' more than the held-out donor
lines differ from each other — in which case the whitening operator should be refitted on a
basis that includes the contexts' own control-cell covariance, which is the one piece of
target-side structure we have never used.

**What not to spend time on.** The magnitude-reweighting family is exhausted (exp 030f: row and
column sharpening are two routes to the same +0.04, and iterating them makes it worse). The
graveyard stays closed (exp 022). `mse` must not drive a shared parameter — but "unreachable"
is now known to be wrong by a factor of 1.25 in transfer cosine, see Addendum 4. The target's basal profile has failed three times and defines the gene universe and
nothing else.

**The one thing that would raise the ceiling rather than the efficiency is more somatic donor
lines**, and equal weighting winning means each new line now helps directly instead of being
penalised.

---

# Addendum 3, 2026-08-22 — after v6 landed at 0.1369 (rank 3/96)

## Where the score actually is now

| metric | ours | best on board | gap |
|---|---|---|---|
| `pds` | **+0.701** | **+0.701 (us)** | — |
| `mse` | 0.000 | +0.073 (Ibrahim Mansour) | +0.073 |
| `nmae` | +0.065 | +0.131 (Vivai) | +0.066 |
| `reach` | +0.090 | +0.166 | +0.076 |
| `fid` | −0.037 | −0.005 | +0.032 |
| `jac` | +0.002 | +0.015 | +0.013 |

Best-in-class on all six is **0.180**. We hold the best `pds` on the board and every remaining
point is in the five metrics we have barely touched.

## The strategic tension that now defines the problem

Whitening buys `pds` **by sacrificing cosine** — measured: −32% cosine for +9.5% PDS. `mse` and
`nmae` both reward cosine. So our best lever is actively foreclosing two of the five remaining
metrics, and the teams above us on those two (Ibrahim's `mse` +0.073 implies cos ≈ 0.20 against
our ~0.13) are the ones *not* making that trade.

The trade is still overwhelmingly correct — dropping whitening costs 0.168 scaled `pds` to buy
0.021 `nmae` and nothing on `mse` — but it means:

- **`mse` is closed for us specifically**, not for everyone. Shrinking alpha to reach it does
  not work at our cosine: alpha 0.5 gives raw `mse` 1.052 and alpha 0.35 with spread gives
  1.030, both still above the 0.989 baseline, while the board's own frontier shows the teams
  that shrink that far (larry, Vu Mai, raw `mse` ~1.003) sit at `reach` 0.090 against our 0.158.
- **`nmae` has no lever of its own** and is the largest genuinely accessible gap.

## What is validated and ready

**v8** = v6 + exact-total quasi emission at `spread 0.25`. Verified on the 21 GB artefacts:
pseudobulk **bit-identical** to v6 on 18,533/18,533 genes, so `pds`, `mse` and `nmae` cannot
move at all. `reach` +14%, `fid` +0.004, `jac` −0.004 on the cell-count-matched harness.
≈ **+0.007 overall**, downside bounded to a metric worth 0.002 to us.

Note: the `.vcc` is **3.67 GB against v6's 936 MB** — the even emitter's base/base+1 counts are
pathologically compressible and realistic ones are not. Upload will take ~17 minutes rather
than 4. Not a limit, but plan for it.

## Method note worth keeping: which harness to believe

They disagree, predictably, and the split follows how they were built:

| question | trust | why |
|---|---|---|
| `fid` yield, `reach` depth — anything through `n_real` | **h1-fast** | 400 reference cells, as the competition has |
| which genes carry signal — `pds`, `jac`, the CPM gate | **somatic HCT116** | the actual VCC panel, in a somatic line |

Orion runs at ~154 cells/perturbation, so the somatic reference is underpowered, `n_real` is
small, and a reduced `n_pred` looks free. It reported `fid` climbing to 0.595 where `h1-fast`
shows the same configuration collapsing to 0.372 on the yield cap.

## Next candidate after v8: extend the whitening operator to the full gene axis

The operator is currently fitted on genes **every** donor measures — 7,312 of 18,533. At the
CPM>5 gate the scorer tests 9,923 genes in context A, of which **6,729 are on the operator axis
and 3,131 are not** (Orion-only). So roughly **a third of the tested gene universe never gets
whitened at all**, and whitening is our single largest lever (+0.070 macro PDS).

Experiment 025 fixed the *scaling* consequence of that (`--whiten-renorm preserve`, +0.006), but
the off-axis genes are still passed through untouched rather than transformed.

**The extension.** Orion covers 18,106 of the 18,533 target genes — a superset of the
intersection — and carries 17,767 + 17,977 perturbations. Fitting the operator on **Orion's gene
axis** using both Orion lines' effects would whiten ~98% of the tested universe instead of 76%.
K562-gwps contributes its 7,679 genes to the *prediction* either way; it just would not
contribute to the *covariance* outside the intersection.

**Cost.** An 18,106 x 18,106 eigendecomposition: ~5e13 flops, roughly 10-20 minutes, with ~2.6 GB
for the Gram and the same for the eigenvectors. Feasible on this machine.

**Why it is not trivially testable in LODO.** The folds that would show the effect are the ones
holding out an Orion line, where the scored axis is Orion's full 38,325 genes — a 38k x 38k
decomposition, hours and ~11 GB. Testing it needs the axis subsampled to ~15,000 genes, which is
a different operator from the one production would use. Worth doing carefully rather than
quickly.

**Expected size.** Unknown, and the honest prior is "less than proportional": the off-axis genes
are Orion-only, so they have no cross-donor consensus behind them and are individually noisier
than the intersection genes. But a third of the universe is a large fraction to be leaving
untransformed under our biggest lever.


---

# Addendum 4, 2026-08-23

> ## ⚠ RETRACTION, 2026-08-23 (later the same day)
>
> **The transfer-cosine numbers below are unfounded. Do not use them.**
>
> They were derived by applying `min_α mse = 1 − cos²θ` to the board's
> `expr_mse_unbiased_capped_norm` column. That identity holds only when the
> metric is `||pred − true||² / ||true||²`, i.e. when predicting *nothing*
> scores exactly `1.0`. Fitting `mse(α) = c0 − 2αA + α²B` to exp077's three
> alpha points gives **`c0 = 1.3388`**, not `1.0`, so the board column uses a
> different normalization — most likely against a reference prediction rather
> than against zero. Under the fitted model the cosine would be
> `A/√(B·c0) = 0.508`, not the `0.0846` published here.
>
> Both numbers are untrustworthy. Three points fit a three-parameter quadratic
> **exactly**, so the "perfect fit" cited as validation had zero degrees of
> freedom and proved nothing, and the metric is *capped*, which can bend the
> parabola away from a true quadratic.
>
> **Specifically retracted:** our cos θ ≈ `0.0846`; Ibrahim's `cos θ ≥ 0.3008`;
> "3.6x our directional accuracy"; the `1.25x` improvement target; the
> `cos > 0.1314` break-even for low alpha.
>
> **What still stands, being direct observation rather than derivation:**
> our Orion floor is `0.99286` at `α* = 0.392` against a `0.989` threshold, so
> **we miss by `0.0039`**; exp077's low-alpha arms are negative; `pds` is
> scale-invariant in alpha (`+0.00026` across a 2.7x change), so the real trade
> is **nmae-for-mse, not pds-for-mse**; and Ibrahim's `0.9095` is the only
> sub-1.0 raw mse on the board.
>
> Resolving the normalization needs a fourth alpha point — see experiment 088.

 — after v9 landed at 0.14440 (rank 2/146)

## The single number that describes this competition

`expr_mse_unbiased_capped_norm` is `||pred − true||² / ||true||²`, so **predicting no effect at
all scores exactly 1.0**, and the baseline to beat is `0.989`. Read the board's raw column:

| team | raw `mse` | implied cos θ |
|---|---:|---|
| **Ibrahim Mansour** | **`0.9095`** | **`≥ 0.3008`** |
| us (v9) | `1.7661` | no bound |
| Cell-tics | `1.8365` | no bound |
| kriagr | `3.7001` | no bound |
| VIVAI | `5.9465` | no bound |
| pop | `6.7378` | no bound |
| Remigiusz Kinas | `10.7532` | no bound |

**Every team except one is worse than predicting nothing.** Since any α lands at or above the
`1 − cos²θ` floor, Ibrahim's `0.9095` forces `cos θ ≥ 0.3008` with no further assumption. That
bound is the most informative single fact we have about anyone else's model.

## How far we actually are

Experiment 077 swept alpha down and measured our floor directly:

| alpha | raw `mse` | scaled |
| ---: | ---: | ---: |
| `1.0796` (production) | `2.0537529` | `0` |
| `0.35` | `0.9968821` | `0` |
| `0.40` | **`0.9929697`** | `0` |

`cos θ = √(1 − 0.9929697) = 0.0838`, independently confirming the `~0.08` above. **We miss the
threshold by `0.00397`** and need `cos θ > 0.1049` — a **1.25x** improvement, not a doubling.
Ibrahim is at `3.6x` ours.

## Why this reframes everything

We are **first on the board at pds** (`pds_cosine 0.8278`) and fifth-worst on raw `mse`. Ibrahim
is *behind us on pds* (`0.8024`) and leads anyway. Our whole approach buys pds by spending
cosine — whitening is documented above as lowering transfer cosine 32% while raising pds — and
v9 shows that trade is close to exhausted.

The next real gain is not another pds lever. It is raising the **underlying transfer cosine**,
which would improve pds *and* open mse simultaneously. Every one of exp063, exp077 and the
Addendum-2 whitening analysis now points at the same variable from a different side.

## What is now known not to work

- Lowering alpha (exp 077): reaches the floor, still scores 0, costs `nmae`. Negative.
- Blending back the raw donor-consensus residual (exp 063): floors at `0.995`. Negative.
- Both fail for the same reason: they rescale a direction instead of improving it.

## Harness note: how to project a board score

Do **not** read h1-test's `avg6`; it is scaled against h1-test's own baseline. Predict *raw*
deltas, then apply the v8→v9 measured realizations — pds `1.12x`, nmae `0.64x`, reach `0.23x`,
jac `0.18x` — and then the board's raw→scaled amplification: pds `2.24x`, nmae `−2.16x`, fid
`3.46x`, reach `1.14x`, jac `2.59x`. h1-test is trustworthy on pds and overstates the DE metrics
by 4–5x. Full derivation in `docs/V9-RESULT.md`.

---

# Addendum 5, 2026-08-25 — the board is one member wide, and it is not `pds`

Addendum 4 said the open direction was raising the transfer cosine to open `mse`.
Both halves of that are now wrong, and the correction came from two places: a
deliberate board probe (v12), and solving the scoring function against our own
submissions.

## 1. The scoring function, solved

`overall = (1/6) Σ sᵢ`, `sᵢ = (uᵢ − bᵢ)/(rᵢ − bᵢ)`. Fitting `b` and `r` from our
two board points reproduces v11 (`0.15020`) and v12 (`0.13997`) to five
decimals, and every fitted anchor lands inside the spec's published range.
`nmae` and `reach` recover **identical `b` from both submissions to six
decimals**. The linear model is not an approximation of our board — it is our
board.

Consequence: **+0.01 overall = +0.06 on any single scaled member**, and the
total pot per member is fixed and knowable:

| member | our scaled | total overall points available, ever |
| --- | ---: | ---: |
| `pds` | +0.736 | **+0.062** |
| `reach` | +0.104 | +0.156 |
| `mse` | 0.000 | +0.167 |
| `nmae` | +0.093 | +0.258 |
| `fid` | −0.036 | +0.272 |
| `jac` | +0.003 | **+0.442** |

**`pds` is finished.** It is our best member and it contains +0.062 overall in
total even at raw 1.0. Whitening, sharpening, donor weighting, the union axis,
PerturbMatch, CD4 — the entire body of work in this repo — competes for a pot
smaller than one seventh of what sits untouched in `jac`.

This is the single most important number in the project and it was never
computed. Every gap analysis before this one compared us to *other teams* on
each metric, which measures what is currently being done, not what is
available.

## 2. `mse` is a solved constraint, not a live one

v11's raw `mse` of 1.5801 unclamped would score **−0.6200**. The `[0, 1]` clamp
is therefore worth **+0.1033 overall to us right now**, and above raw `0.98789`
the marginal cost of any further degradation is **exactly zero**.

So the constraint that has shaped this project since experiment 001 — "`mse`
must not drive a shared parameter", "`mse` is structurally unreachable",
"shrinkage is a losing trade" — was real but is now spent. We are already
paying the maximum possible price for `mse`, which is to say none.

Everything the DE members want from a larger amplitude is **unopposed**.

## 3. Alpha's marginal direction is up

Three points, including the exact `nmae(α=0) = 1.0` anchor that holds by
construction:

| α | raw nmae |
| ---: | ---: |
| 0 | 1.0000 |
| 0.40 | 0.9702 |
| 1.0796 | 0.9447 |

Slopes `−0.0745` then `−0.0375`: convex and **still descending at production
alpha**. The fitted parabola puts the `nmae` optimum at **α ≈ 1.29**. `reach`
and `fid` both moved the same way across the v11→v12 change, and `pds` is
exactly scale-invariant (exp 001, bit-identical across α 0.01→100).

The leader supplies the independent confirmation. Jurassic Park's jump from
0.12595 to 0.16588 was a **≈2.07× amplitude increase and nothing else** —
`pds_cosine` 0.786473 → 0.786700, unchanged to three decimals. Their raw `mse`
went 1.0733 → 3.6867, which cost them nothing because it was already clamped,
and bought `reach` ×1.79 (**+0.0212 overall, 53% of the jump**), `nmae` +0.0098,
`fid` +0.0084. They had earlier tried the `mse` route, landed at 1.0733, missed
the credit threshold by 0.068, banked zero, and reversed.

**Against this stands exp088's h1-test sweep**, which measured α=1.2 at −0.0037
and α=2.0 at −0.0610 with "nmae collapsing each time". That is the harness
family that was wrong by 4× on `nmae` in exp031 and by 0.029 raw on `mse` in
exp077. The board is the instrument; the proxy is the thing being calibrated.
Resolved by measurement, not argument — see v13.

## 4. Where the remaining score actually is: the call set

`fid = (k/n_pred)·min(1, n_pred/n_real)`. The coverage term saturates at
`n_pred = n_real`, so the function **peaks exactly there**, and it is sharply
asymmetric: at 60% sign accuracy, overshooting by 20% costs `0.009` overall
while undershooting by 20% costs `0.066`. Aim high, never low.

`jac = |R ∩ R̂| / |R ∪ R̂|` charges every false positive forever, with no
saturation — precisely where `fid`'s coverage term stops charging and `reach`
never charges at all.

Both want the same thing, and neither is a modelling question. From
`donor_transfer.py:261`:

> The reason a minimum-variance emitter calls essentially every gene is not the
> effect sizes — it is that 400 near-identical values sit at the control
> distribution's *mean*, which for skewed counts is above its *median*, so the
> rank-sum against 18,400 real controls is significant even at zero effect.

And from the same file, on `mode="even"`: *"Over-calling costs `jac`, which
divides by the union — the trade is deliberate and must be checked against the
scorer, not assumed."* **It was never checked.**

`mode="quasi"` holds each gene's panel total bit-exact, so this channel cannot
move `pds`, `mse` or `nmae` by construction. `emit_counts` already accepts
`spread` as a **per-gene vector**. `make_submission.py:260` passes a scalar.
The control channel is built, wired, and switched off.

Modelled at 60% sign accuracy with `n_pred ≈ n_real`: **+0.054 on `fid`,
+0.03–0.07 on `jac`.** That is an order of magnitude more than anything else
currently on the table, and it requires no change to the predicted mean profile.

**The blocking measurement is `n_real`, `n_pred` and `k` per perturbation.**
They are reference-only or prediction-only quantities, computable offline on the
HCT116 harness. Until they exist, §4 is a hypothesis.

## 5. What this retires

- "Raise the transfer cosine to open `mse`" (Addendum 4's open direction) —
  `mse` is clamped and free; there is nothing to open.
- The implied-cosine board table (`STRATEGY.md` §Addendum 4) — it inverts board
  `mse` through a normalization the retraction already showed is wrong, and the
  measured floor is 0.0384 away from what exp077 reported.
- Per-perturbation magnitude calibration (exp012) — its oracle was worth 0.25%
  squared-error reduction, and that measurement stands. It is not the lever.
- Lineage-matched donors — `nadig_jurkat`, `rpe1` and `k562_essential` all cover
  **0/300** of the panel. Closed on data availability, not on principle.

---

# Addendum 6, 2026-08-27 — the emission knobs are done; the instrument was lying

Addendum 5 got two things right and one thing badly wrong. Right: the scoring function is
linear and solvable, and `pds` is finished. Wrong: **"alpha's marginal direction is up".**

## 1. Amplitude is retired, on board evidence

Fitting `scaled = (raw − b)/(r − b)` on **780 board entries** (union of all snapshots,
363 teams, 2026-08-20 → 2026-08-27) reproduces four members at `R² > 0.9999`:

| member | `b` | `r` | span | raw Δ for +0.01 overall |
| --- | ---: | ---: | ---: | ---: |
| `pds` | 0.50020 | 0.95228 | 0.452 | +0.0271 |
| `fid` | 0.51239 | 0.80964 | **0.297** | **+0.0178** |
| `reach` | 0.07902 | 0.96414 | 0.885 | +0.0531 |
| `jac` | 0.03038 | 0.39727 | 0.367 | +0.0220 |

`score_mse` and `score_nmae` are **not** linear in the reported aggregate raw values —
they are per-perturbation normalized then averaged. Do not do arithmetic on them.

Then, on the same union:

* Between teams, `log(raw mse)` explains **9.7%** of `reach` variance. Our reach residual
  against that fit is **−0.0056 (z = −0.15)**. We are on the line; there is +0.001 there.
* Within teams, over **288 consecutive submission pairs**, `Δscore_avg` is **uncorrelated
  with `Δlog mse`** (r = −0.07, **p = 0.24**).
* `Δpds_cosine` does not move with amplitude at all (R² = 0.003, p = 0.32).
* **We527** ships `STATE baseline` at raw mse `1.3577` — within 2.6% of our `1.3236` — and
  scores reach **0.2529** against our **0.1815**.
* **Aginglab** went `GeroAI_v9 → v11` by *halving* amplitude (`1.822 → 0.959`) and *gaining*
  **+0.100** reach, with `nmae` `0.9216 → 0.8923` (board best) and `pds` slightly down.

Addendum 5 §3's parabola through three `nmae(α)` points, and its reading of Jurassic
Park's jump as "amplitude and nothing else", are both superseded. **v13 (α=2.0) is built
and must stay unsubmitted.**

## 2. Every scalar knob is at or adjacent to its measured optimum

| knob | verdict | evidence |
| --- | --- | --- |
| `alpha` | 1.0796 stands | exp088: 1.2 → −0.0037, 2.0 → −0.0610; board: Δscore ⊥ Δamplitude |
| `min_cpm` | **sharp peak at 10** | exp088b: 8 → +0.0028, 10 → +0.00455, 15 → −0.0056, 30 → −0.0129 |
| global `spread` | 0.25 stands | exp034 plateau 0.20–0.30; exp088 h1-test 0.15 ≈ 0.25 |
| confidence-weighted spread (`tight_k`) | +0.00075, mostly `jac` | exp088 Track A on h1-test |
| call-set truncation | **negative on the board** | v11 `0.150183` > v15 `0.149436` > v14 `0.147273` |
| `basal_beta` | wrong direction for us | exp088 §4: buys `pds`, which we already lead |
| `mse` via low `alpha` or low `gamma` | **closed** | exp088 §1: −0.387 `pds` to buy +0.034 `mse` |

## 3. Why the proxies kept lying: the harness is depth-confounded

**`data/harness/orion-hct116.h5ad` has control cells 38% deeper than its perturbation
cells, in 132/132 perturbations**, because `build_orion_harness.py` fills the control pool
greedily from the deepest early batches. The Wilcoxon is a rank test and CPM does not fix
the per-cell zero fraction, so the reference's DE calls skew **90.8% negative** against
`0.538` on h1-test, `0.497` on Replogle K562 and `0.553` on RPE1. Thinning proves it:
invert the depth mismatch and the skew inverts to `0.275`; equalize it and it lands at
`0.541`.

**This is why experiment 093 measured +0.0202 and the board returned −0.00291** with `fid`
at the wrong sign. Its own BOARD-RESULT blamed cells-per-perturbation; the real cause is
that the harness's `n_real`, its precision curve and its direction distribution are all
wrong together. Everything in `experiments/` calibrated on that file needs re-reading:
031, 061, 065, 090, 091, 092, 093, and the HCT116 half of 088. Full diagnosis:
`docs/HARNESS-DEPTH-CONFOUND.md`.

**Use `data/harness/h1-test.h5ad`.** exp088 §5 established it as reach-representative
(raw reach 0.23–0.24 against the board's 0.169; HCT116 is 0.015) and its depth ratio is
0.985.

## 4. What the real challenge artefact looks like

Auditing `submissions/v11.vcc` against the real `data/challenge/context_*.h5ad` controls —
no ground truth required:

* **`n_pred` ≈ 6,000** per perturbation, **61% of the entire gated universe**.
* **75–81% of our significant calls carry `|lfc| < 0.10`**, and ~35% are internally
  contradictory: the rank-sum earns its p-value up while the reported lfc says down.
* Our emission is **under-dispersed ~2.3×** — per-gene CV 0.38–0.40 against 0.84–0.92 in
  the real controls.
* In **linear** CPM space **100%** of gated genes have mean > median and the control's own
  mean sits at the **58.3rd percentile**, so emitting 400 cells at the predicted mean is a
  free 8.3-percentile "up" at zero effect.
* `Spearman(p_adj rank, |lfc| rank) = +0.015`. `--keep-top` ranks on `|effect|`, so it
  never had any purchase on the order the scorer actually walks.

## 5. Where the remaining score is

| member | headroom to board max, in overall points |
| --- | ---: |
| `mse` | +0.0268 — but costs `0.387 pds` to buy `0.034`. Closed. |
| **`reach`** | **+0.0162** — pure confidence-ordering quality |
| **`nmae`** | **+0.0142** — magnitude accuracy on genes that really respond |
| `fid` | +0.0116 — we are *below* the 0.51239 no-information baseline |
| `jac` | +0.0024 |
| `pds` | **0.0000 — we are the board maximum** |

`reach` is ~80–90% unexplained by every other reported metric, so it is a nearly orthogonal
axis. The board's strongest pattern: **every entry naming STATE lands at reach 0.21–0.27
across raw mse spanning 1.36 to 1,774** — amplitude-invariant and model-driven.

**The honest conclusion: our `reach` is about what our per-gene directional confidence
supports, and there is no free re-ordering trick.** The next real gain is a better
directional signal or a better operator — not an emission parameter. That makes the
operator the live surface, which is why experiment 091's pruned basis ships today.

---

# Addendum 7, 2026-08-27 — the architecture the metrics actually reward

Addendum 6 closed every scalar knob. This is what is left, and it follows from a structural
fact about the submission format that took this project six days to state plainly.

## The submission has two independent channels

We emit integer counts. The scorer derives all six metrics from them, but through two
disjoint paths:

| channel | set by | metrics it moves | metrics it cannot move |
| --- | --- | --- | --- |
| **the profile** — each gene's panel total | the model's predicted effect | `pds`, `mse`, `nmae` | — |
| **the emission** — how those totals are scattered over 400 cells | the per-gene `spread` vector | `fid`, `jac`, `reach` | `pds`, `mse`, `nmae` |

The second row is not a hypothesis. `emit_counts(mode="quasi")` repairs each gene's panel
total **exactly**, so the pseudobulk is bit-identical whatever `spread` does — and the board
confirmed it in v14, where `pds` moved by **exactly 0.000000** and `nmae` by 1e-5 across a
full-range spread change.

**So the confidence ordering that `reach`, `fid` and `jac` are computed from is a free
parameter, orthogonal to the predicted expression profile.** We have been setting it by
accident: whatever p-value ordering falls out of a minimum-variance emitter, which
correlates with our own model's magnitudes at `Spearman +0.07`.

## What that is worth

| | overall points |
| --- | ---: |
| our `reach` today | −0.0004 scaled, ≈ board rank 94/363 |
| amplitude-matched expectation | +0.001 |
| a crude `\|effect\|`-keyed promotion (exp095, K=2000) | **+0.0058** |
| matching Aginglab's reach | **+0.0162** |
| oracle ordering (correct signs first) | **+0.103** |

And `fid`, whose raw value of `0.4967` sits **below the 0.51239 no-information baseline**,
is the cheapest member on the board: `+0.0178` raw buys `+0.01` overall.

## The architecture this implies

We hold **the best `pds` of 363 teams** and a bottom-third `reach`. Those two facts are not
in tension — they are measured through the two different channels above. So the shape of the
next real gain is:

> **keep our profile, and drive the emission from a better per-gene directional confidence.**

The board says such a signal exists and is not ours. **Every entry naming STATE lands at
reach 0.21–0.27 across raw mse spanning 1.36 to 1,774** — amplitude-invariant, model-driven,
and roughly where Aginglab and Evidoa also sit. Arc's own State model is explicitly allowed
and explicitly "not able to win", which makes it usable as a *component* rather than a
submission: its per-gene ranking, our profile.

Candidate confidence keys, cheapest first:

1. **Cross-donor sign agreement per gene.** exp013 measured `√(|A|·|B|) × sign-agreement` at
   +24% reach over magnitude ordering on a LODO proxy; exp088 Track A got only +0.00075 from
   it on h1-test, but that arm was confounded with a spread change and had no `k=0` control.
   Cheap, already implemented in `experiments/031-somatic/run.py` as `tight_k`.
2. **A learned per-gene direction classifier** trained on public Perturb-seq to predict
   `sign(lfc)` given (target, gene, donor deltas, basal expression). Rules-legal: training on
   public data is explicitly allowed. Its *calibrated probability* is exactly the ordering
   key `reach` wants, and it does not touch the profile.
3. **STATE-family rankings** as a component, per above.

## The constraint that kills naive versions

`fid = precision × min(1, n_pred/n_real)` and `jac` charges every false positive, while
`reach` only cares about order. Any policy that improves the ordering by *shrinking* the
call set pays `fid`'s coverage term — that is exactly how v14 and v15 lost. The channel to
use is the one that **re-orders without shrinking**: lower `spread` on the genes we trust
rather than raising it on the ones we do not. That is what `--keep-spread` was added for.

## What to measure next, and on what

Not on `data/harness/orion-hct116.h5ad` — see `docs/HARNESS-DEPTH-CONFOUND.md`. Use
`h1-test` for absolute DE numbers, `orion-hct116-thinned.h5ad` for controlled A/Bs, and gate
any harness with `scripts/check_harness_depth.py` first. The board remains the only
instrument that has never lied, and it is noiseless for a fixed artefact — v11 and a
byte-identical rebuild score identically, so a single submission measures a change exactly.

---

# Addendum 8, 2026-08-27 (later) — Addendum 7 was half right, and the half that was wrong matters

Addendum 7, written five hours earlier the same day, said the confidence ordering the DE
metrics are computed from is "a free parameter, orthogonal to the predicted expression
profile", and pointed at driving it from a better per-gene signal. Experiment 095 then
measured that claim properly and it splits in two.

**Still true: the ordering is free.** Moving it cannot touch `pds`, `mse` or `nmae` — the
quasi emitter's exact-total repair guarantees it and the board verified it in v14, where
`pds` moved by exactly `0.000000`.

**False: that the ordering is *steerable*.** With the emitter we ship, `spread` is a **weak
lever on order and a strong lever on membership**. Nine policies, realised through the
unmodified production `emit_counts` and scored on the exact DE path with a **5e-5**
emission-noise floor:

| policy | Δoverall (realised) | `n_pred` |
| --- | ---: | ---: |
| block K=500 … K=6000 | −0.0001 … −0.0007 | 4,023 … 4,922 |
| graded `\|effect\|`, 0.00 → 0.50 | −0.00077 | 3,918 |
| graded **agreement**, 0.00 → 0.50 | −0.00193 | 3,911 |
| block K=2000, tail 0.60 | **+0.00042** | **2,499** |
| graded 0.00 → 1.00 | +0.00024 | 2,594 |

Every policy that leaves the call set near production is at or below zero. The three that
are positive all shrink it by 15–40%, and the board has priced that trade twice: v14 and v15
lost `0.0029` and `0.0007` doing exactly that.

**And the counterfactuals were badly wrong.** The K=2000 block promotion measured `+0.00583`
when the DE table was reordered by fiat and **−0.00021** when actually emitted. The reason is
mechanical and general: predicted p-values span hundreds of orders of magnitude, `spread`
moves them by a few, so a promoted gene that was not already near the front does not travel.

## The methodological lesson, which is the durable part

This project has now been burned three times in the same shape:

| | proxy said | board said |
| --- | ---: | ---: |
| exp093 truncation (counterfactual, then re-emitted on a confounded harness) | +0.0202 | −0.0029 |
| exp091 basis prune (same-harness A/B, confounded) | +0.00088 | +0.00013 |
| exp095 ordering (counterfactual on two harnesses) | +0.00583 | *not spent* — caught first |

The third was caught because the counterfactual was forced through the real emitter before a
slot was spent on it. **That check — `experiments/095-clean-harness/realize.py`, with a
same-policy second seed as its noise control — is now the gate any emission-channel claim
must pass.** A reordering that cannot be produced by a spread vector is not a result, it is
a wish.

## Where that leaves the search

Closed by measurement, not by discouragement: amplitude, `min_cpm`, global spread, call-set
truncation, `basal_beta`, `mse` via low alpha or low gamma, basis-prune tuning
(λ₂ falls 56% and buys `+0.000126`), and now the whole emission ordering channel.

What is left is what Addendum 7 got right for the wrong reason. `reach`'s oracle is `+0.103`
and cross-donor agreement orders the pool at `+0.00678` **as an ordering** — the signal is
real, we simply cannot express it through `spread`. Expressing it needs either

* a change to the emitted **distribution shape** — the gamma-Poisson extension exp093 §7
  probed and never validated, which moves p-values by orders of magnitude rather than a few;
  or
* a better **profile**, which is not free, but which is also where `nmae` (`+0.0142`) and
  `fid` (`+0.0116`) live, and where every team ahead of us is winning.

We hold the board's best `pds` and are bottom-third on `reach`. That gap is a model gap, not
an emission gap, and today established it is not reachable with a knob.

---

# Addendum 9, 2026-08-30 — STATE is closed; new measured biology is the queue

The fixed step-20,000 `state_sm` reproduction is not a usable perturbation
predictor. On held HCT116/HEK293T, pure STATE PDS was `0.50252 / 0.50028`
against the incumbent's `0.80917 / 0.79454`; a 10% direction blend lost
`-0.02574` macro PDS. Experiment 113 then tested the same frozen checkpoint in
its H1 training-domain context: four independent control banks all remained
near random (`0.4940–0.5133`), pooled PDS was `0.50448`, and the incumbent was
`0.86345`. Higher mean cosine did not encode perturbation identity.

Experiment 110 also stops before two additional trainings. Its
perfect-information HEK293T disagreement-routing ceiling is only `+0.00686`,
below the frozen `+0.015` gate, and its first scaffold relied on a
reference-significant mask unavailable at deployment.

The live strategy is therefore:

1. finish experiment 107's manifest-bound Jurkat acquisition and integrity
   gates;
2. fit only the frozen, identifiable target-level/channel-held-out model;
3. require independent Nadig Jurkat target-generalization plus clean
   held-context gains before any production build;
4. if 107 fails, advance experiment 109's gauge-invariant guide-nuisance
   runner;
5. do not submit a15 or any new arm without a completed, preregistered
   promotion gate.

Experiment 107 is not an ordinary fourth donor. Its median 13-guide high-MOI
design estimates marginal effects in a combinatorial background, and its
reliable response axis is only about 349 PDS-visible coordinates after target
masking. The first integration must therefore be a masked post-whitening
residual with bit-exact fallback, not a full-vote donor average. Context-A-only
routing remains forbidden until independently validated; basal identity alone
has already failed as a selector.

---

# Addendum 10, 2026-08-30 — panel corrections are bounded; exploit the new H1 screen

Experiments 122, 125, 127, and 128 isolate a real but limited panel-only
channel. The exact integer-budget transport is mechanically sound, but the
oracle ceiling shows why it cannot be the main strategy: only about 2.28% of
reference DE membership and absolute LFC mass lies on mutable VCC-panel
coordinates. Even perfect panel correction is worth only `0.024665` raw NMAE
on the equal-surface development average, and its HCT116 Jaccard ceiling is
`0.000579` raw.

The search priority therefore returns to genome-wide measured biology, but
with a corrected evidence boundary. Experiment 115's GSE295214 H1 screen is
broad—2,354 post-QC perturbations and 16,290 covered canonical genes—but it is
not untouched. Its direct residual was already tested in Experiment 117 and
its 64 challenge rows supplied Experiment 118's Pan confirmation. Those
residual, mask, dose, and confidence choices are closed. The source may still
support a materially different learner trained and evaluated only on its more
than 1,300 non-panel overlaps with RPE1 and Nadig Jurkat; Experiment 129 is
restricted to that role and must stop if it cannot define a genuinely new
mechanism.

No public source currently combines at least 40 exact challenge targets, broad
response coverage, and complete independence from consumed outcomes. If
another source-only panel candidate passes, the best fresh panel-response
confirmation is IGVF TeloHAEC: 65 targets and 17,464 response symbols. Its
broad direct-donor route failed Experiment 105, but the panel coordinates were
not opened. HCASMC is a more orthogonal secondary confirmation with only 27
targets. Both require a new frozen panel-only protocol before any values are
accessed.

Experiment 125 does not supply that candidate. Its posterior confidence
improved top-20% precision by `+0.01217`, but failed both null thresholds, and
its posterior mean reduced broad signed credit by `-0.01505` with all three
families negative. This is a clean NO-GO: the transport was exact and support
was unchanged. Do not spend either confirmation source on it.

The DepMap/STRING weak-row replacement is closed: its frozen primary arm lost
`-0.120168` raw PDS on the authorized HEK293T diagnostic, and both component
arms were also strongly negative. STATE, guide nuisance, direct added-donor
voting, STRING, regulons, and post-hoc wrong-sign flipping remain closed by
their existing experiments.

---

# Addendum 11, 2026-08-30 — preserve the incumbent; data and constrained residuals are the frontier

Experiments 129 and 132 close two tempting ways to use auxiliary biology.
Pan-neighborhood transfer improved directional credit but reduced PDS by
`-0.0110` to `-0.0209`; DepMap PCA50/Ridge10 replacement lost `-0.204710`
macro PDS. Both failures have the same structure: a generic biological prior
contains real signal but discards or dilutes the much stronger same-target
response already present in the incumbent.

That establishes the next design rule:

1. never replace the same-target source mean with target features, neighbors,
   or a generic learned response;
2. express new information only as a low-capacity residual around the exact
   incumbent;
3. require nonnegative complete-gallery PDS in every outer context, not merely
   better cosine or signed credit;
4. keep identity as the deterministic fallback and leave a fresh whole-target
   test sealed through selection;
5. do not spend a scarce held biological surface on a source-only family that
   cannot clear those gates.

The differentiated-cell source audit clarifies the data bottleneck. HUES8 ESC,
HUES8 endoderm, WTC11 hepatocyte, and WTC11 cardiomyocyte have useful target
coverage, but the broad canonical signed outputs are access-restricted.
Target-level energy distances are not substitutes for gene-level directions,
and the small public cardiomyocyte table is the same sparse family already
rejected by Experiment 081.

Experiment 080 is therefore the only active heavy model test. Its listwise
objective directly optimizes retrieval while retaining an identity fallback,
which makes it better aligned with the competition than the recent credit-
first learners. A margin-constrained directional residual is the most
promising follow-up in principle, but it is not an experiment until it has a
fresh target split and an independently authorized proposal vector. A basal-
conditioned residual is likewise deferred until a fresh context has a valid
control pseudobulk.

---

# Addendum 12, 2026-08-30 — the fixed panel atlas fails independent confirmation

Experiment 135 tested the exact Experiment 122 `equal_mean_agreement` panel
atlas once on 65 previously unaggregated TeloHAEC panel targets. The evaluator
used the exact staged-a15 comparator and the mechanically validated Experiment
128 integer-budget transport.

The atlas carried a small broad signed signal: signed credit improved by
`+0.004066` and exceeded both permutation-null p95 values. That was not enough
to pass. The gain missed the frozen `+0.005` minimum, its bootstrap lower bound
was negative, top-20% sign precision worsened by `-0.005245`, and median
integer-transport sign retention was only `0.890`.

This closes the atlas, not merely its nominal dose. Do not use the consumed
TeloHAEC result to rescale, subset, rerank, or choose targets, and do not spend
a submission on an a15-plus-panel arm. The result also sharpens the strategic
conclusion from Experiment 127: panel-only corrections are too small and too
unstable to be the main route to first place.

Experiment 080 remains the sole active high-value test because it directly
optimizes complete-gallery retrieval while preserving the incumbent fallback.
After it reports, any follow-up must be materially distinct, preserve the
same-target incumbent, use a fresh source-target split, and clear
complete-gallery PDS in every outer context before consuming another held
biological surface.

Experiment 136 rules out the obvious DepMap Ridge version of that follow-up.
Every nonzero residual dose traded PDS away for signed credit, source-rank
preservation did not transfer, and exact source-margin preservation collapsed
to identity. Even a per-target held-outcome oracle stayed below the `+0.003`
PDS gate. A margin-constrained successor therefore requires a genuinely
different residual proposal; a more elaborate router cannot rescue this one.

---

# Addendum 13, 2026-08-30 — listwise learning does not transfer; use a measured residual

Experiment 080 passed every source-only pre-gate with a nonidentity rank-8
adapter and a large `+0.023376` untouched non-panel PDS gain. On the final
held panel, that gain vanished: HCT116 was `-0.0006884`, HEK293T was
`+0.0006435`, and the primary macro was `-0.00002245`. KOLF was positive but
small at `+0.0006868`. The bootstrap crossed zero and the observed macro
failed the permutation threshold. This is a terminal NO-GO.

The result closes a broad learned global geometry correction, even when its
objective is aligned exactly to complete-gallery retrieval. The failure is
not excessive angular movement—the median displacement stayed below one
degree in every final fold—but lack of transferable target-specific benefit.

The next candidate should therefore use a measured, target-specific proposal
whose prior evidence is already close to the gate. Experiment 107's frozen
high-MOI Jurkat source model is the strongest available option: it previously
gained `+0.004246` complete-gallery PDS and passed bootstrap/permutation
evidence, missing promotion only on one slightly negative gallery and
directional credit. A new experiment may test a tangent-plane/spherical
projection of that fixed proposal around the exact incumbent, with identity
fallback and a fresh 79-target source-only split. It must be preregistered and
independently reviewed before any effect values from that split are opened.
