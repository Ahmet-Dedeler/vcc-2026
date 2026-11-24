# Scientific review — September 9, 2026

Scope: continuation of Ahmet's interrupted Claude review, performed by the main
agent only. No subagents, new training, reserve access, cloud compute rental,
submission, or manuscript publication. Board arithmetic and core-input hashes
were recomputed; historical experiment results were reviewed, not rerun.

## 1. Where we actually are

At **04:27:12 UTC on September 9**, our published score was
**0.21865973944346428**, **5th of 821**. This is the best in the 18 distinct
entries recovered from the retained board snapshots, and agrees with the
recorded best in the experiment history. It is not a new score regression.
`TOBE` has overtaken us. The leader is **0.2620924638410216**, not 0.26.21;
its latest entry was submitted **September 8 at 19:53:54 UTC**.

Source: `data/leaderboard/2026-09-09T0427Z.json`; independent authenticated
`vcc --json status xuDelEik7zGWMXdNu0R4`; computed evidence at
`docs/review-2026-09-09/state.json`.

The large earlier climb had identifiable contributors. Relative to a15
(0.1740865), these successive improvements brought us to cd4svd5 (0.2160944):

| Change | Observed overall gain |
|---|---:|
| Rank-controlled emission, rc_hybrid | +0.0093325 |
| Call-set change, rc_callset | +0.0084448 |
| CD4 dose to f=0.34 | +0.0130557 |
| H1 union residual | +0.0069114 |
| CD4 deflation | +0.0042634 |
| Subsequent best normrel/cap construction, vs cd4svd5 | +0.0025653 |

Source: `docs/BOARD-LEDGER.md`, the captures it cites, and 319's new
`BOARD-RESULT.md`. These are successive changes, not independent effects that
can be recombined arbitrarily. **“Best score ever” was correct; “best arm ever
was +0.0016” was not.** The largest individual gain in this sequence was about
eight times that. Small recent gains do not mathematically bound future model
improvements.

## 2. The current gap is not “entirely cosine”

Compute `(leader.score_member − ours.score_member) / 6` directly:

| Member | Contribution to leader's +0.0434327 advantage |
|---|---:|
| PDS | +0.0050646 |
| MSE | **−0.0108983** |
| NMAE | **+0.0231882** |
| fidelity | +0.0088578 |
| reach | **+0.0189888** |
| Jaccard | **−0.0017685** |

We are better on MSE and Jaccard. The leader's NMAE and reach advantages
together are +0.0421770 overall, offset in part by its MSE disadvantage.
This is a **score decomposition**, not proof of which hidden modeling choice
caused it. Our member ranks are PDS 4, MSE 21, NMAE 92, fidelity 129, reach 12,
Jaccard 1. That supports examining downstream response/sign accuracy without
throwing away the call-set machinery that already works.

The leader's change since the September 8 12:53Z snapshot is **+0.0219077**,
with positive contributions on all six members. The largest are MSE
(+0.007184), fidelity (+0.006580), and NMAE (+0.005217). That trajectory does
not identify PromoterAI, chromatin data, or any architecture as its mechanism.
A public model name is a clue, not a methods disclosure.

Two invalid shortcuts must stop:

1. **MSE is not a universal cosine meter.** Its numerator depends on prediction
   norms and jackknife credit as well as alignment. A fitted inversion within
   our emission family cannot identify a rival's cosine. Requiring positive
   scaled MSE is necessary to avoid one clamp error, but matching a raw-MSE
   range is not sufficient to identify the other terms.
2. **PDS is not own-target cosine either.** For a fixed gallery with truths
   `(1,0)` and `(0.8,0.6)`, a prediction with own cosine 0.9 can be closer to the
   decoy and score PDS 0; another with own cosine 0.5 can score PDS 1.
   The runnable counterexample is in `scripts/audit_takeover_state.py`.
   Thus the old READ-FIRST claim that better PDS makes substantially worse
   own-target cosine impossible was mathematically false. Whitening already
   demonstrated this distinction empirically in this repo.

Primary definition: `docs/reference/vcc2026-metrics-BRIEF.md` §§1–2.
The previous READ-FIRST is preserved, not silently rewritten as if it had
always said this.

## 3. What recent experiments establish, and what they do not

**312/316/319: a useful engineering improvement, not an identified noise
mechanism.** Cap beats plain by +0.0009409 at exponent −1.5. But cap changes
which rows move, donor shares, effective dose, delivered directions and
amplitudes. It moves weights about half as far at the median. The score
contrast does not isolate “small donor rows are noisy.”

316's norm match is exact **before whitening**. An anisotropic operator does
not preserve that identity; 318 directly measured delivered redistribution
in both arms. The pure mathematical example is `diag(2,1)` applied to two
unit-length rows: output norms are 2 and 1. A fixture regression proves its
stated invariant at its stated stage, not after the rest of the pipeline.

**319d3 is informative but not urgent.** Its already-frozen 2×2 comparison
(plain/cap × −1.5/−3.0) can test whether the dose response differs. It cannot
by itself identify noise suppression, and 318 did not rule out all magnitude
explanations for it. Keep the existing score bins unchanged. No automatic
midnight scheduling; decide when the result will change the next research
decision. See 319's September 9 coordinator review.

**317/318 provide local negatives, not a theorem closing amplitude or DE
prediction.** Their known-answer baseline is the historical `rc_callset` H1
prediction (`048b9286…`), not a current 319 A/B/C-equivalent model. Their
scorer reproduction is valuable; transfer to a different base and domain is
not supplied by that reproduction.

318's rank-mapping strengthens the correlation with inverse predicted norm
from about −0.30 to −1. It is not a mathematical upper bound on useful
redistribution: imposing a stronger rank relation can make a wrong rule worse.
Its three shuffles yield a descriptive spread, not a reliable significance
calibration; formal permutation resolution with three randomizations is at
best 1/4 under the usual plus-one rule. Its emitted whole-matrix norms also
differ, and its uniform-effect subtraction is a local approximation. Keep
the tested rule NO-GO; do not claim every per-target magnitude model is dead.

Likewise, stable aggregate PDS/MSE is **consistent with**, but does not prove,
unchanged predicted profiles. Two aggregate numbers cannot establish identity
of a high-dimensional prediction. For our own packages, use actual count or
effect digests; for rivals, label it an inference.

**322: do not rediscover KOLF breadth, or misstate it.** The saved bank has
273 panel rows, but the upstream atlas contains 10,167 perturbations; 322
actually re-extracted it and still found no useful PDS gain in the tested
constructions. The narrower bank was an extraction choice. Two synthetic
panel draws and particular source constructions do not prove that all future
breadth/representation changes are useless, nor that all rivals must lose the
same amount in October.

**F16: retain the measurement, withdraw the discovery framing.** 320/321
already did the literature and wording correction. Its very small headline
p resamples genes; it is not independent construct replication or an estimate
of the 2026 panel's failure rate. Repeating known specificity biology in a
new sample is useful, but does not establish a new mechanism. No new wet-lab
spend follows from this review.

Evidence: experiments 312/313/316/317/318/319/320/321/322 `RESULTS.md`;
`docs/FINDINGS-2026-09.md`; `experiments/294-writeup/CLAIMS.md`. I did not
independently rerun their biological analyses or claim a fresh exhaustive
literature search.

## 4. The next scientific work I would prioritize

These are **proposals, not launched experiments or claimed discoveries**.

**Continuation update:** the first proposal was implemented as experiment
323 on September 9. Construction and all-six-metric regression gates pass,
but none of four small historical contrasts has a resolving paired interval.
It is not a standalone selector. See
`experiments/323-current-family-validation/RESULTS.md`; the remaining proposals
below retain their stated feasibility and authority gates.

### First: make a small, current-family, decision-capable validation panel

The bottleneck is not a shortage of possible knobs. We keep asking an old,
different-domain baseline to decide whether a new construction has value.
Start with already-used data, without opening any reserve, and measure:

- Which held study/context/construct contrasts have repeatable downstream
  response at the cell/depth budget being evaluated?
- Which are matched to the challenge's target/effect-strength regime, and
  where is the mismatch unavoidable?
- Can the exact current-family composition be replayed without the held
  context/targets contaminating the donor fit or whitening basis?
- Does it reproduce known positives **and negatives** in their stated scopes?

Keep H1 union **off** on H1 validation. Do not pool target–gene pairs into
an easier statistic, and do not claim a perfect proxy if the necessary
held-safe data does not exist. The deliverable is an eligibility/coverage
table and held-safe regression panel, not a new model. If no available
validation design resolves the intended effect, report that before paying
for a broad model sweep.

This builds on 267/273/278/285 and the repairs already made; it does **not**
repeat their searches for more cells or declare them insufficient by fiat.
The difference is the current composition, all-six-member decision, and
explicit construction-level transfer test.

### Second: test a response predictor, not another scalar confidence knob

A plausible larger improvement has to change **which downstream genes are
predicted to respond and in which direction**, conditional on target and
recipient, while preserving the emission/call-set machinery that already wins.
A modest shared-response model plus a context-conditioned residual head is a
candidate class; it is not a new idea or a promised breakthrough.

Before implementing it, map its exact differences from the failed per-gene
gain, basal ridge, utility router and context-residual constructions
(132/161/224/268/269/283). No rerun under a new name. Demand a testable new
ingredient, such as a reproducible cross-study residual label rather than a
noisy selected donor oracle. Conservation predicting *difficulty* alone is
not enough: 277 already showed that such a ranking need not improve outputs.

Freeze one small model and controls: no-recipient head, shuffled recipient,
shuffled target, capacity-matched shared model, and a repeatability-matched
label control. Hold out the whole recipient context/study; group related
targets where relevant; fit preprocessing, operators and selection only on
training data. Report exact scorer metrics plus downstream sign calibration
and residual reproducibility, with target-level uncertainty and per-context
effects. Any model that merely raises cosine while sacrificing PDS or the DE
members must be priced as that trade, not called progress.

**Stop condition:** no held-context advantage over the shared model, or gains
matched by target/recipient shuffles, stops the candidate before production
packaging. A positive on reused development data earns a separately frozen
confirmation design, not a discovery claim.

### Third: pursue the reproducible unexplained component as the science track

284's residual conservation signal across disjoint sources/lines is a stronger
starting point than another generic off-target story: the residue was
reproducible even when annotation features did not predict it. Its phenotype
is an observation to explain, not proof of an unknown mechanism.

The next question is whether that residual survives a crossed
study/line/construct analysis and predicts an independent response property,
rather than a pipeline-specific agreement score. Audit available independent
units from metadata first. Declare the novelty comparison and the smallest
effect of interest before reading new values. Because the prior reserves
were consumed, further analysis of them is exploratory; a confirmatory claim
requires genuinely new independent evidence.

**Stop condition:** failure to separate the effect from shared measurement,
selection or batch structure means no mechanistic claim. No assay or cloud
budget should be spent just to make an already-reported result look more new.

## 5. Process changes that matter more than more agents

- Use one coordinator and serial work for this task, as requested.
- The current-state page contains current facts and pointers; long
  histories stay archived. A result has one authoritative record, not five
  contradictory “current” summaries.
- Every number is labelled by **population, stage, statistic, units,
  reference, and validation status**. This prevents the fixture/production,
  rank/cosine, bank/source and gene-/construct-sampling errors above.
- An experiment should name the decision it can change. Daily allowance
  neither forces a submission nor justifies a predetermined queue.
- Separate engineering readiness, scientific readiness, novelty, and
  authority. A package can pass integrity checks and still lack a compelling
  next scientific use.
- A score of 1 is replicate-normalized performance, not “100% solved” or
  the maximum possible score (metric brief §8). Near-perfect prediction is a
  legitimate aspiration, but this record does not establish that it is
  reachable from the supplied controls and public donors.

The immediate recommendation is to finish the current-family validation
design and its identifiability checks, then choose one response-model
hypothesis. Preserve 319 as the best known artifact; retain 319d3 as an
optional bounded information experiment, not a timer-driven obligation.
